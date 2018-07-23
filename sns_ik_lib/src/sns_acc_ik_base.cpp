/** @file sns_acc_ik_base.cpp
 *
 * @brief The file provides the basic implementation of the SNS-IK acceleration solver
 *
 * @author Matthew Kelly
 * @author Andy Park
 *
 *    Copyright 2018 Rethink Robotics
 *
 *    Licensed under the Apache License, Version 2.0 (the "License");
 *    you may not use this file except in compliance with the License.
 *    You may obtain a copy of the License at
 *    http://www.apache.org/licenses/LICENSE-2.0
 *
 *    Unless required by applicable law or agreed to in writing, software
 *    distributed under the License is distributed on an "AS IS" BASIS,
 *    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *    See the License for the specific language governing permissions and
 *    limitations under the License.
 */
#include <sns_ik/sns_acc_ik_base.hpp>

#include <ros/console.h>
#include <limits>

namespace sns_ik {

/*************************************************************************************************
 *                                 Public Methods                                                *
 *************************************************************************************************/

SnsAccIkBase::uPtr SnsAccIkBase::create(int nJnt)
{
  if (nJnt <= 0) {
    ROS_ERROR("Bad Input: ddqLow.size(%d) > 0 is required!", nJnt);
    return nullptr;
  }
  Eigen::ArrayXd ddqLow = NEG_INF*Eigen::ArrayXd::Ones(nJnt);
  Eigen::ArrayXd ddqUpp = POS_INF*Eigen::ArrayXd::Ones(nJnt);
  return create(ddqLow, ddqUpp);
}

/*************************************************************************************************/

SnsAccIkBase::uPtr SnsAccIkBase::create(const Eigen::ArrayXd& ddqLow, const Eigen::ArrayXd& ddqUpp)
{
  // Input validation
  int nJnt = ddqLow.size();
  if (nJnt <= 0) {
    ROS_ERROR("Bad Input: ddqLow.size(%d) > 0 is required!", nJnt);
    return nullptr;
  }

  // Create an empty solver
  SnsAccIkBase::uPtr accIk(new SnsAccIkBase(nJnt));

  // Set the joint limits:
  if (!accIk->setBounds(ddqLow, ddqUpp)) { ROS_ERROR("Bad Input!"); return nullptr; };

  return accIk;
}

/*************************************************************************************************/

SnsIkBase::ExitCode SnsAccIkBase::solve(const Eigen::MatrixXd& J, const Eigen::VectorXd& dJdq,
                                           const Eigen::VectorXd& ddx, Eigen::VectorXd* ddq, double* taskScale)
{
  // Input validation
  if (!ddq) { ROS_ERROR("ddq is nullptr!"); return ExitCode::BadUserInput; }
  if (!taskScale) { ROS_ERROR("taskScale is nullptr!"); return ExitCode::BadUserInput; }
  size_t nTask = ddx.size();
  if (nTask <= 0) {
    ROS_ERROR("Bad Input: ddx.size() > 0 is required!");
    return ExitCode::BadUserInput;
  }
  if (size_t(J.rows()) != nTask) {
    ROS_ERROR("Bad Input: J.rows() == ddx.size() is required!");
    return ExitCode::BadUserInput;
  }
  if (size_t(J.cols()) != getNrOfJoints()) {
    ROS_ERROR("Bad Input: J.cols() == nJnt is required!");
    return ExitCode::BadUserInput;
  }
  size_t dJdqDim = dJdq.size();
  if (dJdqDim != nTask) {
    ROS_ERROR("Bad Input: dJdq.size() == nTask is required!");
    return ExitCode::BadUserInput;
  }

  // Local variable initialization:
  Eigen::MatrixXd W = Eigen::MatrixXd::Identity(getNrOfJoints(), getNrOfJoints());  // null-space selection matrix
  Eigen::VectorXd ddqNull = Eigen::VectorXd::Zero(getNrOfJoints());  // acceleration in the null-space
  *taskScale = 1.0;  // task scale (assume feasible solution until proven otherwise)

  // Temp. variables to store the best solution
  double bestTaskScale = 0.0;  // temp variable to track the lower bound on the task scale between iterations
  Eigen::MatrixXd bestW;  // temp variable to track W before it has been accepted
  Eigen::VectorXd bestDdqNull;  // temp variable to track dqNull between iterations

  // Set the linear solver for this iteration:
  if(setLinearSolver(J*W) != ExitCode::Success) {
    ROS_ERROR("Solver failed to set linear solver!");
    return ExitCode::InternalError;
  }

  // Keep track of which joints are saturated:
  std::vector<bool> jointIsFree(getNrOfJoints(), true);

  // Main solver loop:
  double resErr;  // residual error in the linear solver
  for (size_t iter = 0; iter < getNrOfJoints() * MAXIMUM_SOLVER_ITERATION_FACTOR; iter++) {

    // Compute the joint acceleration given current saturation set:
    if (solveProjectionEquation(J, dJdq, ddqNull, ddx, ddq, &resErr) != ExitCode::Success) {
      ROS_ERROR("Failed to solve projection equation!");
      return ExitCode::InternalError;
    }
    if (resErr > LIN_SOLVE_RESIDUAL_TOL) { // check that the solver found a feasible solution
      ROS_ERROR("Task is infeasible!  resErr: %e > tol: %e", resErr, LIN_SOLVE_RESIDUAL_TOL);
      return ExitCode::InfeasibleTask;
    }

    // Check to see if the solution satisfies the joint limits
    if (checkBounds(*ddq)) { // Done! solution is feasible and task scale is at maximum value
      return ExitCode::Success;
    }  //  else joint acceleration is infeasible: saturate joint and then try again

    // Compute the task scaling factor
    double tmpScale;
    int jntIdx;
    ExitCode taskScaleExit = computeTaskScalingFactor(J, ddx, *ddq, jointIsFree, &tmpScale, &jntIdx, &resErr);
    if (resErr > LIN_SOLVE_RESIDUAL_TOL) { // check that the solver found a feasible solution
      ROS_ERROR("Failed to compute task scale!  resErr: %e > tol: %e", resErr, LIN_SOLVE_RESIDUAL_TOL);
      return ExitCode::InfeasibleTask;
    }
    if (taskScaleExit != ExitCode::Success) {
      ROS_ERROR("Failed to compute task scale!");
      return taskScaleExit;
    }
    if (tmpScale < MINIMUM_FINITE_SCALE_FACTOR) { // check that the solver found a feasible solution
      ROS_ERROR("Task is infeasible! scaling --> zero");
      return ExitCode::InfeasibleTask;
    }

    if (tmpScale > 1.0) {
      ROS_ERROR("Task scale is %f, which is more than 1.0", tmpScale);
      return ExitCode::InternalError;
    }

    // If the task scale exceeds previous, then cache the results as "best so far"
    // Also if the current best so far solution violates the limits, update bestTakeScale
    Eigen::VectorXd ddxScaledTmp = (ddx.array() * bestTaskScale).matrix();
    Eigen::VectorXd ddqTmp;
    if (solveProjectionEquation(J, dJdq, ddqNull, ddxScaledTmp, &ddqTmp, &resErr) != ExitCode::Success) {
      ROS_ERROR("Failed to solve projection equation!");
      return ExitCode::InternalError;
    }
    if (resErr > LIN_SOLVE_RESIDUAL_TOL) { // check that the solver found a feasible solution
      ROS_ERROR("Task is infeasible!  resErr: %e > tol: %e", resErr, LIN_SOLVE_RESIDUAL_TOL);
      return ExitCode::InfeasibleTask;
    }
    if (tmpScale > bestTaskScale || !checkBounds(ddqTmp)) {
      bestTaskScale = tmpScale;
      bestW = W;
      bestDdqNull = ddqNull;
    }

    // Saturate the most critical joint
    W(jntIdx, jntIdx) = 0.0;
    jointIsFree[jntIdx] = false;
    if ((*ddq)(jntIdx) > (getUpperBounds())(jntIdx)) {
      ddqNull(jntIdx) = (getUpperBounds())(jntIdx);
    } else if ((*ddq)(jntIdx) < (getLowerBounds())(jntIdx)) {
      ddqNull(jntIdx) = (getLowerBounds())(jntIdx);
    } else {
      ROS_ERROR("Internal error in computing task scale!  ddq(%d) = %f", jntIdx, (*ddq)(jntIdx));
      return ExitCode::InternalError;
    }

    // Update the linear solver
    if (setLinearSolver(J*W) != ExitCode::Success) {
      ROS_ERROR("Solver failed to set linear solver!");
      return ExitCode::InternalError;
    }

    // Test the rank:
    if (getLinSolverRank() < nTask) { // no more degrees of freedom: scale the task
      (*taskScale) = bestTaskScale;
      W = bestW;
      ddqNull = bestDdqNull;

      // Update the linear solver
      if (setLinearSolver(J * W) != ExitCode::Success) {
        ROS_ERROR("Solver failed to set linear solver!");
        return ExitCode::InternalError;
      }

      // Compute the joint acceleration given current saturation set:
      Eigen::VectorXd ddxScaled = (ddx.array() * (*taskScale)).matrix();
      if (solveProjectionEquation(J, dJdq, ddqNull, ddxScaled, ddq, &resErr) != ExitCode::Success) {
        ROS_ERROR("Failed to solve projection equation!");
        return ExitCode::InternalError;
      }
      if (resErr > LIN_SOLVE_RESIDUAL_TOL) { // check that the solver found a feasible solution
        ROS_ERROR("Task is infeasible!  resErr: %e > tol: %e", resErr, LIN_SOLVE_RESIDUAL_TOL);
        return ExitCode::InfeasibleTask;
      }

      return ExitCode::Success;  // DONE

    } // end rank test

  }  // end main solver loop

  ROS_ERROR("Internal Error: reached maximum iteration in solver main loop!");
  return ExitCode::InternalError;
}

/*************************************************************************************************
 *                               Protected Methods                                               *
 *************************************************************************************************/

}  // namespace sns_ik
