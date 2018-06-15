/** @file sns_vel_ik_base.cpp
 *
 * @brief The file provides the basic implementation of the SNS-IK velocity solver
 *
 * @author Matthew Kelly
 *
 * This file provides a set of functions that return simple kinematic chains that are used for the
 * unit tests. This allows unit tests to run quickly without depending on external URDF files.
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
#include <sns_ik/sns_vel_ik_base.hpp>

#include <ros/console.h>
#include <limits>

namespace sns_ik {

/*
 * The code of the SNS-IK solver relies on a linear solver. If the linear system is infeasible,
 * then the solver will return the minimum-norm solution with a non-zero (positive) residual.
 * This corresponds to a situation where there is no feasible solution to the task objective, for
 * example, a one-link robot arm where the task velocity is along the axis of the joint. This
 * threshold is used to determine if the task is infeasible.
 */
static const double LIN_SOLVE_RESIDUAL_TOL = 1e-8;

/*
 * If all goes well, then the solver should always terminate the main loop. If something crazy
 * happens, then we want to prevent an infinite loop. The maximum iteration count is given by
 * this factor multiplied by the number of joints in the task jacobian.
 */
static const int MAXIMUM_SOLVER_ITERATION_FACTOR = 100;

// Short names for useful constants
static const double POS_INF = std::numeric_limits<double>::max();
static const double NEG_INF = std::numeric_limits<double>::lowest();

// Any number smaller than this is considered to be zero for numerical purposes
static const double MINIMUM_FINITE_SCALE_FACTOR = 1e-10;

// Any number larger than this is considered to be infinite for numerical purposes
static const double MAXIMUM_FINITE_SCALE_FACTOR = 1e10;

// Tolerance for checks on the velocity limits
static const double VELOCITY_BOUND_TOLERANCE = 1e-8;

/*************************************************************************************************
 *                                 Public Methods                                                *
 *************************************************************************************************/

SnsVelIkBase::uPtr SnsVelIkBase::create(int nJnt)
{
  if (nJnt <= 0) {
    ROS_ERROR("Bad Input: dqLow.size(%d) > 0 is required!", nJnt);
    return nullptr;
  }
  Eigen::ArrayXd dqLow = NEG_INF*Eigen::ArrayXd::Ones(nJnt);
  Eigen::ArrayXd dqUpp = POS_INF*Eigen::ArrayXd::Ones(nJnt);
  return create(dqLow, dqUpp);
}

/*************************************************************************************************/

SnsVelIkBase::uPtr SnsVelIkBase::create(const Eigen::ArrayXd& dqLow, const Eigen::ArrayXd& dqUpp)
{
  // Input validation
  int nJnt = dqLow.size();
  if (nJnt <= 0) {
    ROS_ERROR("Bad Input: dqLow.size(%d) > 0 is required!", nJnt);
    return nullptr;
  }

  // Create an empty solver
  SnsVelIkBase::uPtr velIk(new SnsVelIkBase(nJnt));

  // Set the joint limits:
  if (!velIk->setVelBnd(dqLow, dqUpp)) { ROS_ERROR("Bad Input!"); return nullptr; };

  return velIk;
}

/*************************************************************************************************/

bool SnsVelIkBase::setVelBnd(const Eigen::ArrayXd& dqLow, const Eigen::ArrayXd& dqUpp)
{
  int nJnt = dqLow.size();
  if (nJnt <= 0) {
    ROS_ERROR("Bad Input: dqLow.size(%d) > 0 is required!", nJnt);
    return false;
  }
  if (dqLow.size() != dqUpp.size()) {
    ROS_ERROR("Bad Input: dqLow.size(%d) == dqUpp.size(%d) is required!",
              int(dqLow.size()), int(dqUpp.size()));
    return false;
  }
  nJnt_ = nJnt;
  dqLow_ = dqLow;
  dqUpp_ = dqUpp;
  return true;
}

/*************************************************************************************************/

SnsVelIkBase::ExitCode SnsVelIkBase::solve(const Eigen::MatrixXd& J, const Eigen::VectorXd& dx,
                               Eigen::VectorXd* dq, double* taskScale)
{
  // Input validation
  if (!dq) { ROS_ERROR("dq is nullptr!"); return ExitCode::BadUserInput; }
  if (!taskScale) { ROS_ERROR("taskScale is nullptr!"); return ExitCode::BadUserInput; }
  int nTask = dx.size();
  if (nTask <= 0) {
    ROS_ERROR("Bad Input: dx.size() > 0 is required!");
    return ExitCode::BadUserInput;
  }
  if (int(J.rows()) != nTask) {
    ROS_ERROR("Bad Input: J.rows() == dx.size() is required!");
    return ExitCode::BadUserInput;
  }
  if (int(J.cols()) != nJnt_) {
    ROS_ERROR("Bad Input: J.cols() == nJnt is required!");
    return ExitCode::BadUserInput;
  }

  /*
   * W is a diagonal selection matrix which indicates free joints.
   * The entry of 1 indicates the corresponding joint is free for the task,
   * and 0 indicates the corresponding joint is saturated.
   */
  Eigen::MatrixXd W = Eigen::MatrixXd::Identity(nJnt_, nJnt_);  // null-space selection matrix
  Eigen::VectorXd dqNull = Eigen::VectorXd::Zero(nJnt_);  // velocity in the null-space
  *taskScale = 1.0;  // task scale (assume feasible solution until proven otherwise)

  // Temp. variables to store the best solution
  double bestTaskScale = 0.0;  // temp variable to track the lower bound on the task scale between iterations
  Eigen::MatrixXd bestW;  // temp variable to track W before it has been accepted
  Eigen::VectorXd bestDqNull;  // temp variable to track dqNull between iterations

  // Set the linear solver for this iteration:
  Eigen::MatrixXd JW = J*W;
  linSolver_ = SnsLinearSolver(JW);
  if(linSolver_.info() != Eigen::ComputationInfo::Success) {
    ROS_ERROR("Solver failed to decompose the combined sparse matrix!");
    return ExitCode::InternalError;
  }

  // Keep track of which joints are saturated:
  std::vector<bool> jointIsFree(nJnt_, true);

  // Main solver loop:
  double resErr;  // residual error in the linear solver
  for (int iter = 0; iter < nJnt_ * MAXIMUM_SOLVER_ITERATION_FACTOR; iter++) {

    // Compute the joint velocity given current saturation set:
    if (!solveProjectionEquation(J, JW, dqNull, dx, dq, &resErr)) {
      ROS_ERROR("Failed to solve projection equation!");
      return ExitCode::InternalError;
    }
    if (resErr > LIN_SOLVE_RESIDUAL_TOL) { // check that the solver found a feasible solution
      ROS_ERROR("Task is infeasible!  resErr: %e > tol: %e", resErr, LIN_SOLVE_RESIDUAL_TOL);
      return ExitCode::InfeasibleTask;
    }

    // Check to see if the solution satisfies the joint limits
    if (checkVelBnd(*dq)) { // Done! solution is feasible and task scale is at maximum value
      return ExitCode::Success;
    }  //  else joint velocity is infeasible: saturate joint and then try again

    // Compute the task scaling factor
    double tmpScale;
    int jntIdx;
    ExitCode taskScaleExit = computeTaskScalingFactor(J, JW, dx, *dq, jointIsFree, &tmpScale, &jntIdx, &resErr);
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
    if (tmpScale > bestTaskScale) {
      bestTaskScale = tmpScale;
      bestW = W;
      bestDqNull = dqNull;
    }

    // Saturate the most critical joint
    W(jntIdx, jntIdx) = 0.0;
    jointIsFree[jntIdx] = false;
    if ((*dq)(jntIdx) > dqUpp_(jntIdx)) {
      dqNull(jntIdx) = dqUpp_(jntIdx);
    } else if ((*dq)(jntIdx) < dqLow_(jntIdx)) {
      dqNull(jntIdx) = dqLow_(jntIdx);
    } else {
      ROS_ERROR("Internal error in computing task scale!  dq(%d) = %f", jntIdx, (*dq)(jntIdx));
      return ExitCode::InternalError;
    }

    // Update the linear solver
    JW = J*W;
    linSolver_ = SnsLinearSolver(JW);
    if(linSolver_.info() != Eigen::ComputationInfo::Success) {
      ROS_ERROR("Solver failed to decompose the combined sparse matrix!");
      return ExitCode::InternalError;
    }

    // Test the rank:
    if (linSolver_.rank() < nTask) { // no more degrees of freedom: scale the task
      (*taskScale) = bestTaskScale;
      W = bestW;
      dqNull = bestDqNull;

      // Update the linear solver
      JW = J*W;
      linSolver_ = SnsLinearSolver(JW);
      if(linSolver_.info() != Eigen::ComputationInfo::Success) {
        ROS_ERROR("Solver failed to decompose the combined sparse matrix!");
        return ExitCode::InternalError;
      }

      // Compute the joint velocity given current saturation set:
      Eigen::VectorXd dxScaled = (dx.array() * (*taskScale)).matrix();
      if (!solveProjectionEquation(J, JW, dqNull, dxScaled, dq, &resErr)) {
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

bool SnsVelIkBase::checkVelBnd(const Eigen::VectorXd& dq)
{
  if (dq.size() != nJnt_) {
    ROS_ERROR("Bad Input:  dq.size(%d) == nJnt(%d) is required!", int(dq.size()), nJnt_);
    return false;
  }
  for (int i = 0; i < nJnt_; i++) {
    if (dq(i) < dqLow_(i) - VELOCITY_BOUND_TOLERANCE) { return false; }
    if (dq(i) > dqUpp_(i) + VELOCITY_BOUND_TOLERANCE) { return false; }
  }
  return true;
}

/*************************************************************************************************/

bool SnsVelIkBase::solveProjectionEquation(const Eigen::MatrixXd& J, const Eigen::MatrixXd& JW,
                                           const Eigen::VectorXd& dqNull, const Eigen::VectorXd& dx,
                                           Eigen::VectorXd* dq, double* resErr)
{
  // Input validation:
  if (J.cols() != JW.cols()) { ROS_ERROR("Bad Input!  J.cols() != JW.cols()"); return false; }
  if (J.rows() != JW.rows()) { ROS_ERROR("Bad Input!  J.rows() != JW.rows()"); return false; }
  if (dx.size() != J.rows()) { ROS_ERROR("Bad Input!  dx.size() != J.rows()"); return false; }
  if (J.cols() != dqNull.rows()) { ROS_ERROR("Bad Input!  J.cols() != dqNull.rows()"); return false; }
  if (!dq) { ROS_ERROR("Bad Input!  dq is nullptr!"); return false; }
  if (!resErr) { ROS_ERROR("Bad Input!  resErr is nullptr!"); return false; }

  // Solve the linear system
  Eigen::MatrixXd B = dx - J*dqNull;
  Eigen::MatrixXd X = linSolver_.solve(B);
  if(linSolver_.info() != Eigen::ComputationInfo::Success) {
    ROS_ERROR("Failed to solve projection equations!");
    return false;
  }

  // Solve for dq
  *dq = X + dqNull;
  *resErr = (JW*X - B).squaredNorm();
  return true;
}

/*************************************************************************************************/

SnsVelIkBase::ExitCode SnsVelIkBase::computeTaskScalingFactor(const Eigen::MatrixXd& J, const Eigen::MatrixXd& JW,
                                            const Eigen::VectorXd& dx, const Eigen::VectorXd& dq,
                                            const std::vector<bool>& jntIsFree,
                                            double* taskScale, int* jntIdx, double* resErr)
{
  if (J.cols() != JW.cols()) { ROS_ERROR("Bad Input!  J.cols() != JW.cols()"); return ExitCode::BadUserInput; }
  if (J.rows() != JW.rows()) { ROS_ERROR("Bad Input!  J.rows() != JW.rows()"); return ExitCode::BadUserInput; }
  if (dx.size() != J.rows()) { ROS_ERROR("Bad Input!  dx.size() != J.rows()"); return ExitCode::BadUserInput; }
  if (dq.size() != J.cols()) { ROS_ERROR("Bad Input!  dq.size() != J.cols()"); return ExitCode::BadUserInput; }
  if (!taskScale) { ROS_ERROR("taskScale is nullptr!"); return ExitCode::BadUserInput; }
  if (!jntIdx) { ROS_ERROR("jntIdx is nullptr!"); return ExitCode::BadUserInput; }
  if (!resErr) { ROS_ERROR("resErr is nullptr!"); return ExitCode::BadUserInput; }

  // Compute "a" and "b" from the paper.   (J*W*a = dx)
  Eigen::VectorXd a = linSolver_.solve(dx);
  if(linSolver_.info() != Eigen::ComputationInfo::Success) {
    ROS_ERROR("Failed to solve projection equations!");
    return ExitCode::InternalError;
  }
  *resErr = (JW*a - dx).squaredNorm();
  Eigen::ArrayXd b = (dq - a).array();

  // Compute the task scale associated with each joint
  Eigen::ArrayXd jntScaleFactorArr(nJnt_);
  Eigen::ArrayXd lowMargin = (dqLow_ - b);
  Eigen::ArrayXd uppMargin = (dqUpp_ - b);
  for (int i = 0; i < nJnt_; i++) {
    if (jntIsFree[i]) {
      jntScaleFactorArr(i) = SnsVelIkBase::findScaleFactor(lowMargin(i), uppMargin(i), a(i));
    } else {  // joint is constrained
      jntScaleFactorArr(i) = POS_INF;
    }
  }

  // Compute the most critical scale factor and corresponding joint index
  *jntIdx = 0;  // index of the most critical joint
  *taskScale = jntScaleFactorArr(*jntIdx);  // minimum value of jntScaleFactorArr()
  for (int i = 1; i < nJnt_; i++) {
    if (jntScaleFactorArr(i) < *taskScale) {
      *jntIdx = i;
      *taskScale = jntScaleFactorArr(i);
    }
  }
  return ExitCode::Success;
}

/*************************************************************************************************/

double SnsVelIkBase::findScaleFactor(double low, double upp, double a)
{
  if (std::abs(a) > MAXIMUM_FINITE_SCALE_FACTOR) {
    return 0.0;
  }
  if (a < 0.0 && low < 0.0) {
    if (a < low) {
      return low / std::min(a, -MINIMUM_FINITE_SCALE_FACTOR);
    } else {
      return 1.0;  // feasible without scaling!
    }
  } else if (a > 0.0 && upp > 0.0) {
    if (upp < a) {
      return upp / std::max(a, MINIMUM_FINITE_SCALE_FACTOR);
    } else {
      return 1.0;  // feasible without scaling!
    }
  } else {
    return 0.0;  // infeasible
  }
}

/*************************************************************************************************/

}  // namespace sns_ik
