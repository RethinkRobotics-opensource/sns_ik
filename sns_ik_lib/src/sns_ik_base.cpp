/** @file sns_ik_base.cpp
 *
 * @brief This file provides the common functions shared by SNS IK base solvers
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
#include <sns_ik/sns_ik_base.hpp>

#include <ros/console.h>
#include <limits>

namespace sns_ik {

const double SnsIkBase::LIN_SOLVE_RESIDUAL_TOL = 1e-8;
const double SnsIkBase::PINV_TOL = 1e-10;
const int SnsIkBase::MAXIMUM_SOLVER_ITERATION_FACTOR = 100;
const double SnsIkBase::POS_INF = std::numeric_limits<double>::max();
const double SnsIkBase::NEG_INF = std::numeric_limits<double>::lowest();
const double SnsIkBase::MINIMUM_FINITE_SCALE_FACTOR = 1e-10;
const double SnsIkBase::MAXIMUM_FINITE_SCALE_FACTOR = 1e10;
const double SnsIkBase::BOUND_TOLERANCE = 1e-8;
const Eigen::IOFormat SnsIkBase::EigArrFmt(4, 0, ", ", "\n", "[", "]");

/*************************************************************************************************
 *                                 Public Methods                                                *
 *************************************************************************************************/

bool SnsIkBase::setBounds(const Eigen::ArrayXd& qLow, const Eigen::ArrayXd& qUpp)
{
  int nJnt = qLow.size();
  if (nJnt <= 0) {
    ROS_ERROR("Bad Input: qLow.size(%d) > 0 is required!", nJnt);
    return false;
  }
  if (qLow.size() != qUpp.size()) {
    ROS_ERROR("Bad Input: qLow.size(%d) == qUpp.size(%d) is required!",
              int(qLow.size()), int(qUpp.size()));
    return false;
  }
  nJnt_ = nJnt;
  qLow_ = qLow;
  qUpp_ = qUpp;
  return true;
}

/*************************************************************************************************
 *                               Protected Methods                                               *
 *************************************************************************************************/

bool SnsIkBase::checkBounds(const Eigen::VectorXd& q)
{
  if (q.size() != nJnt_) {
    ROS_ERROR("Bad Input:  q.size(%d) == nJnt(%d) is required!", int(q.size()), nJnt_);
    return false;
  }
  for (int i = 0; i < nJnt_; i++) {
    if (q(i) < qLow_(i) - BOUND_TOLERANCE) { return false; }
    if (q(i) > qUpp_(i) + BOUND_TOLERANCE) { return false; }
  }
  return true;
}

/*************************************************************************************************/

SnsIkBase::ExitCode SnsIkBase::setLinearSolver(const Eigen::MatrixXd& JW)
{
  linSolver_ = SnsLinearSolver(JW);
  if(linSolver_.info() != Eigen::ComputationInfo::Success) {
    ROS_ERROR("Solver failed to decompose the matrix!");
    return ExitCode::InternalError;
  }
  JW_ = JW;  // store the matrix that was decomposed - used for computing the residual error
  return ExitCode::Success;
}

/*************************************************************************************************/

SnsIkBase::ExitCode SnsIkBase::solveLinearSystem(const Eigen::MatrixXd& rhs, Eigen::VectorXd* q, double* resErr)
{
  if (!q) {
    ROS_ERROR("q is nullptr!");
    return ExitCode::BadUserInput;
  }
  if (JW_.size() == 0) {
    ROS_ERROR("Cannot solve an empty system! Have you called setLinearSolver()?");
    return ExitCode::BadUserInput;
  }
  if (rhs.rows() != JW_.rows()) {
    ROS_ERROR("Invalid matrix dimensions! rhs.rows() == JW_.rows(). Linear system is inconsistent.");
    return ExitCode::BadUserInput;
  }
  *q = linSolver_.solve(rhs);
  if(linSolver_.info() != Eigen::ComputationInfo::Success) {
    ROS_ERROR("Failed to solve linear system!");
    return ExitCode::InfeasibleTask;
  }
  if (resErr) {
    *resErr = (JW_*(*q) - rhs).squaredNorm();
  }
  return ExitCode::Success;
}

/*************************************************************************************************/

SnsIkBase::ExitCode SnsIkBase::solveProjectionEquation(const Eigen::MatrixXd& J, const Eigen::VectorXd& dqNull,
                                                       const Eigen::VectorXd& dx, Eigen::VectorXd* dq, double* resErr)
{
  // Input validation:
  if (dx.size() != J.rows()) { ROS_ERROR("Bad Input!  dx.size() != J.rows()"); return ExitCode::BadUserInput; }
  if (J.cols() != dqNull.rows()) { ROS_ERROR("Bad Input!  J.cols() != dqNull.rows()"); return ExitCode::BadUserInput; }
  if (!dq) { ROS_ERROR("Bad Input!  dq is nullptr!"); return ExitCode::BadUserInput; }
  if (!resErr) { ROS_ERROR("Bad Input!  resErr is nullptr!"); return ExitCode::BadUserInput; }

  // Solve the linear system
  Eigen::MatrixXd B = dx - J*dqNull;
  ExitCode result = solveLinearSystem(B, dq, resErr);
  if (result != ExitCode::Success) {
    ROS_ERROR("Failed to solve linear system!");
  }

  // Solve for dq
  *dq = *dq + dqNull;
  return result;
}


/*************************************************************************************************/

SnsIkBase::ExitCode SnsIkBase::solveProjectionEquation(const Eigen::MatrixXd& J, const Eigen::VectorXd& dJdq,
                                                       const Eigen::VectorXd& ddqNull, const Eigen::VectorXd& ddx,
                                                       Eigen::VectorXd* ddq, double* resErr)
{
  /// Input validation:
  if (ddx.size() != J.rows()) { ROS_ERROR("Bad Input!  ddx.size() != J.rows()"); return ExitCode::BadUserInput; }
  if (ddx.size() != dJdq.rows()) { ROS_ERROR("Bad Input!  ddx.size() != dJdq.rows()"); return ExitCode::BadUserInput; }
  if (J.cols() != ddqNull.rows()) { ROS_ERROR("Bad Input!  J.cols() != ddqNull.rows()"); return ExitCode::BadUserInput; }
  if (!ddq) { ROS_ERROR("Bad Input!  ddq is nullptr!"); return ExitCode::BadUserInput; }
  if (!resErr) { ROS_ERROR("Bad Input!  resErr is nullptr!"); return ExitCode::BadUserInput; }

  // Solve the linear system
  Eigen::MatrixXd B = ddx - dJdq - J*ddqNull;
  ExitCode result = solveLinearSystem(B, ddq, resErr);
  if (result != ExitCode::Success) {
    ROS_ERROR("Failed to solve linear system!");
  }

  // Solve for ddq
  *ddq = *ddq + ddqNull;
  return result;
}

/*************************************************************************************************/

SnsIkBase::ExitCode SnsIkBase::computeTaskScalingFactor(const Eigen::MatrixXd& J,
                                            const Eigen::VectorXd& desiredTask, const Eigen::VectorXd& jointOut,
                                            const std::vector<bool>& jntIsFree,
                                            double* taskScale, int* jntIdx, double* resErr)
{
  if (desiredTask.size() != J.rows()) { ROS_ERROR("Bad Input!  desiredTask.size() != J.rows()"); return ExitCode::BadUserInput; }
  if (jointOut.size() != J.cols()) { ROS_ERROR("Bad Input!  jointOut.size() != J.cols()"); return ExitCode::BadUserInput; }
  if (!taskScale) { ROS_ERROR("taskScale is nullptr!"); return ExitCode::BadUserInput; }
  if (!jntIdx) { ROS_ERROR("jntIdx is nullptr!"); return ExitCode::BadUserInput; }
  if (!resErr) { ROS_ERROR("resErr is nullptr!"); return ExitCode::BadUserInput; }

  // Compute "a" and "b" from the paper.   (J*W*a = dx)
  Eigen::VectorXd a;
  ExitCode result = solveLinearSystem(desiredTask, &a, resErr);
  if (result != ExitCode::Success) {
    ROS_ERROR("Failed to solve linear system!");
    return result;
  }
  Eigen::ArrayXd b = (jointOut - a).array();

  // Compute the task scale associated with each joint
  Eigen::ArrayXd jntScaleFactorArr(nJnt_);
  Eigen::ArrayXd lowMargin = (qLow_ - b);
  Eigen::ArrayXd uppMargin = (qUpp_ - b);
  for (int i = 0; i < nJnt_; i++) {
    if (jntIsFree[i]) {
      jntScaleFactorArr(i) = SnsIkBase::findScaleFactor(lowMargin(i), uppMargin(i), a(i));
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

double SnsIkBase::findScaleFactor(double low, double upp, double a)
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