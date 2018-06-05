/** @file sns_linear_solver.cpp
 *
 * @brief Linear Solver, used by the SNS-IK velocity IK solvers.
 * @author: Matthew Kelly
 */

/**
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

/*
 *
 * This class is used to solve linear systems. It is a wrapper for two different internal solvers:
 *   --> psuedo-inverse solver: included for legacy support on Eigen 3.2.0
 *   --> direct linear solver: prefered solver when available. Requries Eigen 3.3.4
 */

#include "sns_linear_solver.hpp"

#if EIGEN_VERSION_AT_LEAST(3,3,4)  //  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
// Eigen version is newer than 3.3.4: CompleteOrthogonalDecomposition is defined
// We're done - just use the implementation from Eigen

#else  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
// Eigen version is older than 3.3.4: CompleteOrthogonalDecomposition is not defined
// Implement much of the API for CompleteOrthogonalDecomposition, but use the same back-end as
// the original implementation of the SNS-IK solver.
// https://eigen.tuxfamily.org/dox/classEigen_1_1CompleteOrthogonalDecomposition.html

#include "sns_ik_math_utils.hpp"
#include <ros/console.h>

namespace sns_ik {

static const double DEFAULT_THRESHOLD = 1e-8;

SnsLinearSolver::SnsLinearSolver()
  : info_(Eigen::Success), rank_(0)
{
  setThreshold(Eigen::Default);
}

SnsLinearSolver::SnsLinearSolver(int n, int m)
  : info_(Eigen::Success), A_(n, m), invA_(m, n), rank_(0)
{
  setThreshold(Eigen::Default);
}

SnsLinearSolver::SnsLinearSolver(const Eigen::MatrixXd& A)
  : SnsLinearSolver(A.rows(), A.cols())
{
  setThreshold(Eigen::Default); compute(A);
}

Eigen::MatrixXd SnsLinearSolver::solve(const Eigen::MatrixXd& b)
{
  static bool PRINT_WARNING = true;
  if (PRINT_WARNING) {
    ROS_WARN("Using legacy code. Upgrade to at least Eigen 3.3.4 for improved linear solver.");
    PRINT_WARNING = false;  // only print the warning once.
  }
  return invA_ * b;
}

void SnsLinearSolver::compute(const Eigen::MatrixXd& A)
{
  A_ = A;
  if (!sns_ik::pseudoInverse(A_, tol_, &invA_, &rank_)) {
    info_ = Eigen::InvalidInput;
  }
}

void SnsLinearSolver::setThreshold(Eigen::Default_t tol) {
  setThreshold(DEFAULT_THRESHOLD);
}

}  // namespace sns_ik

#endif  // EIGEN_VERSION_AT_LEAST(3,3,4)  //- - - - - - - - - - - - - - - - - - - - - - - - - - //
