/** @file sns_linear_solver.hpp
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

#ifndef SNS_IK_LIB__SNS_LINEAR_SOLVER_H_
#define SNS_IK_LIB__SNS_LINEAR_SOLVER_H_

#include <Eigen/Dense>

namespace sns_ik {

#if EIGEN_VERSION_AT_LEAST(3,3,4)  //  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
// Eigen version is newer than 3.3.4: CompleteOrthogonalDecomposition is defined
typedef Eigen::CompleteOrthogonalDecomposition<Eigen::MatrixXd> SnsLinearSolver;

#else  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
// Eigen version is older than 3.3.4: CompleteOrthogonalDecomposition is not defined
// Implement much of the API for CompleteOrthogonalDecomposition, but use the same back-end as
// the original implementation of the SNS-IK solver.
// https://eigen.tuxfamily.org/dox/classEigen_1_1CompleteOrthogonalDecomposition.html

class SnsLinearSolver {

public:

  /*
   * Create a default linear solver
   */
  SnsLinearSolver();

  /*
   * Create a default linear solver with memory preallocation
   */
  SnsLinearSolver(int n, int m);

  /*
   * Create a linear solver for the matrix A*x = b
   * @param A: matrix of interest
   */
  SnsLinearSolver(const Eigen::MatrixXd& A);

  /*
   * Set and decompose the matrix in the linear system.
   */
  void compute(const Eigen::MatrixXd& A);

  /*
   * Find x to minimize:  ||A*x - b||^2
   * @param b: right hand side of the linear system
   * @return: x = solution to the optimization problem
   */
  Eigen::MatrixXd solve(const Eigen::MatrixXd& b);

  /*
   * @return: status of the solver
   */
  Eigen::ComputationInfo info() const { return info_; };

  /*
   * Set the threshold that is used for computing rank and pseudoinverse
   * @param tol: threshold used for decomposing matrix and computing rank
   */
  void setThreshold(Eigen::Default_t tol);

  /*
   * Set the threshold that is used for computing rank and pseudoinverse
   * @param tol: threshold used for decomposing matrix and computing rank
   */
  void setThreshold(double tol) { tol_ = tol; };

  /*
   * @return: the rank of the matrix A
   */
  unsigned int rank() const { return rank_; };

private:

  // status of the most recent solve operation
  Eigen::ComputationInfo info_;

  // the "A" matrix in A*x = b
  Eigen::MatrixXd A_;

  // the pseudo-inverse of "A"
  Eigen::MatrixXd invA_;

  // tolerance that is used for checking for non-trivial singular values
  double tol_;

  // rank of "A"
  int rank_;

};

#endif  // EIGEN_VERSION_AT_LEAST(3,3,4)  //- - - - - - - - - - - - - - - - - - - - - - - - - - //

}  // namespace sns_ik

#endif // SNS_IK_LIB__SNS_LINEAR_SOLVER_H_
