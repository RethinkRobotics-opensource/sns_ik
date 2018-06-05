/*! \file sns_ikl_math_utils.hpp
 * \brief Math utilities for the SNS IK solvers
 * \author Fabrizio Flacco
 * \author Forrest Rogers-Marcovitz
 * \author Matthew Kelly
 */
/*
 *    Copyright 2016-2018 Rethink Robotics
 *
 *    Copyright 2012-2016 Fabrizio Flacco
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

#ifndef SNS_IKL_MATH_UTILS
#define SNS_IKL_MATH_UTILS

#include <Eigen/Dense>

#include "sns_ik_math_utils.hpp"

namespace sns_ik {

static const double INF = std::numeric_limits<double>::max();

/*
 * FIXME:  Is it possible to avoid doing all of these inverse operations? It is far better
 *         numerically to solve a linear equation, rather than compute an inverse and then
 *         perform a matrix multiply.
 */

/*
 * FIXME:  There is no input validation on any of this code. If a matrix is input with the wrong
 *         size then it simply hits an assertion and crashes the thread. This is not desirable
 *         in a core solver that is running on a robot.
 */

/*
 * Compute the pseudoinverse of A using an algorithm based on singular value decomposition.
 * @param A: input matrix
 *   FIXME: if A.rows() >= A.cols() causes a failed assertion
 * @param[out] invA: pseudoinverse of A
 * @param[opt] eps: singular values smaller than this will be set to zero
 * @return: true if A is full rank, false if A is rank deficient
 */
bool pinv(const Eigen::MatrixXd &A, Eigen::MatrixXd *invA, double eps = 1e-6);

/*
 * Compute the pseudoinverse of A along with the nullspace projector matrix.
 * @param A: input matrix
 *   FIXME: if A.rows() >= A.cols() causes a failed assertion, but this should be valid input
 * @param[out] invA: pseudoinverse of A
 * @param[in/out] P: the nullspace projector matrix:  P = (P - pinv(A)*A)
 *                   P.rows() == P.cols() = A.cols() is required
 *                   P should be initialized with the identity matrix
 * @param[opt] eps: singular values smaller than this will be set to zero
 * @return: true if A is full rank, false if A is rank deficient
 */
bool pinv_P(const Eigen::MatrixXd &A, Eigen::MatrixXd *invA, Eigen::MatrixXd *P, double eps = 1e-6);

/*
 * Compute the pseudoinverse of A along with the nullspace projector matrix.
 * @param A: input matrix
 *   FIXME: if A.rows() >= A.cols() causes a failed assertion, but this should be valid input
 * @param[out] invA: pseudoinverse of A
 * @param[in/out/opt] P: the nullspace projector matrix:  P = (P - pinv(A)*A)
 *                   P.rows() == P.cols() = A.cols() is required
 *                   P should be initialized with the identity matrix
 * @param[opt] lambda_max: damping parameter for the damped pseudoinverse
 * @param[opt] eps: singular values smaller than this will be set to zero
 * @return: true if A is full rank, false if A is rank deficient
 */
bool pinv_damped_P(const Eigen::MatrixXd &A, Eigen::MatrixXd *invA, Eigen::MatrixXd *P = nullptr,
                   double lambda_max = 1e-6, double eps = 1e-6);

/*
 * Compute the pseudoinverse of A using an algorithm based on QR decomposition
 * @param A: input matrix
 *   FIXME: if A.rows() >= A.cols() causes a failed assertion, but this should be valid input
 * @param[out] invA: pseudoinverse of A
 * @param[opt] eps: singular values smaller than this will be set to zero
 * @return: true if A is full rank, false if A is rank deficient
 */
bool pinv_QR(const Eigen::MatrixXd &A, Eigen::MatrixXd *invA, double eps = 1e-6);

/*
 * This code is used to compute the equations 6-10 in the main SNS-IK paper.
 * The variable names are taken to match those in the paper, letting k = 1.
 * It seems that the output of this code is dependent on the particular type of QR factorization
 * routine in Eigen, in this case: HouseholderQR.
 *
 * Here are the relationships between the inputs and the outputs, written in Matlab pseudocode
 * (J1*Za0)' = [Ya1, Z1] * [Ra1; 0];  % QR decomposition
 * Jstar = Za0*Ya1*inv(Ra1');  % projected inverse task jacobian
 * Za1 = Za0*Z1;  % updated nullspace projector
 *
 * @param J1:  task jacobian with size [m, n], with m <= n
 * @param Za0:  previous nullSpaceProjector with size [n, n]
 * @param[out] Jstar:  (not named in the paper)
 * @param[out] Za1:  new nullSpaceProjector
 * @param[opt] lambda_max: damping parameter for one of the inverses
 * @param[opt] eps: singular values smaller than this will be set to zero
 */
bool pinv_QR_Z(const Eigen::MatrixXd &J1, const Eigen::MatrixXd &Za0, Eigen::MatrixXd *Jstar,
               Eigen::MatrixXd *Za1, double lambda_max = 1e-6, double eps = 1e-6);

/*
 * This function computes the inverse of the projection of the P matrix onto the dimensions that
 * are selected by W. One way to think about this would be to reorder the dimensions such that
 * the matrix W is block diagonal, with the upper left block the identify matrix and the lower
 * right block all zeros. Then compute the inverse of the block in the reordered version of
 * P, while setting all other entries to zero. The final step is to reorder the dimesions again.
 *
 * This function is related to the concepts in equations 16-17 of the main paper.
 *
 * These concepts are perhaps best expressed in Matlab code:
 *    % This code is optimized for readability - this is a bad numerical implementation
 *    n = 10;  % number of dimensions in P
 *    s = rand(1, n)>0.4;  % selection vector (true == use this joint)
 *    W = diag(s);  % selection matrix (square)
 *    P = randn(length(s));  % random square input matrix
 *    K = W(s, :);  % sub-space selection matrix (rectangular)
 *    M = K*P*K';  % sub-matrix of P, in active joints
 *    A = inv(M);  % inverse of the sub-matrix  % note: avoid doing this!  (prefer linear solve)
 *    B = K'*A*K;  % project B back into the original space
 *    zero = eye(sum(s)) - K*P*B*K';  % condition to test
 *
 *   FIXME: remove direct inverse, replace with linear solve
 *   FIXME: remove hard-coded constants embedded in code
 *   FIXME: is it possible to replace W with a boolean vector?
 *
 * @param W: selection matrix, diagonal entries are either zero or one, other entries are zero.
      W(i, i) == 1 indicates that dimension i should be used, otherwise ignore dimension i
 * @param P: project matrix, must be square. P.rows() == P.cols() == W.cols()
 * @param[out] B: the inverse of the sub-matrix of P that is selected by W
 *      if return false, then B is zeros
 *      M = select(P); // smaller square matrix, size == number of non-zero diagonal entries in W
 *      K = backProject(inverse(M)); // compute the inverse of M and project back to size of P
 *      C = P*K;
 * @return: true iff inversePbar is invertible
 */
bool pinv_forBarP(const Eigen::MatrixXd &W, const Eigen::MatrixXd &P, Eigen::MatrixXd *C);

/*
 * @return true iff the diagonal elements are near unity
 *   FIXME: if A.rows() >= A.cols() causes a failed assertion
 *   FIXME: if A(i,i) > 1.0, then this function returns the wrong answer
 */
bool isIdentity(const Eigen::MatrixXd &A);

/*
 * Compute the pseudo-inverse of a matrix.
 * Note: do not use this to solve a linear system: use solveLinearSystem() instead.
 * @param A: matrix of interest
 * @param eps: small parameter, used for singular value threshold and damping
 * @param[out] invA: pseudo-inverse of the matrix
 * @param[out, opt] rank: the rank of matrix A
 * @param[out, opt] damped: true if a damped pseudo-inverse was used
 * @return: true iff successful
 */
bool pseudoInverse(const Eigen::MatrixXd& A, double eps, Eigen::MatrixXd* invA,
                   int* rank = nullptr, bool* damped = nullptr);

/*
 * Compute the solution to the linear system: A*x = b
 * @param A: matrix of size [n, m]
 * @param b: matrix of size [n, p]
 * @param[out] x: matrix of size [m, p]
 * @param[out, opt] rank: rank of the A
 * @param[out, opt] err: residual error in solution: (A*x-b).squaredNorm()
 * @return: true iff successful
 */
bool solveLinearSystem(const Eigen::MatrixXd& A, const Eigen::MatrixXd& b,
                       Eigen::MatrixXd* x,
                       int* rank = nullptr, double* err = nullptr);

}  // namespace sns_ik

#endif
