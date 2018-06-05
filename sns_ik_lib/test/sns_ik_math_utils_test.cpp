/*! \file sns_ik_math_utils_test.cpp
 * \brief Unit Test: sns_ik_math_utils
 * \author Matthew Kelly
 */
/*
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

#include <gtest/gtest.h>

#include <Eigen/Dense>
#include <ros/console.h>

#include "rng_utilities.hpp"
#include "sns_ik_math_utils.hpp"

/*************************************************************************************************
 *                               Utilities Functions                                             *
 *************************************************************************************************/

/*
 * Element-wise check that two matrices are equal to within some tolerance
 * @param A: first input matrix
 * @param B: second input matrix
 * @param tol: tolerance on equality checks
 */
void checkEqualMatrices(const Eigen::MatrixXd& A, const Eigen::MatrixXd& B, double tol)
{
  ASSERT_EQ(A.rows(), B.rows());
  ASSERT_EQ(A.cols(), B.cols());
  int n = A.rows();
  int m = A.cols();
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < m; j++) {
      ASSERT_NEAR(A(i,j), B(i,j), tol);
    }
  }
}

/*************************************************************************************************/

/*
 * Check that two matrices are related by a Moore-Penrose pseudoinverse, to within tol.
 * X = pinv(A)
 * @param A: matrix A
 * @param X: the pseudoinverse of matrix A
 * @param tol: tolerance on element-wise matrix equality check
 */
void checkPseudoInverse(const Eigen::MatrixXd& A, const Eigen::MatrixXd& X, double tol)
{
  checkEqualMatrices(A*X*A, A, tol);
  checkEqualMatrices(X*A*X, X, tol);
  Eigen::MatrixXd XA = X*A;
  Eigen::MatrixXd AX = A*X;
  checkEqualMatrices(XA, XA.transpose(), tol);
  checkEqualMatrices(AX, AX.transpose(), tol);
}

/*************************************************************************************************/

/*
 * Compute the matrices that form the block-matrix QR decomposition of A = [Y, Z]*[R; 0]
 * @param A: matrix of size [n x m] to perform decomposition of
 * @param[out] Y: left block of Q in the QR decomposition of A
 * @param[out] Z: right block of Q in the QR decomposition of A
 * @param[out] R: upper (non-zero) block of R in the QR decomposition of A
 */
void qrBlockDecompose(const Eigen::MatrixXd& A,
                      Eigen::MatrixXd* Y, Eigen::MatrixXd* Z, Eigen::MatrixXd* R)
{
  int n = A.rows();
  int m = A.cols();
  Eigen::HouseholderQR<Eigen::MatrixXd> qr(A);  // Tests fail with other QR decompositions
  *Y = ((Eigen::MatrixXd) qr.householderQ()).leftCols(m);
  *Z = ((Eigen::MatrixXd) qr.householderQ()).rightCols(n-m);
  *R = qr.matrixQR().topLeftCorner(m, m).template triangularView<Eigen::Upper>();
}


/*************************************************************************************************/

/*
 * Check that the solution to the linear system: A*x = b is valid
 * @param A: matrix of size [n, m]
 * @param b: matrix of size [n, p]
 * @param x: matrix of size [m, p]
 * @param rank: rank of the A
 * @param err: residual error in solution: (A*x-b).squaredNorm()
 * @param tol: tolerance for equality checks
 */
void checkLinearSolve(const Eigen::MatrixXd& A, const Eigen::MatrixXd& b,
                     const Eigen::MatrixXd& x, int rank, double err, double tol)
{
  // check sizes:
  ASSERT_EQ(A.cols(), x.rows());
  ASSERT_EQ(A.rows(), b.rows());
  ASSERT_EQ(x.cols(), b.cols());

  // Check residual:
  ASSERT_NEAR(err, (A*x-b).squaredNorm(), tol);

  // Check linear solve:
  int n = A.rows();
  int m = A.cols();
  ASSERT_LE(rank, std::min(n, m));  // upper bound of rank is based on matrix size
  if (m >= n && rank == std::min(n, m)) {  // there is a feasible solution
    checkEqualMatrices(A*x, b, tol);
  }  //  else: the solution is infeasible, x should minimize the norm
}

/*************************************************************************************************
 *                                        Tests                                                  *
 *************************************************************************************************/

/*
 * Unit test for the method sns_ik::pinv()
 * TODO: pass arbitrary size matrices into pinv() once the code is fixed to support it
 */
TEST(sns_ik_math_utils, pinv_test)
{
  // Run standard tests
  double tolSvd = 1e-8;  // svd minimum value tolerance
  double tolMat = 1e-6;  // tolerance for matrix equality check
  Eigen::MatrixXd X;  // pinv(A)
  int seed = 97460;
  double low = -2.0;  double upp = 2.0;  // bounds on values in the A matrix
  for (int i = 0; i < 25; i++) {
    // generate test data
    seed++;
    int nCol = sns_ik::rng_util::getRngInt(seed + 50322, 3, 9);
    int nRow = sns_ik::rng_util::getRngInt(seed + 87288, 2, (nCol-1));  // TODO: any number of rows
    // Decide whether to use a full-rank or rank degenerate case
    bool fullRank = sns_ik::rng_util::getRngBool(seed + 96717, 0.85);  // usually full rank
    if (fullRank) {  // generate and test for a full-rank input
      Eigen::MatrixXd A = sns_ik::rng_util::getRngMatrixXd(seed + 81068, nRow, nCol, low, upp);
      ASSERT_TRUE(sns_ik::pinv(A, &X, tolSvd));
      checkPseudoInverse(A, X, tolMat);
    } else {  // generate an test for a rank-deficient input
      int nRank = std::min(nRow, nCol) - 1;
      Eigen::MatrixXd A = sns_ik::rng_util::getRngMatrixXdRanked(seed + 81068, nRow, nCol, nRank);
      ASSERT_FALSE(sns_ik::pinv(A, &X, tolSvd));
    }
  }
}

/*************************************************************************************************/

/*
 * Unit test for pinv_P()
 * TODO: pass arbitrary size matrices into pinv() once the code is fixed to support it
 */
TEST(sns_ik_math_utils, pinv_P_test)
{
  // Run standard tests
  double tolSvd = 1e-8;  // svd minimum value tolerance
  double tolMat = 1e-6;  // tolerance for matrix equality check
  Eigen::MatrixXd X;  // pinv(A)
  int seed = 84658;
  double low = -2.0;  double upp = 2.0;  // bounds on values in the A and P matrix
  Eigen::MatrixXd A;
  for (int iTest = 0; iTest < 25; iTest++) {
    // generate test data
    seed++;
    int nCol = sns_ik::rng_util::getRngInt(seed + 80226, 3, 9);
    int nRow = sns_ik::rng_util::getRngInt(seed + 49166, 2, (nCol-1));  // TODO: any number of rows
    Eigen::MatrixXd P = sns_ik::rng_util::getRngMatrixXd(seed + 70730, nCol, nCol, low, upp);
    Eigen::MatrixXd PP = P;  // Projected P matrix
    // Decide whether to use a full-rank or rank degenerate case
    bool fullRank = sns_ik::rng_util::getRngBool(seed + 11052, 0.85);  // usually full rank
    if (fullRank) {  // generate and test for a full-rank input
      A = sns_ik::rng_util::getRngMatrixXd(seed + 75011, nRow, nCol, low, upp);
      ASSERT_TRUE(sns_ik::pinv_P(A, &X, &PP, tolSvd));
      checkPseudoInverse(A, X, tolMat);
      checkEqualMatrices(PP, P - X*A, tolMat);
    } else {  // generate an test for a rank-deficient input
      int nRank = std::min(nRow, nCol) - 1;
      A = sns_ik::rng_util::getRngMatrixXdRanked(seed + 91579, nRow, nCol, nRank);
      ASSERT_FALSE(sns_ik::pinv_P(A, &X, &PP, tolSvd));
    }
  }
}

/*************************************************************************************************/

/*
 * Unit test for pinv_damped_P()
 * TODO: pass arbitrary size matrices into pinv_damped_P() once the code is fixed to support it
 */
TEST(sns_ik_math_utils, pinv_damped_P_test)
{
  // Run standard tests
  double tolSvd = 1e-6;  // svd minimum value tolerance
  double tolMat = 1e-5;  // tolerance for matrix equality check
  double lambda = 1e-6;  // damping parameter for damped pseudoinverse
  Eigen::MatrixXd X;  // pinv(A)
  int seed = 81012;
  double low = -2.0;  double upp = 2.0;  // bounds on values in the A and P matrix
  Eigen::MatrixXd A;
  for (int iTest = 0; iTest < 25; iTest++) {
    // generate test data
    int nCol = sns_ik::rng_util::getRngInt(seed + 82956, 3, 9);
    int nRow = sns_ik::rng_util::getRngInt(seed + 44438, 2, (nCol-1));  // TODO: any number of rows
    Eigen::MatrixXd P = sns_ik::rng_util::getRngMatrixXd(seed + 54852, nCol, nCol, low, upp);
    Eigen::MatrixXd PP = P;  // Projected P matrix
    // Decide whether to use a full-rank or rank degenerate case
    bool fullRank = sns_ik::rng_util::getRngBool(seed + 86532, 0.85);  // usually full rank
    if (fullRank) {  // generate and test for a full-rank input
      A = sns_ik::rng_util::getRngMatrixXd(seed + 75312, nRow, nCol, low, upp);
      ASSERT_TRUE(sns_ik::pinv_damped_P(A, &X, &PP, lambda, tolSvd));
      checkPseudoInverse(A, X, tolMat);
    } else {  // generate an test for a rank-deficient input
      int nRank = std::min(nRow, nCol) - 1;
      A = sns_ik::rng_util::getRngMatrixXdRanked(seed + 88433, nRow, nCol, nRank);
      ASSERT_FALSE(sns_ik::pinv_damped_P(A, &X, &PP, lambda, tolSvd));
    }
    checkEqualMatrices(PP, P - X*A, tolMat);
  }
}

/*************************************************************************************************/

TEST(sns_ik_math_utils, pinv_forBarP_test)
{
  double tolMat = 1e-6;  // tolerance for matrix equality check
  int seed = 36230;
  double low = -2.0;  double upp = 2.0;  // bounds on values in the A and P matrix
  Eigen::MatrixXd W, C, K, P;
  for (int iTest = 0; iTest < 40; iTest++) {
    seed++;

    // Compute the selection matrix (W) and data matrix (P)
    int nDim = sns_ik::rng_util::getRngInt(seed + 85004, 2, 11);
    std::vector<bool> selectionVector(nDim);
    for (int iDim = 0; iDim < nDim; iDim++) {
      selectionVector[iDim] = sns_ik::rng_util::getRngBool(0, 0.6);
    }
    int nnz = 0; // number of non-zero elements
    for (int iDim = 0; iDim < nDim; iDim++) {
      if (selectionVector[iDim]) { nnz++; }
    }
    if (nnz == 0) { selectionVector[0] = true;  nnz = 1; }
    W = Eigen::MatrixXd::Zero(nDim, nDim);  // selection matrix
    for (int iDim = 0; iDim < nDim; iDim++) {
      if (selectionVector[iDim]) { W(iDim, iDim) = 1.0; }
    }

    // Decide whether to use a full-rank or rank degenerate case
    bool fullRank = sns_ik::rng_util::getRngBool(seed + 86532, 0.7);  // usually full rank
    bool smallNullSpace = nDim - nnz < 2;  // true if the nullspace has fewer than two dimensions
    if (fullRank || smallNullSpace) {  // generate and test for a full-rank input
      P = sns_ik::rng_util::getRngMatrixXd(seed + 47524, nDim, nDim, low, upp);

      // Compute the sub-selection matrix:
      K = Eigen::MatrixXd::Zero(nnz, nDim);
      int idx = 0;
      for (int iDim = 0; iDim < nDim; iDim++) {
        if (selectionVector[iDim]) { K(idx, iDim) = 1.0; idx++; }
      }

      // Call the function to be tested:
      ASSERT_TRUE(sns_ik::pinv_forBarP(W, P, &C));
      // Check that the result is valid:
      checkEqualMatrices(Eigen::MatrixXd::Identity(nnz, nnz), K*C*(K.transpose()), tolMat);

    } else { // use the rank-deficient case
        int nRank = 1;
        P = sns_ik::rng_util::getRngMatrixXdRanked(seed + 20880, nDim, nDim, nRank);

        // Call the function to be tested:
        ASSERT_FALSE(sns_ik::pinv_forBarP(W, P, &C));
        // Check that the result is valid:
        checkEqualMatrices(Eigen::MatrixXd::Zero(nDim, nDim), C, tolMat);
    }
  }
}

/*************************************************************************************************/

/*
 * Unit test for the method sns_ik::pinv_QR()
 * TODO: pass arbitrary size matrices into pinv_QR() once the code is fixed to support it
 */
TEST(sns_ik_math_utils, pinv_QR_test)
{
  // Run standard tests
  double tolSvd = 1e-8;  // svd minimum value tolerance
  double tolMat = 1e-6;  // tolerance for matrix equality check
  Eigen::MatrixXd X;  // pinv(A)
  int seed = 64799;
  double low = -2.0;  double upp = 2.0;  // bounds on values in the A matrix
  for (int iTest = 0; iTest < 25; iTest++) {
    // generate test data
    seed++;
    int nCol = sns_ik::rng_util::getRngInt(seed + 59010, 3, 9);
    int nRow = sns_ik::rng_util::getRngInt(seed + 73052, 2, (nCol-1));  // TODO: any number of rows
    // Decide whether to use a full-rank or rank degenerate case
    bool fullRank = sns_ik::rng_util::getRngBool(seed + 51348, 0.85);  // usually full rank
    if (fullRank) {  // generate and test for a full-rank input
      Eigen::MatrixXd A = sns_ik::rng_util::getRngMatrixXd(seed + 23253, nRow, nCol, low, upp);
      ASSERT_TRUE(sns_ik::pinv_QR(A, &X, tolSvd));
      checkPseudoInverse(A, X, tolMat);
    } else {  // generate an test for a rank-deficient input
      int nRank = std::min(nRow, nCol) - 1;
      Eigen::MatrixXd A = sns_ik::rng_util::getRngMatrixXdRanked(seed + 19669, nRow, nCol, nRank);
      ASSERT_FALSE(sns_ik::pinv_QR(A, &X, tolSvd));
    }
  }
}

/*************************************************************************************************/

/*
 * Unit test for pinv_QR_Z()
 */
TEST(sns_ik_math_utils, pinv_QR_Z_test)
{
  // Run standard tests
  double tolSvd = 1e-6;  // svd minimum value tolerance
  double tolMat = 1e-5;  // tolerance for matrix equality check
  double lambda = 1e-6;  // damping parameter for damped pseudoinverse
  Eigen::MatrixXd Jstar, Za1;  // outputs of pinv_QR_Z()
  int seed = 98606;
  double low = -2.0;  double upp = 2.0;  // bounds on values in the A and P matrix
  Eigen::MatrixXd A;
  for (int iTest = 0; iTest < 25; iTest++) {
    // generate test data
    int nJoint = sns_ik::rng_util::getRngInt(seed + 98138, 4, 12);
    int nTask = sns_ik::rng_util::getRngInt(seed + 90906, 2, nJoint);
    bool fullRank = sns_ik::rng_util::getRngBool(seed + 51348, 0.85);  // usually full rank
    Eigen::MatrixXd J1;
    Eigen::MatrixXd Za0 = sns_ik::rng_util::getRngMatrixXd(seed + 10218, nJoint, nJoint, low, upp);
    if (fullRank) {  // generate and test for a full-rank input
      J1 = sns_ik::rng_util::getRngMatrixXd(seed + 65903, nTask, nJoint, low, upp);
      ASSERT_TRUE(sns_ik::pinv_QR_Z(J1, Za0, &Jstar, &Za1, lambda, tolSvd));
    } else {  // test the rank-deficient case
      int nRank = sns_ik::rng_util::getRngInt(seed + 82189, 1, nTask-1);
      Eigen::MatrixXd J1 = sns_ik::rng_util::getRngMatrixXdRanked(seed + 48185, nTask, nJoint, nRank);
      ASSERT_FALSE(sns_ik::pinv_QR_Z(J1, Za0, &Jstar, &Za1, lambda, tolSvd));
    }
    // generate matrices for checking the result
    Eigen::MatrixXd A = (J1 * Za0).transpose();
    Eigen::MatrixXd Ya1, Z1, Ra1;  // QR decomposition block matrices
    qrBlockDecompose(A, &Ya1, &Z1, &Ra1);  // perform QR decomposition (for checking result only)
    // checks:
    checkEqualMatrices(Za1, Za0 * Z1, tolMat);
    checkEqualMatrices(Jstar * (Ra1.transpose()), Za0 * Ya1, tolMat);
  }
}

/*************************************************************************************************/

// Unit test for isIdentity()
TEST(sns_ik_math_utils, isIdentity_test)
{
  Eigen::MatrixXd M = Eigen::MatrixXd::Identity(3, 3);
  ASSERT_TRUE(sns_ik::isIdentity(M));
  M(1,1) = 0.3;
  ASSERT_FALSE(sns_ik::isIdentity(M));
}

/*************************************************************************************************/

/*
 * Unit test for pseudoInverse() with full rank A matrix
 *  -- this is primarily a regression test, confirming that the new implementation of the pseudo-
 *     inverse gives the same results as the version that was originally implemented in this code.
 */
TEST(sns_ik_math_utils, pseudoInverse_fullRank_test)
{
  // test parameters
  double tolSvd = 1e-8;  // svd minimum value tolerance
  double tolMat = 1e-6;  // tolerance for matrix equality check
  int nTest = 20;
  int seed = 84658;
  double low = -2.0;  double upp = 2.0;  // bounds on values in the A matrix
  // test setup
  Eigen::MatrixXd X;  // pinv(A)
  Eigen::MatrixXd invA;  // pseudoInverse(A)
  Eigen::MatrixXd A;  // test matrix
  int rank; // rank of invA, computed by pseudoInverse(A)
  bool damped;  // true iff invA was computed with a damped pseudo-inverse
  for (int iTest = 0; iTest < nTest; iTest++) {
    seed++;
    int nCol = sns_ik::rng_util::getRngInt(seed + 63541, 3, 9);
    int nRow = sns_ik::rng_util::getRngInt(0, 2, (nCol-1));
    A = sns_ik::rng_util::getRngMatrixXd(seed + 68409, nRow, nCol, low, upp);
    ASSERT_TRUE(sns_ik::pinv_damped_P(A, &X, nullptr, tolSvd, tolSvd));
    ASSERT_TRUE(sns_ik::pseudoInverse(A, tolSvd, &invA, &rank, &damped));
    checkPseudoInverse(A, X, tolMat);
    checkPseudoInverse(A, invA, tolMat);
    ASSERT_FALSE(damped);
    ASSERT_EQ(rank, nRow);
    checkEqualMatrices(X, invA, tolMat);
  }
}

/*************************************************************************************************/

/*
 * Unit test for pseudoInverse() for rank-deficient inputs
 *  -- this is primarily a regression test, confirming that the new implementation of the pseudo-
 *     inverse gives the same results as the version that was originally implemented in this code.
 */
TEST(sns_ik_math_utils, pseudoInverse_damped_test)
{
  double tolSvd = 1e-8;  // svd minimum value tolerance
  double tolMat = 1e-6;  // tolerance for matrix equality check
  int nTest = 15;
  int seed = 76407;
  Eigen::MatrixXd X;  // pinv(A)
  Eigen::MatrixXd invA;  // pseudoInverse(A)
  Eigen::MatrixXd A;
  int rank; // rank of invA, computed by pseudoInverse(A)
  bool damped;  // true iff invA was computed with a damped pseudo-inverse
  for (int iTest = 0; iTest < nTest; iTest++) {
    // generate test data
    seed++;
    int nCol = sns_ik::rng_util::getRngInt(seed + 75019, 6, 9);
    int nRow = sns_ik::rng_util::getRngInt(seed + 94846, 4, (nCol-1));
    int nRank = sns_ik::rng_util::getRngInt(seed + 94846, 1, (nRow-1));
    A = sns_ik::rng_util::getRngMatrixXdRanked(seed + 43203, nRow, nCol, nRank);
    ASSERT_FALSE(sns_ik::pinv_damped_P(A, &X, nullptr, tolSvd, tolSvd));
    ASSERT_TRUE(sns_ik::pseudoInverse(A, tolSvd, &invA, &rank, &damped));
    checkPseudoInverse(A, X, tolMat);
    checkPseudoInverse(A, invA, tolMat);
    ASSERT_TRUE(damped);
    ASSERT_EQ(rank, nRank);
    checkEqualMatrices(X, invA, tolMat);
  }
}

/*************************************************************************************************/

TEST(sns_ik_math_utils, solveLinearSystem_fullRank_test)
{
  // test parameters
  double tol = 1e-10;  // tolerance for matrix equality check
  int nTest = 25;
  int s1 = 84658;  // seed for the integer RNG
  int s2 = 63541;  // seed for the floating point RNG
  double low = -2.0;  double upp = 2.0;  // bounds on values in the A matrix
  // test setup
  Eigen::MatrixXd A, b;  // test linear system
  Eigen::MatrixXd x;  // result
  int rank; // rank of A, computed by linear solver
  double err;  // residual error in the solution
  for (int iTest = 0; iTest < nTest; iTest++) {
    int n = sns_ik::rng_util::getRngInt(s1, 3, 9);
    int m = sns_ik::rng_util::getRngInt(0, 3, 9);
    int p = sns_ik::rng_util::getRngInt(0, 1, 3);
    A = sns_ik::rng_util::getRngMatrixXd(s2, n, m, low, upp);
    b = sns_ik::rng_util::getRngMatrixXd(s2, n, p, low, upp);
    ASSERT_TRUE(sns_ik::solveLinearSystem(A, b, &x, &rank, &err));
    checkLinearSolve(A, b, x, rank, err, tol);
    s1 = s2 = 0; // let the RNG automatically increment after the first iteration
  }
}

/*************************************************************************************************/

TEST(sns_ik_math_utils, solveLinearSystem_rankDeficient_test)
{
  // test parameters
  double tol = 1e-10;  // tolerance for matrix equality check
  int nTest = 20;
  int s1 = 84658;  // seed for the integer RNG
  int s2 = 63541;  // seed for the floating point RNG
  double low = -2.0;  double upp = 2.0;  // bounds on values in the A matrix
  // test setup
  Eigen::MatrixXd A, b;  // test linear system
  Eigen::MatrixXd x;  // result
  int rank; // rank of A, computed by linear solver
  double err;  // residual error in the solution
  for (int iTest = 0; iTest < nTest; iTest++) {
    int n = sns_ik::rng_util::getRngInt(s1, 4, 9);  // number of equations
    int m = sns_ik::rng_util::getRngInt(0, 4, 9); // number of decision variables
    int p = sns_ik::rng_util::getRngInt(0, 1, 3);
    int r = sns_ik::rng_util::getRngInt(0, 1, std::min(n, m));  // rank
    A = sns_ik::rng_util::getRngMatrixXdRanked(s2, n, m, r);
    b = sns_ik::rng_util::getRngMatrixXd(s2, n, p, low, upp);
    ASSERT_TRUE(sns_ik::solveLinearSystem(A, b, &x, &rank, &err));
    ASSERT_EQ(r, rank);  // check the rank
    checkLinearSolve(A, b, x, rank, err, tol);
    s1 = s2 = 0; // let the RNG automatically increment after the first iteration
  }
}


// Run all the tests that were declared with TEST()
int main(int argc, char **argv){
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
