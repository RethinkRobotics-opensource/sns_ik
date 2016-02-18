/*! \file sns_ikl_math_utils.hpp
 * \brief Math utilities for the SNS IK solvers
 * \author Fabrizio Flacco
 * \author Forrest Rogers-Marcovitz
 */
/*
 *    Copyright 2016 Rethink Robotics
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

using namespace Eigen;

namespace sns_ikl {

#define _USE_DOUBLE_  
/*! \def _USE_DOUBLE_  
 * use values with double precision.
 * This is used to have the possibility to switch (at compile time) between float and double precision in order to be faster or more accurate.
 * Inside the \b IKL we will use the follow definition:
 * - \b MatrixD  a dynamic sized Matrix of real numbers
 * - \b VectorD a column vector of real numbers with dynamic length 
 * - \b Scalar a real number
 */

#ifdef _USE_DOUBLE_
  /*! \typedef MatrixD
   *  dynamic sized matrix of \em float values
   *  (\em double if \b _USE_DOUBLE_ is defined)
   */
  typedef Eigen::Matrix<double, Dynamic, Dynamic> MatrixD;
  /*! \typedef VectorD
   *  dynamic sized column vector of \em float values
   *  (\em double if \b _USE_DOUBLE_ is defined)
   */
  typedef Eigen::Matrix<double, Dynamic, 1> VectorD;
  /*! \typedef Scalar
   *  a \em float value
   *  (\em double if \b _USE_DOUBLE_ is defined)
   */
  typedef double Scalar;
  #define EPS 0.15
  #define LAMBDA_MAX 0.3
  #define EPSQ 1e-10
#else
  typedef Eigen::Matrix<float,Dynamic,Dynamic> MatrixD;
  typedef Eigen::Matrix<float,Dynamic,1> VectorD;
  typedef float Scalar;
  #define EPS 0.15
  #define LAMBDA_MAX 0.3
  #define EPSQ 1e-15
#endif

#define INF 1e20;  // Should this be the std INF?

// compute the pseudoinverse of A: it return 0 if A is (row) rank deficient
bool pinv(const MatrixD &A, MatrixD *invA, Scalar eps = EPS);
// compute the pseudoinverse of A: it return 0 if A is (row) rank deficient
// return also the null space projector P=(P-pinv(A)A)
bool pinv_P(const MatrixD &A, MatrixD *invA, MatrixD *P, Scalar eps = EPS);
// compute the pseudoinverse of A: it return 0 if A is (row) rank deficient
bool pinv_damped(const MatrixD &A, MatrixD *invA, Scalar lambda_max = LAMBDA_MAX,
                 Scalar eps = EPS);
bool pinv_damped_P(const MatrixD &A, MatrixD *invA, MatrixD *P,
                   Scalar lambda_max = LAMBDA_MAX, Scalar eps = EPS);
bool pinv_forBarP(const MatrixD &W, const MatrixD &P, MatrixD *inv);
bool pinv_QR_Z(const MatrixD &A, const MatrixD &Z0, MatrixD *invA, MatrixD *Z,
               Scalar lambda_max = LAMBDA_MAX, Scalar eps = EPS);
bool pinv_QR(const MatrixD &A, MatrixD *invA, Scalar eps = EPS);

bool isIdentity(const MatrixD &A);

}  // namespace sns_ikl

#endif
