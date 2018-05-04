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

namespace sns_ik {

namespace {
  const double EPS = 1e-6;
  const double LAMBDA_MAX = 1e-6; //0.3
}

static const double INF = std::numeric_limits<double>::max();

// compute the pseudoinverse of A: it return 0 if A is (row) rank deficient
bool pinv(const Eigen::MatrixXd &A, Eigen::MatrixXd *invA, double eps = EPS);
// compute the pseudoinverse of A: it return 0 if A is (row) rank deficient
// return also the null space projector P=(P-pinv(A)A)
bool pinv_P(const Eigen::MatrixXd &A, Eigen::MatrixXd *invA, Eigen::MatrixXd *P, double eps = EPS);
// compute the pseudoinverse of A: it return 0 if A is (row) rank deficient
bool pinv_damped(const Eigen::MatrixXd &A, Eigen::MatrixXd *invA, double lambda_max = LAMBDA_MAX,
                 double eps = EPS);
bool pinv_damped_P(const Eigen::MatrixXd &A, Eigen::MatrixXd *invA, Eigen::MatrixXd *P,
                   double lambda_max = LAMBDA_MAX, double eps = EPS);
bool pinv_forBarP(const Eigen::MatrixXd &W, const Eigen::MatrixXd &P, Eigen::MatrixXd *inv);
bool pinv_QR_Z(const Eigen::MatrixXd &A, const Eigen::MatrixXd &Z0, Eigen::MatrixXd *invA,
               Eigen::MatrixXd *Z, double lambda_max = LAMBDA_MAX, double eps = EPS);
bool pinv_QR(const Eigen::MatrixXd &A, Eigen::MatrixXd *invA, double eps = EPS);

bool isIdentity(const Eigen::MatrixXd &A);

}  // namespace sns_ikl

#endif
