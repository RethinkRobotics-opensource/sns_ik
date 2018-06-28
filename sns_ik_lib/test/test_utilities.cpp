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

#include "test_utilities.hpp"

namespace sns_ik {
namespace test_util {

/*************************************************************************************************/

void checkEqualVector(const Eigen::VectorXd& A, const Eigen::VectorXd& B, double tol)
{
  ASSERT_EQ(A.size(), B.size());
  int n = A.size();
  for (int i = 0; i < n; i++) {
    ASSERT_NEAR(A(i), B(i), tol);
  }
}

/*************************************************************************************************/

void checkVectorLimits(const Eigen::VectorXd& low, const Eigen::VectorXd& val, const Eigen::VectorXd& upp, double tol)
{
  ASSERT_EQ(low.size(), val.size());
  ASSERT_EQ(upp.size(), val.size());
  int n = val.size();
  for (int i = 0; i < n; i++) {
    ASSERT_LE(val(i), upp(i) + tol);
    ASSERT_GE(val(i), low(i) - tol);
  }
}

/*************************************************************************************************/

} // namespace test_util
} // namespace sns_ik