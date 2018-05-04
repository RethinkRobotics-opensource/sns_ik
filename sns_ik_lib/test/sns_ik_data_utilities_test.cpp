/*! \file fosns_velocity_ik.hpp
 * \brief Unit Test: sns_ik_math_utils
 * \author Matthew Kelly
 *
 * Unit tests for the sns_ik_data_utilities suite of functions.
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

#include <sns_ik/sns_ik_data_utilities.hpp>

/*************************************************************************************************/

/*
 * Unit test for the method sns_ik::data_util::getRngDouble_test()
 */
TEST(sns_ik_data_utils, getRngDouble_test)
{
  // Check repeatability
  double tol = 1e-12;
  EXPECT_NEAR(sns_ik::data_util::getRngDouble(15381, -3, 8), 6.0196095670444016, tol);
  EXPECT_NEAR(sns_ik::data_util::getRngDouble(22247, -3, 8), 2.8762760491777524, tol);
  EXPECT_NEAR(sns_ik::data_util::getRngDouble(75855, -3, 8), -2.7505936842008896, tol);
  EXPECT_NEAR(sns_ik::data_util::getRngDouble(83586, -3, 8), 4.3001619259988875, tol);
  EXPECT_NEAR(sns_ik::data_util::getRngDouble(0, -3, 8), 0.10624858533393144, tol);
  EXPECT_NEAR(sns_ik::data_util::getRngDouble(0, -3, 8), 4.6314850369164899, tol);
  EXPECT_NEAR(sns_ik::data_util::getRngDouble(0, -3, 8), 6.0602696997106875, tol);

  // Check bounds:
  for (int i = 0; i < 49; i++) {
    double low = sns_ik::data_util::getRngDouble(14469 + i, -5, 20);
    double upp = low + sns_ik::data_util::getRngDouble(25209 + i, 1e-6, 10);
    ASSERT_LT(low, upp);
    double x = sns_ik::data_util::getRngDouble(0, low, upp);
    ASSERT_LE(low, x);
    ASSERT_LE(x, upp);
  }
}

/*************************************************************************************************/

// Unit test for sns_ik::data_util::getRngInt_test()
TEST(sns_ik_data_utils, getRngInt_test)
{
  // Check repeatability
  EXPECT_EQ(sns_ik::data_util::getRngInt(81656, 0, 8), 7);
  EXPECT_EQ(sns_ik::data_util::getRngInt(20001, 1, 1), 1);  // edge case test
  EXPECT_EQ(sns_ik::data_util::getRngInt(87349, 2, 1), 2);  // edge case test
  EXPECT_EQ(sns_ik::data_util::getRngInt(46785, 0, 80), 33);
  EXPECT_EQ(sns_ik::data_util::getRngInt(0, -2, 8), 7);
  EXPECT_EQ(sns_ik::data_util::getRngInt(0, 2, 7), 6);
  EXPECT_EQ(sns_ik::data_util::getRngInt(0, 2, 7), 5);

  // Check bounds:
  for (int i = 0; i < 49; i++) {
    int low = sns_ik::data_util::getRngInt(85987 + i, -5, 20);
    int upp = low + sns_ik::data_util::getRngInt(40075 + i, 0, 10);
    ASSERT_LE(low, upp);
    int x = sns_ik::data_util::getRngInt(0, low, upp);
    ASSERT_LE(low, x);
    ASSERT_LE(x, upp);
  }
}

/*************************************************************************************************/
// Unit test for sns_ik::data_util::getRngMatrixXd
TEST(sns_ik_data_utilities, getRngMatrixXd_test)
{
  // Check repeatability:
  double tol = 1e-12;
  Eigen::MatrixXd X;
  X = sns_ik::data_util::getRngMatrixXd(46215, 3, 9, -0.2, 0.9);
  EXPECT_NEAR(X(2, 5), 0.58648152100546125, tol);
  EXPECT_NEAR(X(1, 4), 0.75799866897218138, tol);
  X = sns_ik::data_util::getRngMatrixXd(0, 3, 6, -0.9, 3.9);
  EXPECT_NEAR(X(2, 3), 0.36170151893839109, tol);
  EXPECT_NEAR(X(1, 2), 1.3053972293469402, tol);
  X = sns_ik::data_util::getRngMatrixXd(0, 3, 6, -0.9, 3.9);
  EXPECT_NEAR(X(2, 3), 1.0321424117752471, tol);
  EXPECT_NEAR(X(1, 2), 1.3353510589525368, tol);
  // Check bounds and size:
  int seed = 26569;
  for (int i = 0; i < 15; i++) {
    // Generate the parameters for the matrix
    seed++;
    int nRows = sns_ik::data_util::getRngInt(seed + 20817, 1, 9);
    int nCols = sns_ik::data_util::getRngInt(seed + 66461, 1, 9);
    int nRank = sns_ik::data_util::getRngInt(seed + 77126, 1, 9);
    double low = sns_ik::data_util::getRngDouble(seed + 30185, -5, 20);
    double upp = low + sns_ik::data_util::getRngDouble(seed + 45177, 1e-6, 10);
    X = sns_ik::data_util::getRngMatrixXd(seed + 10883, nRows, nCols, low, upp);
    ASSERT_EQ(X.rows(), nRows);
    ASSERT_EQ(X.cols(), nCols);
    ASSERT_GE(X.minCoeff(), low);
    ASSERT_LE(X.maxCoeff(), upp);
  }
}

/*************************************************************************************************/

// Unit test for sns_ik::data_util::getRngMatrixXdRanked()
TEST(sns_ik_data_utilities, getRngMatrixXdRanked_test)
{
  int seed = 95314;
  for (int i = 0; i < 15; i++) {

    // Generate the parameters for the matrix
    seed++;
    int nRows = sns_ik::data_util::getRngInt(seed + 94522, 1, 9);
    int nCols = sns_ik::data_util::getRngInt(seed + 67168, 1, 9);
    int nRank = sns_ik::data_util::getRngInt(seed + 86769, 1, 9);

    // Generate the matrix
    Eigen::MatrixXd X = sns_ik::data_util::getRngMatrixXdRanked(seed + 12880, nRows, nCols, nRank);

    // Check the rank:
    nRank = std::min({nRows, nCols, nRank});
    Eigen::ColPivHouseholderQR<Eigen::MatrixXd> decomp(X);
    ASSERT_LE(decomp.rank(), nRank);

    // Check the size:
    ASSERT_EQ(X.rows(), nRows);
    ASSERT_EQ(X.cols(), nCols);
  }
}

/*************************************************************************************************/

// Run all the tests that were declared with TEST()
int main(int argc, char **argv){
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
