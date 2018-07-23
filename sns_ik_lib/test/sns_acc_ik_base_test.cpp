/**  @file sns_acc_ik_base_test.cpp
 *
 *  @brief Unit Test: sns_acc_ik_base solver
 *  @author Matthew Kelly
 *  @author Andy Park
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

#include <gtest/gtest.h>

#include <Eigen/Dense>
#include <ros/console.h>

#include <sns_ik/sns_acc_ik_base.hpp>
#include <sns_ik/sns_ik_base.hpp>
#include "rng_utilities.hpp"
#include "test_utilities.hpp"

/*************************************************************************************************/

/*
 * This test is for the SnsAccIkBase::solve() without any joint limits.
 * This checks whether the solution obtained from SNS IK without limits is valid.
 */
TEST(sns_acc_ik_base, basic_no_limits)
{
  sns_ik::rng_util::setRngSeed(75716, 11487);  // set the initial seed for the random number generators
  int nTest = 100;
  double tol = 1e-10;
  for (int iTest = 0; iTest < nTest; iTest++) {

    // generate a test problem
    int nJoint = sns_ik::rng_util::getRngInt(0, 1, 8);
    int nTask = sns_ik::rng_util::getRngInt(0, 1, nJoint);
    Eigen::MatrixXd J = sns_ik::rng_util::getRngMatrixXd(0, nTask, nJoint, -1.0, 1.0);
    Eigen::MatrixXd dJdq = sns_ik::rng_util::getRngVectorXd(0, nTask, -1.0, 1.0);
    Eigen::VectorXd ddx = sns_ik::rng_util::getRngVectorXd(0, nTask, -1.0, 1.0);
    Eigen::VectorXd tolVec = tol * Eigen::VectorXd::Ones(nJoint);
    Eigen::VectorXd ddq;
    double taskScale;

    // solve
    sns_ik::SnsAccIkBase::uPtr ikSolver = sns_ik::SnsAccIkBase::create(nJoint);
    ASSERT_TRUE(ikSolver.get() != nullptr);
    sns_ik::SnsIkBase::ExitCode exitCode = ikSolver->solve(J, dJdq, ddx, &ddq, &taskScale);
    ASSERT_TRUE(exitCode == sns_ik::SnsIkBase::ExitCode::Success);

    // check requirements
    ASSERT_LE(taskScale, 1.0 + tol);
    ASSERT_GT(taskScale, 0.0);
    sns_ik::test_util::checkEqualVector(taskScale * ddx, J * ddq + dJdq, tol);
  }
}

/*************************************************************************************************/

/*
 * This test is for the SnsAccIkBase::solve() with joint limits.
 * This checks if the solution obtained from SNS IK with limits is valid.
 * Test data including task scale factor is randomly generated and they
 * will be compared against the solution from the solver. A solution is
 * considered valid if it is within the limits while meeting the goal
 * task acceleration with the computed task scale factor.
 */
TEST(sns_acc_ik_base, basic_with_limits)
{
  sns_ik::rng_util::setRngSeed(65444, 48535);  // set the initial seed for the random number generators
  int nTest = 10000;
  double tol = 1e-10;
  int nPass = 0;
  int nFail = 0;
  int nSubOpt = 0;
  double meanSolveTime = 0.0;
  for (int iTest = 0; iTest < nTest; iTest++) {
    // generate a test problem
    int nTask = sns_ik::rng_util::getRngInt(0, 1, 6);
    int nJoint = sns_ik::rng_util::getRngInt(0, nTask, nTask + 4);
    Eigen::MatrixXd J = sns_ik::rng_util::getRngMatrixXd(0, nTask, nJoint, -3.0, 3.0);
    Eigen::ArrayXd ddqLow = sns_ik::rng_util::getRngVectorXd(0, nJoint, -3.0, -1.0);
    Eigen::ArrayXd ddqUpp = sns_ik::rng_util::getRngVectorXd(0, nJoint, 1.0, 3.0);
    Eigen::VectorXd ddqTest = sns_ik::rng_util::getRngArrBndXd(0, ddqLow, ddqUpp).matrix();
    Eigen::VectorXd dJdq = sns_ik::rng_util::getRngVectorXd(0, nTask, -1.0, 1.0);

    // create a task that is feasible with scaling
    Eigen::VectorXd ddxFeas = J*ddqTest + dJdq; // this task acceleration is feasible by definition
    double taskScaleMin = sns_ik::rng_util::getRngDouble(0, 0.2, 1.2);
    taskScaleMin = std::min(1.0, taskScaleMin);  // clamp max value to 1.0
    Eigen::VectorXd ddx = ddxFeas / taskScaleMin;

    // solve
    Eigen::VectorXd ddq;
    double taskScale;
    sns_ik::SnsAccIkBase::uPtr ikSolver = sns_ik::SnsAccIkBase::create(ddqLow, ddqUpp);
    ASSERT_TRUE(ikSolver.get() != nullptr);
    ros::Time startTime = ros::Time::now();
    sns_ik::SnsIkBase::ExitCode exitCode = ikSolver->solve(J, dJdq, ddx, &ddq, &taskScale);
    double solveTime = (ros::Time::now() - startTime).toSec();
    meanSolveTime += solveTime;

    if (exitCode == sns_ik::SnsIkBase::ExitCode::Success) {
      nPass++;
      // check requirements
      ASSERT_LE(taskScale, 1.0 + tol);
      if (taskScale < taskScaleMin - tol) nSubOpt++;
      sns_ik::test_util::checkEqualVector(taskScale * ddx, J * ddq + dJdq, tol);
      sns_ik::test_util::checkVectorLimits(ddqLow, ddq, ddqUpp, tol);
    } else {
      nFail++;
      EXPECT_TRUE(false) << "Solver failed  --  infeasible task?";
    }
  }
  meanSolveTime /= static_cast<double>(nPass + nFail);
  ROS_INFO("Pass: %d  --  Fail: %d  --  nSubOpt: %d  --  Mean solve time: %.4f ms",
           nPass, nFail, nSubOpt, meanSolveTime*1000.0);
}

/*************************************************************************************************/

/*
 * This test is for the SnsAccIkBase::solve() with a secondary goal (as well as joint limits).
 * This checks if the solution obtained from SNS IK with limits is valid.
 * Test data including task scale factor is randomly generated and they
 * will be compared against the solution from the solver. A solution is
 * considered valid if it is within the limits while meeting the goal
 * task acceleration with the computed task scale factor.
 */
TEST(sns_acc_ik_base, basic_with_secondary_goal)
{
  sns_ik::rng_util::setRngSeed(65444, 48535);  // set the initial seed for the random number generators
  int nTest = 10000;
  double tol = 1e-10;
  int nPass = 0;
  int nFail = 0;
  int nSubOpt = 0;
  double meanSolveTime = 0.0;
  double meanTaskScaleCS = 0.0;
  for (int iTest = 0; iTest < nTest; iTest++) {
    // generate a test problem
    int nTask = sns_ik::rng_util::getRngInt(0, 1, 6);
    int nJoint = sns_ik::rng_util::getRngInt(0, nTask, nTask + 4);
    Eigen::MatrixXd J = sns_ik::rng_util::getRngMatrixXd(0, nTask, nJoint, -3.0, 3.0);
    Eigen::ArrayXd ddqLow = sns_ik::rng_util::getRngVectorXd(0, nJoint, -3.0, -1.0);
    Eigen::ArrayXd ddqUpp = sns_ik::rng_util::getRngVectorXd(0, nJoint, 1.0, 3.0);
    Eigen::VectorXd ddqTest = sns_ik::rng_util::getRngArrBndXd(0, ddqLow, ddqUpp).matrix();
    Eigen::VectorXd dJdq = sns_ik::rng_util::getRngVectorXd(0, nTask, -1.0, 1.0);

    // create a task that is feasible with scaling
    Eigen::VectorXd ddxFeas = J*ddqTest + dJdq; // this task acceleration is feasible by definition
    double taskScaleMin = sns_ik::rng_util::getRngDouble(0, 0.2, 1.2);
    taskScaleMin = std::min(1.0, taskScaleMin);  // clamp max value to 1.0
    Eigen::VectorXd ddx = ddxFeas / taskScaleMin;

    // generate a configuration space task as secondary goal
    Eigen::VectorXd ddqCS = sns_ik::rng_util::getRngArrBndXd(0, ddqLow, ddqUpp).matrix();

    // solve
    Eigen::VectorXd ddq;
    double taskScale, taskScaleCS;
    sns_ik::SnsAccIkBase::uPtr ikSolver = sns_ik::SnsAccIkBase::create(ddqLow, ddqUpp);
    ASSERT_TRUE(ikSolver.get() != nullptr);
    ros::Time startTime = ros::Time::now();
    sns_ik::SnsIkBase::ExitCode exitCode = ikSolver->solve(J, dJdq, ddx, ddqCS, &ddq, &taskScale, &taskScaleCS);
    double solveTime = (ros::Time::now() - startTime).toSec();
    meanSolveTime += solveTime;
    meanTaskScaleCS += taskScaleCS;

    if (exitCode == sns_ik::SnsIkBase::ExitCode::Success) {
      nPass++;
      // check requirements
      ASSERT_LE(taskScale, 1.0 + tol);
      if (taskScale < taskScaleMin - tol) nSubOpt++;
      sns_ik::test_util::checkEqualVector(taskScale * ddx, J * ddq + dJdq, tol);
      sns_ik::test_util::checkVectorLimits(ddqLow, ddq, ddqUpp, tol);
    } else {
      nFail++;
      EXPECT_TRUE(false) << "Solver failed  --  infeasible task?";
    }
  }
  meanSolveTime /= static_cast<double>(nPass + nFail);
  meanTaskScaleCS /= static_cast<double>(nPass);
  ROS_INFO("Pass: %d  --  Fail: %d  --  nSubOpt: %d  --  Mean solve time: %.4f ms -- Mean secondary goal scale : %.4f",
           nPass, nFail, nSubOpt, meanSolveTime*1000.0, meanTaskScaleCS);
}

/*************************************************************************************************/

// Run all the tests that were declared with TEST()
int main(int argc, char **argv){
  ros::Time::init();
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
