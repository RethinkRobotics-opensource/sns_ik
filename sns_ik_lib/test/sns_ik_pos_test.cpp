/*! \file sns_ik_pos_test.cpp
 * \brief Unit Test: sns_ik position solver
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

#include <kdl/chainiksolverpos_nr_jl.hpp>
#include <kdl/chainiksolvervel_pinv.hpp>
#include <Eigen/Dense>
#include <kdl/chain.hpp>
#include <kdl/chainfksolvervel_recursive.hpp>
#include <map>
#include <ros/console.h>
#include <ros/duration.h>
#include <ros/time.h>

#include "rng_utilities.hpp"
#include "sawyer_model.hpp"
#include <sns_ik/sns_ik.hpp>

/*
 * Define a common interface to call both SNS and KDL solvers
 */
typedef std::function<int(const KDL::JntArray&, const KDL::Frame&, KDL::JntArray&)> IkSolver;

/*
 * This file provides a unit and regression test for the top-level SNS-IK solver, along with the
 * various solvers that are derived from it.
 */

/*************************************************************************************************
 *                               Position IK Test Params                                         *
 *************************************************************************************************/

// How many position IK tests should be run?
static const int POS_IK_TEST_COUNT = 250;

// Should the position IK tests print the results (beyond pass/fail) to terminal?
static const bool POS_IK_TEST_VERBOSE = true;

// List of perturbations to apply to the solver initialization:
static const std::vector<double> POS_IK_TEST_DELTA_LIST = {1e-6, 1e-4, 1e-2, 0.5};  // radians

// When comparing perturbations, what score to give for successful solve for each delta
static const std::vector<double> POS_IK_TEST_DELTA_SCORE = {0.05, 0.1, 0.25, 0.6};

// Perturbation to apply to the null-space bias
static const double POS_IK_TEST_BIAS_DELTA = 0.1;  // radians

/*************************************************************************************************
 *                               Utilities Functions                                             *
 *************************************************************************************************/

/*
 * A struct to store the result of an individual test
 */
struct PosTestResult {
  ros::Duration solveTime;  // duration spent in the solver on successful solves
  int nPass;  // number of solves that passed
  int nFail;  // number of solves that failed
  double score; // overall score, given which perturbations solved correctly [0.0, 1.0]
};

/*************************************************************************************************/


/*
 * A function to run a benchmarking test on the position IK solvers.
 *
 * Structure:
 * - generate N test poses in joint space
 * - use forward kinematics to compute end-effector pose
 * - generate seeds for a variety of perturbations, from close to far
 * - measure and report solve time and success rate
 * - require each solver to get the correct solution when initialization is very good
 *
 * @param seed: seed to pass to the RNG on each call
 *              if seed == 0, then RNG seed is not updated
 * @param qLow:  lower joint limits
 * @param qUpp:  upper joint limits
 * @param fwdKin:  forward kinematics solver
 * @param ikSolver:  inverse kinematics solver
 * @return: position IK test result
 */
PosTestResult runPosIkSingleTest(int seed, const KDL::JntArray& qLow, const KDL::JntArray& qUpp,
                                KDL::ChainFkSolverPos_recursive& fwdKin,
                                IkSolver& ikSolver)
{
  // Compute a random joint position for the solution and the corresponding test pose
  KDL::JntArray qTest = sns_ik::rng_util::getRngBoundedJoints(seed, qLow, qUpp);
  KDL::Frame zTest;
  fwdKin.JntToCart(qTest, zTest);

  // compute the initialization and nullspace bias:
  int nDel = POS_IK_TEST_DELTA_LIST.size();

  std::vector<KDL::JntArray> initList(nDel);
  std::vector<KDL::JntArray> biasList(nDel);
  for (int iDel = 0; iDel < nDel; iDel++) {
    initList[iDel] = sns_ik::rng_util::getNearbyJoints(0, qTest, POS_IK_TEST_DELTA_LIST[iDel], qLow, qUpp);
    biasList[iDel] = sns_ik::rng_util::getNearbyJoints(0, qTest, POS_IK_TEST_BIAS_DELTA, qLow, qUpp);
  }

  // loop over each solver type
  KDL::JntArray qSoln;
  PosTestResult result;
  result.solveTime = ros::Duration(0);
  result.nPass = 0;
  result.nFail = 0;
  result.score = 0.0;

  // run benchmarking for the set of nullspace biases and joint seeds
  for (int iDel = 0; iDel < nDel; iDel++) {
    ros::Time startTime = ros::Time::now();
    KDL::JntArray init = initList[iDel];
    KDL::JntArray bias = biasList[iDel];

   int exitFlag = ikSolver(init, zTest, qSoln);

    if (exitFlag >= 0) {
      result.solveTime += ros::Time::now() - startTime;
      result.nPass++;
      result.score += POS_IK_TEST_DELTA_SCORE[iDel];
    } else {
      result.nFail++;
    }
  }
  return result;
}

/*************************************************************************************************/

/*
 * @param seed: seed to pass to the RNG on each call
 *              if seed == 0, then seed is ignored
 * @param fwdKin: forward kinematics solver for position
 * @param invKin: inverse kinematics solver for position (to be tested)
 * @param qLow: lower bounds on joints (for generating the test problem)
 * @param qLow: upper bounds on joints (for generating the test problem)
 * @param solverName: solver name, used for logging only
 * @param nTest: number of tests to run
 */
void runGeneralPosIkTest(int seed, KDL::ChainFkSolverPos_recursive& fwdKin, IkSolver& invKin,
                  const KDL::JntArray& qLow, const KDL::JntArray& qUpp, std::string solverName,
                  int nTest = POS_IK_TEST_COUNT)
{
  // Set up the data structure for the results:
  PosTestResult results;  // accumulate results here
  PosTestResult tmp; // individual test output here

  // run the test many times and accumulate the results
  bool success = true;
  for (int iTest = 0; iTest < nTest; iTest++) {
    seed++;
    tmp = runPosIkSingleTest(seed, qLow, qUpp, fwdKin, invKin);
    if (iTest == 0) {
      results = tmp;
    } else {
      results.nPass += tmp.nPass;
      results.nFail += tmp.nFail;
      results.score += tmp.score;
      results.solveTime += tmp.solveTime;
      if (tmp.nPass < 2) { success = false; }
    }
  }

  // Print out the results:
  if (POS_IK_TEST_VERBOSE) {
    int nPass = results.nPass;
    int nFail = results.nFail;
    int nTotal = nPass + nFail;
    double pass = static_cast<double>(nPass);
    double total = static_cast<double>(nTotal);
    ROS_INFO("Position IK Test  -->  pass: %d, fail: %d, score: %f  --  %f ms  --  %s",
              nPass, nFail, results.score / total,
              1000.0 * results.solveTime.toSec() / pass,
              solverName.c_str());
  };

  // Final check:
  ASSERT_TRUE(success);
}

/*************************************************************************************************/

void runSnsPosIkTest(int seed, sns_ik::VelocitySolveType solverType) {
  // Create a sawyer model:
  std::vector<std::string> jointNames;
  KDL::Chain sawyerChain = sns_ik::sawyer_model::getSawyerKdlChain(&jointNames);
  KDL::JntArray qLow, qUpp, vMax, aMax;
  sns_ik::sawyer_model::getSawyerJointLimits(&qLow, &qUpp, &vMax, &aMax);

  // Create a forward-kinematics solver:
  KDL::ChainFkSolverPos_recursive fwdKin(sawyerChain);

  // Create a SNS-IK solver:
  sns_ik::SNS_IK ikSolver(sawyerChain, qLow, qUpp, vMax, aMax, jointNames);
  ikSolver.setVelocitySolveType(solverType);

  // Function template for position IK solver
  IkSolver invKin = [&ikSolver](const KDL::JntArray& qInit, const KDL::Frame& pGoal, KDL::JntArray& qSoln) {
    return ikSolver.CartToJnt(qInit, pGoal, qSoln);
  };

  // Run the test:
  /*
   * Note: This seems to cause a segfault in ubuntu 16.04 with Eigen 3.3.4, but it
   * works fine in Ubuntu 14.04 with Eigen 3.2.0. FIXME
   */
  runGeneralPosIkTest(seed, fwdKin, invKin, qLow, qUpp, sns_ik::toStr(solverType));
}

/*************************************************************************************************
 *                                        Tests                                                  *
 *************************************************************************************************/

/*
 * Run the test on the standard KDL position IK solver.
 */
TEST(sns_ik_pos, KDL_test_1) {
  // Create a sawyer model:
  std::vector<std::string> jointNames;
  KDL::Chain sawyerChain = sns_ik::sawyer_model::getSawyerKdlChain(&jointNames);
  KDL::JntArray qLow, qUpp, vMax, aMax;
  sns_ik::sawyer_model::getSawyerJointLimits(&qLow, &qUpp, &vMax, &aMax);

  // Create a forward-kinematics solver:
  KDL::ChainFkSolverPos_recursive fwdKin(sawyerChain);

  // Create a SNS-IK solver:
  // sns_ik::SNS_IK ikSolver(sawyerChain, qLow, qUpp, vMax, aMax, jointNames);
  // ikSolver.setVelocitySolveType(velSolver);
  KDL::ChainIkSolverVel_pinv invKinVel(sawyerChain);
  KDL::ChainIkSolverPos_NR_JL invKinPos(sawyerChain, qLow, qUpp, fwdKin, invKinVel);

  // Function template for position IK solver
  IkSolver invKin = [&invKinPos](const KDL::JntArray& qInit, const KDL::Frame& pGoal, KDL::JntArray& qSoln) {
    return invKinPos.CartToJnt(qInit, pGoal, qSoln);
  };

  // Run the test:
  int seed = 82025;
  runGeneralPosIkTest(seed, fwdKin, invKin, qLow, qUpp, "KDL_pinv_NR");
}

/*
 * Run the benchmark test on all five versions of the position solver.
 * Note: each test uses the same seed so that the test problems are identical for each solver.
 */
TEST(sns_ik, pos_ik_SNS_test) {
    runSnsPosIkTest(82025, sns_ik::VelocitySolveType::SNS); }
TEST(sns_ik, pos_ik_SNS_Base_test) {
    runSnsPosIkTest(82025, sns_ik::VelocitySolveType::SNS_Base); }
TEST(sns_ik, pos_ik_SNS_Optimal_test) {
    runSnsPosIkTest(82025, sns_ik::VelocitySolveType::SNS_Optimal); }
TEST(sns_ik, pos_ik_SNS_OptimalScaleMargin_test) {
    runSnsPosIkTest(82025, sns_ik::VelocitySolveType::SNS_OptimalScaleMargin); }
TEST(sns_ik, pos_ik_SNS_Fast_test) {
    runSnsPosIkTest(82025, sns_ik::VelocitySolveType::SNS_Fast); }
TEST(sns_ik, pos_ik_SNS_FastOptimal_test) {
    runSnsPosIkTest(82025, sns_ik::VelocitySolveType::SNS_FastOptimal); }

/*************************************************************************************************/
// Run all the tests that were declared with TEST()
int main(int argc, char **argv){
  ros::Time::init();
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
