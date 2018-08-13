/*! \file sns_ik_vel_test.cpp
 * \brief Unit Test: sns_ik velocity solver
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
#include <kdl/chain.hpp>
#include <kdl/chainiksolvervel_pinv.hpp>
#include <kdl/chainfksolvervel_recursive.hpp>
#include <map>
#include <ros/console.h>
#include <ros/duration.h>
#include <ros/time.h>

#include "rng_utilities.hpp"
#include "sawyer_model.hpp"
#include <sns_ik/sns_ik.hpp>

/*
 * Define a common interface to call both SNS and KDL solvers for velocity IK
 */
typedef std::function<int(const KDL::JntArray&, const KDL::Twist&,
                          KDL::JntArray&, double& taskScale)> IkSolver;

/*
 * This file provides a unit and regression test for the top-level SNS-IK solver, along with the
 * various solvers that are derived from it.
 */

/*************************************************************************************************
 *                               Velocity IK Test Params                                         *
 *************************************************************************************************/

// How many velocity IK tests should be run?
static const int VEL_IK_TEST_COUNT = 250;

// Should the velocity IK tests print the results (beyond pass/fail) to terminal?
static const bool VEL_IK_TEST_VERBOSE = true;

// Tolerance for checks on forward kinematics
static const double VEL_IK_TEST_FK_LIN_TOL = 1e-8;  // meters per second
static const double VEL_IK_TEST_FK_ANG_TOL = 1e-8;  // radians per second

/*************************************************************************************************
 *                               Utilities Functions                                             *
 *************************************************************************************************/

struct VelTestResult {
  ros::Duration solveTime;  // duration spent in the solver across all calls
  int exitCode;  // exit code returned by the solver
  double fkLinErr;  // error in the linear speed
  double fkAngErr;  // error in the angular speed
  bool fkValid;  // did the forward-kinematics on the solution meet tolerance?
  int nPass;  // number of tests that passed
  int nFail;  // number of tests that failed
  double taskScale;  // scaling applied to the task
};

/*************************************************************************************************/

/*
 * Computes the maximum absolute error in the forward kinematics at the velocity level.
 * @param q: joint angles
 * @param dq: joint rates
 * @param fwdKin: solver to use for the forward-kinematics
 * @param twist: expected twist at the end-effector
 * @param[out] linErr: linear speed error
 * @param[out] angErr: angular speed error
 */
void checkForwardKinematicVelocity( const KDL::JntArray& q, const KDL::JntArray& dq,
                                    KDL::ChainFkSolverVel_recursive& fwdKin, const KDL::Twist& twistSns,
                                    double* linErr, double* angErr)
{
  if (!linErr) { ROS_ERROR("linErr is nullptr!"); return; }
  if (!angErr) { ROS_ERROR("angErr is nullptr!"); return; }

  // compute the goal endpoint twist
  KDL::JntArrayVel jointData(q, dq);
  KDL::FrameVel frameVel;
  fwdKin.JntToCart(jointData, frameVel);
  KDL::Twist twistFk = frameVel.GetTwist();

  // compute the difference
  KDL::Twist twistErr = twistFk - twistSns;

  // compute the error in each component:
  KDL::Vector rot = twistErr.rot;
  KDL::Vector vel = twistErr.vel;
  *linErr = rot.Norm();
  *angErr = vel.Norm();
}

/*************************************************************************************************/

/*
 * A function to run a benchmarking test on the velocity IK solvers.
 *
 * Structure:
 * - generate N test poses in joint space
 * - use forward kinematics to compute end-effector twist
 * - solve the IK problem using ikSolver
 * - measure and report solve time and success rate
 *
 * @param seed: seed to pass to the RNG on each call
 *              if seed == 0, then RNG seed is not updated
 * @param chain:  kinematic chain for the model system
 * @param qLow:  lower joint limits
 * @param qUpp:  upper joint limits
 * @param vMax:  maximum joint speed
 * @param fwdKin:  forward kinematics solver
 * @param ikSolver:  inverse kinematics solver
 * @return: velocity solver test result
 */
VelTestResult runVelIkSingleTest(int seed, const KDL::JntArray& qLow, const KDL::JntArray& qUpp,
  const KDL::JntArray& vMax, KDL::ChainFkSolverVel_recursive& fwdKin,  IkSolver& ikSolver)
{
  // compute a random joint velocity and velocity for the solution and the corresponding test pose
  KDL::JntArray qTest = sns_ik::rng_util::getRngBoundedJoints(seed, qLow, qUpp);
  seed = 0;  // let the RNG update the sequence automatically after the first call
  KDL::JntArray vMin = vMax;
  for (int i = 0; i < int(vMin.rows()); i++) { vMin(i) = -vMin(i); }
  KDL::JntArray wTest = sns_ik::rng_util::getRngBoundedJoints(seed, vMin, vMax);
  KDL::JntArrayVel qwTest(qTest, wTest);

  // compute the goal endpoint twist
  KDL::FrameVel zTest;
  fwdKin.JntToCart(qwTest, zTest);
  KDL::Twist twist = zTest.GetTwist();

  // initialize the output for the IK solver
  KDL::JntArray dqSolve;
  double taskScale;

  // initialize the test result:
  VelTestResult result;

  // run benchmarking:
  ros::Time startTime = ros::Time::now();
  result.exitCode = ikSolver(qTest, twist, dqSolve, taskScale);
  result.solveTime += ros::Time::now() - startTime;
  result.taskScale = taskScale;  // primary task scale
  KDL::Twist scaledTwist = twist * result.taskScale;
  checkForwardKinematicVelocity(qTest, dqSolve, fwdKin, scaledTwist, &(result.fkLinErr), &(result.fkAngErr));
  result.fkValid = result.fkLinErr <= VEL_IK_TEST_FK_LIN_TOL && result.fkAngErr <= VEL_IK_TEST_FK_ANG_TOL;
  if (result.fkValid && result.exitCode >= 0) {
    result.nPass = 1; result.nFail = 0;
  } else {
    result.nPass = 0; result.nFail = 1;
  }
  return result;
}


/*************************************************************************************************/

/*
 * @param seed: seed to pass to the RNG on each call
 *              if seed == 0, then seed is ignored
 * @param fwdKin: forward kinematics solver for velocity
 * @param ikSolver: inverse kinematics solver for velocity (to be tested)
 * @param qLow: lower bounds on joints (for generating the test problem)
 * @param qLow: upper bounds on joints (for generating the test problem)
 * @param vMax: maximum joint speed (for generating the test problem)
 * @param solverName: solver name, used for logging only
 * @param nTest: number of tests to run
 */
void runGeneralVelIkTest(int seed, KDL::ChainFkSolverVel_recursive& fwdKin, IkSolver& ikSolver,
                  const KDL::JntArray& qLow, const KDL::JntArray& qUpp,
                  const KDL::JntArray& vMax, std::string solverName, int nTest = VEL_IK_TEST_COUNT)
{
  // Set up the data structure for the results:
  VelTestResult results;  // accumulate results here
  VelTestResult tmp; // individual test output here

  // run the test many times and accumulate the results
  bool success = true;
  for (int iTest = 0; iTest < nTest; iTest++) {
    seed++;
    tmp = runVelIkSingleTest(seed, qLow, qUpp, vMax, fwdKin, ikSolver);
    if (iTest == 0) {
      results = tmp;
    } else {
      results.nPass += tmp.nPass;
      results.nFail += tmp.nFail;
      results.solveTime += tmp.solveTime;
      results.fkLinErr += tmp.fkLinErr;
      results.fkAngErr += tmp.fkAngErr;
      ASSERT_LE(tmp.fkLinErr, VEL_IK_TEST_FK_LIN_TOL);
      ASSERT_LE(tmp.fkAngErr, VEL_IK_TEST_FK_ANG_TOL);
    }
  }

  // Print out the results:
  if (VEL_IK_TEST_VERBOSE) {
    int nPass = results.nPass;
    int nFail = results.nFail;
    int nTotal = nPass + nFail;
    double total = static_cast<double>(nTotal);
    ROS_INFO("Velocity IK Test  -->  pass: %d, fail: %d  --  %f ms  --  linErr: %e, angErr: %e  --  %s",
              nPass, nFail, 1000.0 * results.solveTime.toSec() / total,
              results.fkLinErr / total, results.fkAngErr / total,
              solverName.c_str());
  };

  // Final check:
  ASSERT_TRUE(success);
}

/*************************************************************************************************/

void runSnsVelkTest(int seed, sns_ik::VelocitySolveType solverType) {

  // Create a sawyer model:
  std::vector<std::string> jointNames;
  KDL::Chain sawyerChain = sns_ik::sawyer_model::getSawyerKdlChain(&jointNames);
  KDL::JntArray qLow, qUpp, vMax, aMax;
  sns_ik::sawyer_model::getSawyerJointLimits(&qLow, &qUpp, &vMax, &aMax);

  // Create a forward-kinematics solver:
  KDL::ChainFkSolverVel_recursive fwdKin(sawyerChain);

  // Create a SNS-IK solver:
  sns_ik::SNS_IK ikSolver(sawyerChain, qLow, qUpp, vMax, aMax, jointNames);
  ikSolver.setVelocitySolveType(solverType);

  // Function template for velocity IK solver
  IkSolver invKin = [&ikSolver](const KDL::JntArray& qInit, const KDL::Twist& dpGoal,
                                KDL::JntArray& dqSoln, double& taskScale) {
    int exitCode = ikSolver.CartToJntVel(qInit, dpGoal, dqSoln);
    std::vector<double> taskScaleVec;
    ikSolver.getTaskScaleFactors(taskScaleVec);
    if (taskScaleVec.empty()) {
      ROS_ERROR("No task scale!");
      taskScale = -1.0;  // This is what SNS-IK uses as an error flag...
    } else {
      taskScale = taskScaleVec.front();
    }
    return exitCode;
  };

  // Run the test:
  /*
   * Note: This seems to cause a segfault in ubuntu 16.04 with Eigen 3.3.4, but it
   * works fine in Ubuntu 14.04 with Eigen 3.2.0. FIXME
   */
  runGeneralVelIkTest(seed, fwdKin, invKin, qLow, qUpp, vMax, sns_ik::toStr(solverType));
}

/*************************************************************************************************
 *                                        Tests                                                  *
 *************************************************************************************************/

/*
 * Run the benchmark test on all five versions of the velocity solver.
 * Note: each test uses the same seed so that the test problems are identical for each solver.
 */
TEST(sns_ik, vel_ik_SNS_test) {
    runSnsVelkTest(23539, sns_ik::VelocitySolveType::SNS); }
TEST(sns_ik, vel_ik_SNS_Base_test) {
    runSnsVelkTest(23539, sns_ik::VelocitySolveType::SNS_Base); }
// FIXME - uncomment to run test. For details, see
// https://github.com/RethinkRobotics-opensource/sns_ik/issues/87
// TEST(sns_ik, vel_ik_SNS_Optimal_test) {
//    runSnsVelkTest(23539, sns_ik::VelocitySolveType::SNS_Optimal); }
TEST(sns_ik, vel_ik_SNS_OptimalScaleMargin_test) {
    runSnsVelkTest(23539, sns_ik::VelocitySolveType::SNS_OptimalScaleMargin); }
TEST(sns_ik, vel_ik_SNS_Fast_test) {
    runSnsVelkTest(23539, sns_ik::VelocitySolveType::SNS_Fast); }
TEST(sns_ik, vel_ik_SNS_FastOptimal_test) {
    runSnsVelkTest(23539, sns_ik::VelocitySolveType::SNS_FastOptimal); }

/*************************************************************************************************/

// Run all the tests that were declared with TEST()
int main(int argc, char **argv){
  ros::Time::init();
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
