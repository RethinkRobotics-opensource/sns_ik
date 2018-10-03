/*! \file sns_ik_acc_test.cpp
 * \brief Unit Test: sns_ik acceleration solver
 * \author Andy Park
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
 * Define a common interface to call both SNS and KDL solvers for acceleration IK
 */
typedef std::function<int(const KDL::JntArray&, const KDL::JntArray&, const KDL::Twist&,
                          KDL::JntArray&, double& taskScale)> IkSolver;

/*
 * This file provides a unit and regression test for the top-level SNS-IK solver, along with the
 * various solvers that are derived from it.
 */

/*************************************************************************************************
 *                               Acceleration IK Test Params                                         *
 *************************************************************************************************/

// How many acceleration IK tests should be run?
static const int ACC_IK_TEST_COUNT = 250;

// Should the acceleration IK tests print the results (beyond pass/fail) to terminal?
static const bool ACC_IK_TEST_VERBOSE = true;

// Tolerance for checks on forward kinematics
static const double ACC_IK_TEST_FK_LIN_TOL = 1e-8;  // meters per second^2
static const double ACC_IK_TEST_FK_ANG_TOL = 1e-8;  // radians per second^2

// Margins for joint angles and rates with respect to their limits
static const double TEST_ANGLE_MARGIN = 1e-1;
static const double TEST_RATE_MARGIN = 1e-1;

/*************************************************************************************************
 *                               Utilities Functions                                             *
 *************************************************************************************************/

struct AccTestResult {
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

void forwardKinematicsAcceleration(const KDL::JntArray& q, const KDL::JntArray& dq,
                                   const KDL::JntArray& ddq, sns_ik::SNS_IK pSolver, KDL::Twist *accFK)
{
  // compute Jacobian
  Eigen::MatrixXd jacobian;

  if (!pSolver.getJacobian(q, &jacobian))
  {
    ROS_ERROR("Jacobian computation failed!");
  }

  // compute Jacobian dot
  Eigen::MatrixXd jacobianDot;
  if (!pSolver.getJacobianDot(q, dq, &jacobianDot))
  {
    ROS_ERROR("Jacobian dot computation failed!");
  }

  // compute the goal endpoint acc
  Eigen::VectorXd accFKtmp = jacobian*ddq.data + jacobianDot*dq.data;

  accFK->vel(0) = accFKtmp(0);
  accFK->vel(1) = accFKtmp(1);
  accFK->vel(2) = accFKtmp(2);
  accFK->rot(0) = accFKtmp(3);
  accFK->rot(1) = accFKtmp(4);
  accFK->rot(2) = accFKtmp(5);
}


/*
 * Computes the maximum absolute error in the forward kinematics at the acceleration level.
 * @param q: joint angles
 * @param dq: joint rates
 * @param ddq: joint accelerations
 * @param accSns: expected acc at the end-effector
 * @param[out] linErr: linear acc error
 * @param[out] angErr: angular acc error
 */
void checkForwardKinematicAcceleration( const KDL::JntArray& q, const KDL::JntArray& dq, const KDL::JntArray& ddq,
                                    const KDL::Twist& accSns, sns_ik::SNS_IK pSolver, double* linErr, double* angErr)
{
  if (!linErr) { ROS_ERROR("linErr is nullptr!"); return; }
  if (!angErr) { ROS_ERROR("angErr is nullptr!"); return; }

  // get forward kinematics in acc
  KDL::Twist accFk;
  forwardKinematicsAcceleration(q, dq, ddq, pSolver, &accFk);

  // compute the difference
  KDL::Twist twistErr = accSns - accFk;

  // compute the error in each component:
  KDL::Vector rot = twistErr.rot;
  KDL::Vector vel = twistErr.vel;
  *linErr = rot.Norm();
  *angErr = vel.Norm();
}

/*************************************************************************************************/

/*
 * A function to run a benchmarking test on the acceleration IK solvers.
 *
 * Structure:
 * - generate N test poses in joint space
 * - use forward kinematics to compute end-effector acc
 * - solve the IK problem using ikSolver
 * - measure and report solve time and success rate
 *
 * @param seed: seed to pass to the RNG on each call
 *              if seed == 0, then RNG seed is not updated
 * @param chain:  kinematic chain for the model system
 * @param qLow:  lower joint limits
 * @param qUpp:  upper joint limits
 * @param vMax:  maximum joint speed
 * @param aMax:  maximum joint acc
 * @param fwdKin:  forward kinematics solver
 * @param ikSolver:  inverse kinematics solver
 * @return: acceleration solver test result
 */
AccTestResult runAccIkSingleTest(int seed, const KDL::JntArray& qLow, const KDL::JntArray& qUpp,
  const KDL::JntArray& vMax, const KDL::JntArray& aMax, IkSolver& ikSolver, sns_ik::SNS_IK pSolver)
{
  // compute a random joint acceleration and acceleration for the solution and the corresponding test pose
  KDL::JntArray qLowTmp = qLow;
  KDL::JntArray qUppTmp = qUpp;
  for (int i = 0; i < int(qLow.rows()); i++) { qLowTmp(i) = qLowTmp(i) + TEST_ANGLE_MARGIN; }
  for (int i = 0; i < int(qUpp.rows()); i++) { qUppTmp(i) = qUppTmp(i) - TEST_ANGLE_MARGIN; }
  KDL::JntArray qTest = sns_ik::rng_util::getRngBoundedJoints(seed, qLowTmp, qUppTmp);
  seed = 0;  // let the RNG update the sequence automatically after the first call

  KDL::JntArray vMin = vMax;
  for (int i = 0; i < int(vMin.rows()); i++) { vMin(i) = -vMin(i); }

  KDL::JntArray qdLowTmp = vMin;
  KDL::JntArray qdUppTmp = vMax;
  for (int i = 0; i < int(qdLowTmp.rows()); i++) { qdLowTmp(i) = qdLowTmp(i) + TEST_RATE_MARGIN; }
  for (int i = 0; i < int(qdUppTmp.rows()); i++) { qdUppTmp(i) = qdUppTmp(i) - TEST_RATE_MARGIN; }
  KDL::JntArray dqTest = sns_ik::rng_util::getRngBoundedJoints(seed, qdLowTmp, qdUppTmp);

  KDL::JntArray aMin = aMax;
  for (int i = 0; i < int(aMin.rows()); i++) { aMin(i) = -aMin(i); }
  KDL::JntArray ddqTest = sns_ik::rng_util::getRngBoundedJoints(seed, aMin, aMax);

  // get forward kinematics in acc
  KDL::Twist accFkTmp;
  forwardKinematicsAcceleration(qTest, dqTest, ddqTest, pSolver, &accFkTmp);

  // initialize the output for the IK solver
  KDL::JntArray ddqSolve(qTest.rows());
  double taskScale;

  // initialize the test result:
  AccTestResult result;

  // run benchmarking:
  ros::Time startTime = ros::Time::now();
  result.exitCode = ikSolver(qTest, dqTest, accFkTmp, ddqSolve, taskScale);

  result.solveTime += ros::Time::now() - startTime;
  result.taskScale = taskScale;  // primary task scale
  KDL::Twist scaledAcc = accFkTmp * result.taskScale;
  checkForwardKinematicAcceleration(qTest, dqTest, ddqSolve, scaledAcc, pSolver, &(result.fkLinErr), &(result.fkAngErr));
  result.fkValid = result.fkLinErr <= ACC_IK_TEST_FK_LIN_TOL && result.fkAngErr <= ACC_IK_TEST_FK_ANG_TOL;
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
 * @param fwdKin: forward kinematics solver for acceleration
 * @param ikSolver: inverse kinematics solver for acceleration (to be tested)
 * @param qLow: lower bounds on joints (for generating the test problem)
 * @param qLow: upper bounds on joints (for generating the test problem)
 * @param vMax: maximum joint speed (for generating the test problem)
 * @param aMax: maximum joint acc (for generating the test problem)
 * @param solverName: solver name, used for logging only
 * @param nTest: number of tests to run
 */
void runGeneralAccIkTest(int seed, sns_ik::SNS_IK pSolver, IkSolver& ikSolver,
                  const KDL::JntArray& qLow, const KDL::JntArray& qUpp,
                  const KDL::JntArray& vMax, const KDL::JntArray& aMax,
                  std::string solverName, int nTest = ACC_IK_TEST_COUNT)
{
  // Set up the data structure for the results:
  AccTestResult results;  // accumulate results here
  AccTestResult tmp; // individual test output here

  // run the test many times and accumulate the results
  bool success = true;
  for (int iTest = 0; iTest < nTest; iTest++) {
    seed++;
    tmp = runAccIkSingleTest(seed, qLow, qUpp, vMax, aMax, ikSolver, pSolver);
    if (iTest == 0) {
      results = tmp;
    } else {
      results.nPass += tmp.nPass;
      results.nFail += tmp.nFail;
      results.solveTime += tmp.solveTime;
      results.fkLinErr += tmp.fkLinErr;
      results.fkAngErr += tmp.fkAngErr;
      ASSERT_LE(tmp.fkLinErr, ACC_IK_TEST_FK_LIN_TOL);
      ASSERT_LE(tmp.fkAngErr, ACC_IK_TEST_FK_ANG_TOL);
    }
  }

  // Print out the results:
  if (ACC_IK_TEST_VERBOSE) {
    int nPass = results.nPass;
    int nFail = results.nFail;
    int nTotal = nPass + nFail;
    double total = static_cast<double>(nTotal);
    ROS_INFO("Acceleration IK Test  -->  pass: %d, fail: %d  --  %f ms  --  linErr: %e, angErr: %e  --  %s",
              nPass, nFail, 1000.0 * results.solveTime.toSec() / total,
              results.fkLinErr / total, results.fkAngErr / total,
              solverName.c_str());
  };

  // Final check:
  ASSERT_TRUE(success);
}

/*************************************************************************************************/

void runSnsAccIkTest(int seed, sns_ik::VelocitySolveType solverType) {

  // Create a sawyer model:
  std::vector<std::string> jointNames;
  KDL::Chain sawyerChain = sns_ik::sawyer_model::getSawyerKdlChain(&jointNames);
  KDL::JntArray qLow, qUpp, vMax, aMax;
  sns_ik::sawyer_model::getSawyerJointLimits(&qLow, &qUpp, &vMax, &aMax);

  // Create a SNS-IK solver:
  sns_ik::SNS_IK ikSolver(sawyerChain, qLow, qUpp, vMax, aMax, jointNames);
  ikSolver.setVelocitySolveType(solverType);

  // Function template for acceleration IK solver
  IkSolver invKin = [&ikSolver](const KDL::JntArray& qInit, const KDL::JntArray& dqInit,
                                const KDL::Twist& ddpGoal, KDL::JntArray& ddqSoln, double& taskScale) {
    int exitCode = ikSolver.CartToJntAcc(qInit, dqInit, ddpGoal, ddqSoln);
    std::vector<double> taskScaleAcc;
    ikSolver.getTaskScaleFactors(taskScaleAcc);
    if (taskScaleAcc.empty()) {
      ROS_ERROR("No task scale!");
      taskScale = -1.0;  // This is what SNS-IK uses as an error flag...
    } else {
      taskScale = taskScaleAcc.front();
    }
    return exitCode;
  };

  // Run the test:
  /*
   * Note: This seems to cause a segfault in ubuntu 16.04 with Eigen 3.3.4, but it
   * works fine in Ubuntu 14.04 with Eigen 3.2.0. FIXME
   */
  runGeneralAccIkTest(seed, ikSolver, invKin, qLow, qUpp, vMax, aMax, sns_ik::toStr(solverType));
}

/*************************************************************************************************
 *                                        Tests                                                  *
 *************************************************************************************************/

/*
 * Run the benchmark test on the acceleration solvers.
 * Note: each test uses the same seed so that the test problems are identical for each solver.
 */
TEST(sns_ik, acc_ik_SNS_test) {
    runSnsAccIkTest(23539, sns_ik::VelocitySolveType::SNS_Base); }

/*************************************************************************************************/

// Run all the tests that were declared with TEST()
int main(int argc, char **argv){
  ros::Time::init();
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
