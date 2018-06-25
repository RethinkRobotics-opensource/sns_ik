/** @file sns_acc_ik_base.hpp
 *
 * @brief The file provides the basic implementation of the SNS-IK acceleration solver
 *
 * @author Matthew Kelly
 *
 * This file provides a set of functions that return simple kinematic chains that are used for the
 * unit tests. This allows unit tests to run quickly without depending on external URDF files.
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

#ifndef SNS_IK_LIB__SNS_ACC_IK_H_
#define SNS_IK_LIB__SNS_ACC_IK_H_

#include <Eigen/Dense>
#include <memory>

#include "sns_linear_solver.hpp"

namespace sns_ik {


class SnsAccIkBase {

public:

  // Smart pointer typedefs. Note: all derived classes MUST override these smart pointers.
  typedef std::shared_ptr<SnsAccIkBase> Ptr;
  typedef std::unique_ptr<SnsAccIkBase> uPtr;

  enum class ExitCode {
    Success,  // successfully solver the optimization problem
    BadUserInput,  // user input failed basic checks (eg. inconsistent matrix size)
    InfeasibleTask,  // there is no feasible solution to the primary task
    InternalError  // there was an internal error in the solver (should never happen...)
  };

  /**
   * Create a default solver with nJnt joints and no bounds on joint acceleration
   * @param nJnt: number of joints in the robot model (columns in the jacobian)
   * @return: acceleration solver iff successful, nullptr otherwise
   */
  static std::unique_ptr<SnsAccIkBase> create(int nJnt);

  /**
   * Create a default solver with nJnt joints and infinite bounds on the joint acceleration.
   * @param ddqLow: lower bound on the acceleration of each joint
   * @param ddqUpp: upper bound on the acceleration of each joint
   * @return: acceleration solver iff successful, nullptr otherwise
   */
  static std::unique_ptr<SnsAccIkBase> create(const Eigen::ArrayXd& ddqLow,
                                              const Eigen::ArrayXd& ddqUpp);

  // Make sure that class is cleaned-up correctly
  virtual ~SnsAccIkBase() {};

  /**
   * Set the bounds on joint acceleration.
   * Set a bound to infinity (or an arbitrarily large value) to disable it.
   * This method will update the number of joints in the solver
   * Requirements: inputs must be the same size and ddqLow < ddqUpp
   * @param ddqLow: lower bound on the acceleration of each joint
   * @param ddqUpp: upper bound on the acceleration of each joint
   * @return: true iff successful
   */
  bool setAccBnd(const Eigen::ArrayXd& ddqLow, const Eigen::ArrayXd& ddqUpp);

  /**
   * Solve a acceleration IK problem with no null-space bias of joint-space optimization.
   *
   *  This method implements   "Algorithm 1: SNS algorithm"   from the paper:
   *  "Prioritized multi-task motion control of redundant robots under hard joint constraints"
   *   by: Fabrizio Flacco, Alessandro De Luca, Oussama Khatib
   *
   * Solve for joint acceleration ddqUpp and task scale s:
   *
   *  maximize: s
   *  subject to:
   *    s * ddx = J * ddq + dJ*dq
   *    0 < s <= 1
   *    ddqLow <= ddqUpp <= ddqUpp      ( bounds set in constructor or setVelBnd() )
   *
   * @param J: Jacobian matrix, mapping from joint to task space. Size = [nTask, nJoint]
   * @param dJdq: the product of Jacobian derivative and joint velocity. Length = nTask
   * @param dx: task acceleration vector. Length = nTask
   * @param[out] ddqUpp: joint acceleration solution. Length = nJoint
   * @param[out] taskScale: task scale.  fwdKin(q, dq, ddq) = taskScale*ddx
   *                            taskScale == 1.0  --> task was feasible
   *                            taskScale < 1.0  --> task was infeasible and had to be scaled
   * @return: Success: the algorithm worked correctly and satisfied the problem statement
   *          otherwise: something went wrong
   */
  ExitCode solve(const Eigen::MatrixXd& J, const Eigen::VectorXd& dJdq, const Eigen::VectorXd& dx,
                      Eigen::VectorXd* ddqUpp, double* taskScale);


protected:

  /*
   * protected constructor: require factory method to create an object.
   */
  SnsAccIkBase(int nJnt) : nJnt_(nJnt), ddqLow_(nJnt), ddqUpp_(nJnt) {};

  /*
   * Check that ddqLow_ <= ddqUpp <= ddqUpp_
   * @param ddqUpp: joint acceleration to test
   * @return: true iff ddqLow <= ddq <= ddqUpp
   *          if ddqUpp.empty() return false
   */
  bool checkAccBnd(const Eigen::VectorXd& ddq);

  /*
   * Solve the following equation for the variable ddqUpp:
   *    J * W * (ddq - dqNull) = ddx - dJdq - J*ddqNull
   * This equation appears throughout the SNS-IK papers. One example is in block "D" of Algorithm 1:
   *  "Control of Redundant Robots Under Hard Joint Constraint: Saturation in the Null Space"
   *   by: Fabrizio Flacco, Alessandro De Luca, Oussama Khatib
   *
   * PRECONDITION:
   * --> Assumes that linSolver_ has been initialized with J*W
   *
   * @param J: Jacobian matrix, mapping from joint to task space. Size = [nTask, nJoint]
   * @param dJdq: the product of Jacobian derivative and joint velocity. Length = nTask
   * @param JW: J*W = Jacobian projected onto the active joints. Size = [nTask, nJoint]
   * @param ddqNull: null-space joint acceleration. Size = nJoint
   * @param ddx: task acceleration vector. Length = nTask
   * @param[out] ddq: joint acceleration solution. Length = nJoint
   * @param[out, opt] resErr: residual error (norm-squared)
   * @return: true --> success!
   *          false --> invalid input or other error
   */
  bool solveProjectionEquation(const Eigen::MatrixXd& J, const Eigen::VectorXd& dJdq,
                               const Eigen::MatrixXd& JW, const Eigen::VectorXd& dqNull,
                               const Eigen::VectorXd& ddx, Eigen::VectorXd* ddq, double* resErr);

  /*
   * This method implements Algorithm 2 (and a bit of Algorithm 1) from the paper:
   *  "Control of Redundant Robots Under Hard Joint Constraint: Saturation in the Null Space"
   *   by: Fabrizio Flacco, Alessandro De Luca, Oussama Khatib
   *
   * PRECONDITION:
   * --> Assumes that linSolver_ has been initialized with J*W
   *
   * @param J: Jacobian matrix, mapping from joint to task space. Size = [nTask, nJoint]
   * @param dJdq: the product of Jacobian derivative and joint velocity. Length = nTask
   * @param JW: J*W = Jacobian projected onto the active joints. Size = [nTask, nJoint]
   * @param ddx: task acceleration vector. Length = nTask
   * @param ddqUpp: joint acceleration. Length = nJoint
   * @param ddqUpp: joint acceleration. Length = nJoint
   * @param jntIsFree: which joints are free to saturate? Length = nJoint
   * @param[out] taskScale: task scale factor
   * @param[out] jntIdx: index corresponding to the most critical joint that is free
   * @param[out] resErr: residual error (norm-squared) in the linear solve
   * @return: true --> success!
   *          false --> invalid input or other error
   */
  ExitCode computeTaskScalingFactor(const Eigen::MatrixXd& J, const Eigen::MatrixXd& JW, const Eigen::VectorXd& ddx,
                                              const Eigen::VectorXd& ddq, const std::vector<bool>& jntIsFree,
                                              double* taskScale, int* jntIdx, double* resErr);

   /*
    * This algorithm computes the scale factor that is associated with a given joint, but considering
    * both the sensativity of the joint (a) and the distance to the upper and lower limits.
    * @param low: lower margin
    * @param upp: upper margin
    * @param a: margin scale factor
    * @return: joint scale factor
    */
   static double findScaleFactor(double low, double upp, double a);

private:

  int nJnt_; //!< number of joints

  Eigen::ArrayXd ddqLow_;  //!< lower bound on joint acceleration
  Eigen::ArrayXd ddqUpp_;  //!< upper bound on joint acceleration

  SnsLinearSolver linSolver_;  //!< linear solver for the core SNS-IK algorithm

};  // class SnsAccIkBase

}  // namespace sns_ik

#endif  // SNS_IK_LIB__SNS_ACC_IK_H_
