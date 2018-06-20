/** @file sns_vel_ik_base.hpp
 *
 * @brief The file provides the basic implementation of the SNS-IK velocity solver
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

#ifndef SNS_IK_LIB__SNS_VEL_IK_H_
#define SNS_IK_LIB__SNS_VEL_IK_H_

#include <Eigen/Dense>
#include <memory>
#include <vector>

#include "sns_linear_solver.hpp"

namespace sns_ik {

class SnsVelIkBase {

public:

  // Smart pointer typedefs. Note: all derived classes MUST override these smart pointers.
  typedef std::shared_ptr<SnsVelIkBase> Ptr;
  typedef std::unique_ptr<SnsVelIkBase> uPtr;

  enum class ExitCode {
    Success,  // successfully solved the optimization problem
    BadUserInput,  // user input failed basic checks (eg. inconsistent matrix size)
    InfeasibleTask,  // there is no feasible solution to the primary task
    InternalError  // there was an internal error in the solver (should never happen...)
  };

  /**
   * Create a default solver with nJnt joints and no bounds on joint velocity
   * @param nJnt: number of joints in the robot model (columns in the jacobian)
   * @return: velocity solver iff successful, nullptr otherwise
   */
  static std::unique_ptr<SnsVelIkBase> create(int nJnt);

  /**
   * Create a default solver with constant bounds on the joint velocity
   * @param dqLow: lower bound on the velocity of each joint
   * @param dqUpp: upper bound on the velocity of each joint
   * @return: velocity solver iff successful, nullptr otherwise
   */
  static std::unique_ptr<SnsVelIkBase> create(const Eigen::ArrayXd& dqLow,
                                              const Eigen::ArrayXd& dqUpp);

  // Make sure that class is cleaned-up correctly
  virtual ~SnsVelIkBase() {};

  /**
   * Set the bounds on joint velocity.
   * Set a bound to infinity (or an arbitrarily large value) to disable it.
   * This method will update the number of joints in the solver
   * Requirements: inputs must be the same size and dqLow < dqUpp
   * @param dqLow: lower bound on the velocity of each joint
   * @param dqUpp: upper bound on the velocity of each joint
   * @return: true iff successful
   */
  bool setVelBnd(const Eigen::ArrayXd& dqLow, const Eigen::ArrayXd& dqUpp);

  /**
   * Solve a velocity IK problem with no null-space bias of joint-space optimization.
   *
   *  This method implements   "Algorithm 1: SNS algorithm"   from the paper:
   *  "Control of Redundant Robots Under Hard Joint Constraint: Saturation in the Null Space"
   *   by: Fabrizio Flacco, Alessandro De Luca, Oussama Khatib
   *
   * Solve for joint velocity dq and task scale s:
   *
   *  maximize: s
   *  subject to:
   *    s * dx = J * dq
   *    0 < s <= 1
   *    dqLow <= dq <= dqUpp      ( bounds set in constructor or setVelBnd() )
   *
   * @param J: Jacobian matrix, mapping from joint to task space. Size = [nTask, nJoint]
   * @param dx: task velocity vector. Length = nTask
   * @param[out] dq: joint velocity solution. Length = nJoint
   * @param[out] taskScale: task scale.  fwdKin(dq) = taskScale*dx
   *                            taskScale == 1.0  --> task was feasible
   *                            taskScale < 1.0  --> task was infeasible and had to be scaled
   * @return: ExitCode::Success: the algorithm worked correctly and satisfied the problem statement
   *          otherwise: something went wrong, exit code specifics the type of problem
   */
  ExitCode solve(const Eigen::MatrixXd& J, const Eigen::VectorXd& dx,
                      Eigen::VectorXd* dq, double* taskScale);


protected:

  /*
   * protected constructor: require factory method to create an object.
   */
  SnsVelIkBase(int nJnt) : nJnt_(nJnt), dqLow_(nJnt), dqUpp_(nJnt) {};

  /*
   * Check that dqLow_ <= dq <= dqUpp_
   * @param dq: joint velocity to test
   * @return: true iff qLow <= q <= qUpp
   *          if dq.empty() return false
   */
  bool checkVelBnd(const Eigen::VectorXd& dq);

  /*
   * Solve the following equation for the variable dq:
   *    J * W * (dq - dqNull) = dx - J*dqNull
   * This equation appears throughout the SNS-IK papers. One example is in block "D" of Algorithm 1:
   *  "Control of Redundant Robots Under Hard Joint Constraint: Saturation in the Null Space"
   *   by: Fabrizio Flacco, Alessandro De Luca, Oussama Khatib
   *
   * PRECONDITION:
   * --> Assumes that linSolver_ has been initialized with J*W
   *
   * @param J: Jacobian matrix, mapping from joint to task space. Size = [nTask, nJoint]
   * @param JW: J*W = Jacobian projected onto the active joints. Size = [nTask, nJoint]
   * @param dqNull: null-space joint velocity. Size = nJoint
   * @param dx: task velocity vector. Length = nTask
   * @param[out] dq: joint velocity solution. Length = nJoint
   * @param[out, opt] resErr: residual error (norm-squared)
   * @return: true --> success!
   *          false --> invalid input or other error
   */
  ExitCode solveProjectionEquation(const Eigen::MatrixXd& J, const Eigen::MatrixXd& JW,
                                   const Eigen::VectorXd& dqNull, const Eigen::VectorXd& dx,
                                   Eigen::VectorXd* dq, double* resErr);

  /*
   * This method implements Algorithm 2 (and a bit of Algorithm 1) from the paper:
   *  "Control of Redundant Robots Under Hard Joint Constraint: Saturation in the Null Space"
   *   by: Fabrizio Flacco, Alessandro De Luca, Oussama Khatib
   *
   * PRECONDITION:
   * --> Assumes that linSolver_ has been initialized with J*W
   *
   * @param J: Jacobian matrix, mapping from joint to task space. Size = [nTask, nJoint]
   * @param JW: J*W = Jacobian projected onto the active joints. Size = [nTask, nJoint]
   * @param dx: task velocity vector. Length = nTask
   * @param dq: joint velocity. Length = nJoint
   * @param jntIsFree: which joints are free to saturate? Length = nJoint
   * @param[out] taskScale: task scale factor
   * @param[out] jntIdx: index corresponding to the most critical joint that is free
   * @param[out] resErr: residual error (norm-squared) in the linear solve
   * @return: ExitCode::Success: the algorithm worked correctly and satisfied the problem statement
   *          otherwise: something went wrong, exit code specifics the type of problem
   */
  ExitCode computeTaskScalingFactor(const Eigen::MatrixXd& J, const Eigen::MatrixXd& JW,
                                const Eigen::VectorXd& dx, const Eigen::VectorXd& dq,
                                const std::vector<bool>& jntIsFree,
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

  Eigen::ArrayXd dqLow_;  //!< lower bound on joint velocity
  Eigen::ArrayXd dqUpp_;  //!< upper bound on joint velocity

  SnsLinearSolver linSolver_;  //!< linear solver for the core SNS-IK algorithm

};  // class SnsVelIkBase

}  // namespace sns_ik

#endif  // SNS_IK_LIB__SNS_VEL_IK_H_
