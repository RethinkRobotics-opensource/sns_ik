/** @file sns_acc_ik_base.hpp
 *
 * @brief The file provides the basic implementation of the SNS-IK acceleration solver
 *
 * @author Matthew Kelly
 * @author Andy Park
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
#include "sns_ik_base.hpp"

namespace sns_ik {


class SnsAccIkBase : SnsIkBase{

public:

  // Smart pointer typedefs. Note: all derived classes MUST override these smart pointers.
  typedef std::shared_ptr<SnsAccIkBase> Ptr;
  typedef std::unique_ptr<SnsAccIkBase> uPtr;

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
   *    ddqLow <= ddqUpp <= ddqUpp      ( bounds set in constructor or setBounds() )
   *
   * @param J: Jacobian matrix, mapping from joint to task space. Size = [nTask, nJoint]
   * @param dJdq: the product of Jacobian derivative and joint velocity. Length = nTask
   * @param ddx: task acceleration vector. Length = nTask
   * @param[out] ddq: joint acceleration solution. Length = nJoint
   * @param[out] taskScale: task scale.  fwdKin(q, dq, ddq) = taskScale*ddx
   *                            taskScale == 1.0  --> task was feasible
   *                            taskScale < 1.0  --> task was infeasible and had to be scaled
   * @return: Success: the algorithm worked correctly and satisfied the problem statement
   *          otherwise: something went wrong
   */
  ExitCode solve(const Eigen::MatrixXd& J, const Eigen::VectorXd& dJdq, const Eigen::VectorXd& dx,
                      Eigen::VectorXd* ddqUpp, double* taskScale);

  /**
   * Solve a acceleration IK problem with null-space bias task as the secondary goal.
   *
   *  This method implements a simplified version of "Algorithm 4: SNS algorithm for multiple tasks"
   *  where the total number of tasks is two and the secondary goal is a configuration-space velocity
   *  task.
   *
   * Solve for joint acceleration ddqUpp and task scales (s and sCS):
   *
   *  maximize: s
   *  subject to:
   *    s * ddx = J * ddq + dJ*dq
   *    sCS * ddqCS = (I - pinv(J)*J) * ddq;
   *    0 < s <= 1
   *    0 < sCS <= 1
   *    ddqLow <= ddq <= ddqUpp      ( bounds set in constructor or setBounds() )
   *
   * @param J: Jacobian matrix, mapping from joint to task space. Size = [nTask, nJoint]
   * @param dJdq: the product of Jacobian derivative and joint velocity. Length = nTask
   * @param ddx: task acceleration vector. Length = nTask
   * @param ddqCS: configuration space acceleration (secondary goal). Length = nJnt
   * @param[out] ddq: joint acceleration solution. Length = nJoint
   * @param[out] taskScale: task scale.  fwdKin(q, dq, ddq) = taskScale*ddx
   *                            taskScale == 1.0  --> task was feasible
   *                            taskScale < 1.0  --> task was infeasible and had to be scaled
   * @param[out] taskScaleCS: task scale for secondary goal.
   * @return: Success: the algorithm worked correctly and satisfied the problem statement
   *          otherwise: something went wrong
   */
  ExitCode solve(const Eigen::MatrixXd& J, const Eigen::VectorXd& dJdq, const Eigen::VectorXd& ddx,
                 const Eigen::VectorXd& ddqCS, Eigen::VectorXd* ddq, double* taskScale, double* taskScaleCS);

protected:

  /*
   * protected constructor: require factory method to create an object.
   */

private:

  SnsAccIkBase(int nJnt) : SnsIkBase(nJnt) {};

};  // class SnsAccIkBase

}  // namespace sns_ik

#endif  // SNS_IK_LIB__SNS_ACC_IK_H_
