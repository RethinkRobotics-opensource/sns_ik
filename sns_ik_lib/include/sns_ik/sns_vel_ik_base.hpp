/** @file sns_vel_ik_base.hpp
 *
 * @brief The file provides the basic implementation of the SNS-IK velocity solver
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

#ifndef SNS_IK_LIB__SNS_VEL_IK_H_
#define SNS_IK_LIB__SNS_VEL_IK_H_

#include <Eigen/Dense>
#include <memory>
#include <vector>

#include "sns_ik_base.hpp"
// #include "sns_linear_solver.hpp"

namespace sns_ik {

class SnsVelIkBase : public SnsIkBase{

public:

  // Smart pointer typedefs. Note: all derived classes MUST override these smart pointers.
  typedef std::shared_ptr<SnsVelIkBase> Ptr;
  typedef std::unique_ptr<SnsVelIkBase> uPtr;

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
   *    dqLow <= dq <= dqUpp      ( bounds set in constructor or setBounds() )
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


  /**
   * Solve a velocity IK problem with null-space bias task as the secondary goal.
   *
   *  This method implements a simplified version of "Algorithm 4: SNS algorithm for multiple tasks"
   *  where the total number of tasks is two and the secondary goal is a configuration-space velocity
   *  task.
   *
   * Solve for joint velocity dq and task scales (s and sCS):
   *
   *  maximize: s
   *  subject to:
   *    s * dx = J * dq
   *    sCS * dqCS = (I - pinv(J)*J) * dq;
   *    0 < s <= 1
   *    0 < sCS <= 1
   *    dqLow <= dq <= dqUpp      ( bounds set in constructor or setBounds() )
   *
   * @param J: Jacobian matrix, mapping from joint to task space. Size = [nTask, nJoint]
   * @param dx: task space velocity vector (primary goal). Length = nTask
   * @param dqCS: configuration space velocity (secondary goal). Length = nJnt
   * @param[out] dq: joint velocity solution. Length = nJoint
   * @param[out] taskScale: task scale for primary goal.  fwdKin(dq) = taskScale*dx
   *                            taskScale == 1.0  --> task was feasible
   *                            taskScale < 1.0  --> task was infeasible and had to be scaled
   * @param[out] taskScaleCS: task scale for secondary goal.
   * @return: ExitCode::Success: the algorithm worked correctly and satisfied the problem statement
   *          otherwise: something went wrong, exit code specifics the type of problem
   */
  ExitCode solve(const Eigen::MatrixXd& J, const Eigen::VectorXd& dx, const Eigen::VectorXd& dqCS,
                      Eigen::VectorXd* dq, double* taskScale, double* taskScaleCS);

protected:

  /*
   * protected constructor: require factory method to create an object.
   */


private:

  SnsVelIkBase(int nJnt) : SnsIkBase(nJnt) {};

};  // class SnsVelIkBase

}  // namespace sns_ik

#endif  // SNS_IK_LIB__SNS_VEL_IK_H_
