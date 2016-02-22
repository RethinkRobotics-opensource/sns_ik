/*! \file sns_position_ik.hpp
 * \brief Basic SNS Position IK solver
 * \author Forrest Rogers-Marcovitz
 */
/*
 *    Copyright 2016 Rethink Robotics
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

#ifndef SNS_IK_POSITION_IK
#define SNS_IK_POSITION_IK

#include <Eigen/Dense>
#include <kdl/chain.hpp>
#include <kdl/chainjnttojacsolver.hpp>
#include <kdl/chainfksolverpos_recursive.hpp>

#include "sns_ik/sns_ik_math_utils.hpp"
#include "sns_ik/sns_velocity_ik.hpp"

namespace sns_ik {

class SNSPositionIK {
  public:
    SNSPositionIK(KDL::Chain chain, SNSVelocityIK velocity_ik);

    ~SNSPositionIK();

    int CartToJnt(const KDL::JntArray& joint_seed,
                  const KDL::Frame& desired_end_effector_pose,
                  KDL::JntArray* return_joints,
                  const KDL::Twist& tolerances);

    KDL::Chain* getChain() { return &m_chain; }

    SNSVelocityIK* getVelocityIK() { return &m_ikVelSolver; }

  private:
    KDL::Chain m_chain;
    SNSVelocityIK m_ikVelSolver;
    KDL::ChainFkSolverPos_recursive m_positionFK;
    KDL::ChainJntToJacSolver m_jacobianSolver;
};

}  // namespace sns_ikl

#endif
