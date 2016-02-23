/*! \file sns_position_ik.cpp
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

#include <sns_ik/sns_position_ik.hpp>

#include <iostream>

using namespace sns_ik;

SNSPositionIK::SNSPositionIK(KDL::Chain chain, SNSVelocityIK velocity_ik) :
    m_chain(chain),
    m_ikVelSolver(velocity_ik),
    m_positionFK(chain),
    m_jacobianSolver(chain)
{
}

SNSPositionIK::~SNSPositionIK()
{
}

int SNSPositionIK::CartToJnt(const KDL::JntArray& joint_seed,
                             const KDL::Frame& desired_end_effector_pose,
                             KDL::JntArray* return_joints,
                             const KDL::Twist& tolerances)
{
  // TODO: config params
  // TODO: use tolerance twist
  double linearTolerance = 1e-6;
  double angularTolerance = 1e-6;
  double linearMaxStepSize = 0.05;
  double angularMaxStepSize = 0.025;
  int maxInterations = 100;
  double dt = 0.2;

  bool solutionFound = false;
  KDL::JntArray q_i = joint_seed;
  KDL::Frame pose_i, pose_delta;
  int n_dof = joint_seed.rows();
  StackOfTasks sot(1);
  sot[0].desired = VectorD::Zero(6);

  for (int ii = 0; ii < maxInterations; ++ii) {
    if (m_positionFK.JntToCart(q_i, pose_i) < 0)
    {
      // ERROR
      return -1;
    }

    // TODO: need to check math here
    pose_delta = pose_i.Inverse() * desired_end_effector_pose;
    Eigen::Vector3d trans(pose_delta.p.data);
    double L = trans.norm();
    KDL::Vector rotAxis;
    double theta = pose_delta.M.GetRotAngle(rotAxis);
    Eigen::Vector3d rotAxisVec(rotAxis.data);
    double rotSign = 1;
    if (theta < 0.0) {
      rotSign = -1;
      theta = -theta;
    }

    if (L <= linearTolerance && theta <= angularTolerance) {
      solutionFound = true;
      break;
    }

    if (L > linearMaxStepSize) {
      trans = linearMaxStepSize / L * trans;
    }
    if (theta > M_PI) {
      theta = 2 * M_PI - theta;
      rotSign = -rotSign;
    }
    if (theta > angularMaxStepSize) {
      theta = angularMaxStepSize;
    }

    // Calculate the desired Cartesian twist
    sot[0].desired.head<3>() = (1.0/dt) * trans;
    sot[0].desired.head<3>() = rotSign * theta/dt * rotAxisVec;

    KDL::Jacobian jacobian;
    m_jacobianSolver.JntToJac(q_i, jacobian);
    sot[0].jacobian = jacobian.data;

    std::cout << "desired: " << sot[0].desired.transpose() << std::endl;
    std::cout << "jacobian: " << std::endl << sot[0].jacobian << std::endl;
    VectorD q_ii(q_i.data);
    std::cout << "q_ii: " << q_ii.transpose() << std::endl;

    VectorD qDot(n_dof);
    //m_ikVelSolver.getJointVelocity(&qDot, sot, q_ii);

    q_i.data += dt * qDot;
  }

  if (solutionFound) {
      *return_joints = q_i;
      return 1;  // TODO: return success/fail code
    } else {
      return -1;
    }
}
