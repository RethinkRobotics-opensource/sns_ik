/** @file sawyer_model.cpp
 *
 * @brief The file provides simple functions to return robot models for testing.
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

#include "sawyer_model.hpp"

#include <ros/console.h>

namespace sns_ik {
namespace sawyer_model {

/*************************************************************************************************/

KDL::Chain getSawyerKdlChain(std::vector<std::string>* jointNames)
{
  if (jointNames){
    *jointNames = {"right_j0", "right_j1", "right_j2", "right_j3", "right_j4", "right_j5", "right_j6"};
  }

  using namespace KDL;
  Chain chain;

  chain.addSegment(Segment("right_arm_mount", Joint(Joint::None), KDL::Frame(Rotation::RPY(0.0,0.0,0.0), Vector(0.0, 0.0, 0.08))));
  chain.addSegment(Segment("right_l0", Joint(Joint::RotZ), KDL::Frame(Rotation::RPY(-1.57079632679, 1.57079632679, 0), Vector(0.081, 0.05, 0.237))));
  chain.addSegment(Segment("right_l1", Joint(Joint::RotZ), KDL::Frame(Rotation::RPY(1.57079632679, 0, 0), Vector(0.0, -0.14, 0.1425))));
  chain.addSegment(Segment("right_l2", Joint(Joint::RotZ), KDL::Frame(Rotation::RPY(-1.57079632679, 0, 0), Vector(0.0, -0.042, 0.26))));
  chain.addSegment(Segment("right_l3", Joint(Joint::RotZ), KDL::Frame(Rotation::RPY(1.57079632679, 0, 0), Vector(0.0, -0.125, -0.1265))));
  chain.addSegment(Segment("right_l4", Joint(Joint::RotZ), KDL::Frame(Rotation::RPY(-1.57079632679, 0, 0), Vector(0.0, 0.031, 0.275))));
  chain.addSegment(Segment("right_l5", Joint(Joint::RotZ), KDL::Frame(Rotation::RPY(-1.57079632679, -0.17453, 3.1416), Vector(0.0, -0.11, 0.1053))));
  chain.addSegment(Segment("right_l6", Joint(Joint::RotZ), KDL::Frame(Rotation::RPY(0, 0, 1.570796), Vector(0, 0, 0.0245))));
  chain.addSegment(Segment("right_hand", Joint(Joint::None), KDL::Frame(Rotation::RPY(0.0,0.0,0.0), Vector(0.0, 0.0, 0.0))));

  return chain;
}

/*************************************************************************************************/

void getSawyerJointLimits(KDL::JntArray* qLow, KDL::JntArray* qUpp,
                         KDL::JntArray* vMax, KDL::JntArray* aMax)
{
  if (!qLow || !qUpp || !vMax || !aMax) { ROS_ERROR("Bad input!"); return; }
  int nJnt = 7;  // Sawer has seven joints
  *qLow = KDL::JntArray(nJnt);
  *qUpp = KDL::JntArray(nJnt);
  *vMax = KDL::JntArray(nJnt);
  *aMax = KDL::JntArray(nJnt);
  (*qLow)(0) = -3.0503;
  (*qUpp)(0) = 3.0503;
  (*vMax)(0) = 1.6;
  (*aMax)(0) = 8.0;
  (*qLow)(1) = -3.8095;
  (*qUpp)(1) = 2.2736;
  (*vMax)(1) = 1.3;
  (*aMax)(1) = 8.0;
  (*qLow)(2) = -3.0426;
  (*qUpp)(2) = 3.0426;
  (*vMax)(2) = 1.9;
  (*aMax)(2) = 8.0;
  (*qLow)(3) = -3.0439;
  (*qUpp)(3) = 3.0439;
  (*vMax)(3) = 1.9;
  (*aMax)(3) = 8.0;
  (*qLow)(4) = -2.9761;
  (*qUpp)(4) = 2.9761;
  (*vMax)(4) = 3.4;
  (*aMax)(4) = 10.0;
  (*qLow)(5) = -2.9761;
  (*qUpp)(5) = 2.9761;
  (*vMax)(5) = 3.4;
  (*aMax)(5) = 10.0;
  (*qLow)(6) = -3.14;
  (*qUpp)(6) = 3.14;
  (*vMax)(6) = 4.5;
  (*aMax)(6) = 10.0;
}

/*************************************************************************************************/

}  // namespace sns_ik
}  // namespace sawyer_model
