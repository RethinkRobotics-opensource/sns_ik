/* test for the IK library*/

#include <sns_ik/sns_ik_math_utils.hpp>
#include <sns_ik/sns_velocity_ik.hpp>
#include <sns_ik/osns_velocity_ik.hpp>
#include <sns_ik/osns_sm_velocity_ik.hpp>
#include <sns_ik/fsns_velocity_ik.hpp>
#include <sns_ik/sns_position_ik.hpp>

#include <Eigen/Dense>
#include <iostream>
#include <memory>
#include <kdl/chain.hpp>

using namespace Eigen;
using namespace KDL;
using namespace sns_ik;

// run command:
// rosrun sns_ik_examples test_sns_ik

int main(int argc, char** argv) {
  StackOfTasks sot;
  Task task;

  VectorD jointVelocity;

  task.jacobian = MatrixD::Random(3,7);
  task.desired =  MatrixD::Random(3,1);
  sot.push_back(task);
  VectorD joints = VectorD::Random(7);

  std::cout << "desired: " << task.desired.transpose() << std::endl;
  std::cout << "jacobian: " << std::endl << task.jacobian << std::endl;
  std::cout << "joints: " << joints.transpose() << std::endl;
  std::cout << "-----------------------------" << std::endl;

  VectorD l = VectorD::Ones(7);

  std::shared_ptr<SNSVelocityIK> ikVelSolver(new SNSVelocityIK(7, 0.01));
  ikVelSolver->setJointsCapabilities(-3.0*l, 3.0*l, l, 0.5*l);
  ikVelSolver->getJointVelocity(&jointVelocity, sot, joints);

  std::cout << "SNS Velocity IK result: " << std::endl
            << jointVelocity.transpose() << std::endl;
  std::cout << "-----------------------------" << std::endl;

  OSNSVelocityIK ikVelSolver_osns(7, 0.01);
  ikVelSolver_osns.setJointsCapabilities(-3.0*l, 3.0*l, l, 0.5*l);
  ikVelSolver_osns.getJointVelocity(&jointVelocity, sot, joints);

  std::cout << "Optimal SNS Velocity IK result: " << std::endl
      << jointVelocity.transpose() << std::endl;
  std::cout << "-----------------------------" << std::endl;

  OSNS_sm_VelocityIK ikVelSolver_osns_sm(7, 0.01);
  ikVelSolver_osns_sm.setJointsCapabilities(-3.0*l, 3.0*l, l, 0.5*l);
  ikVelSolver_osns_sm.getJointVelocity(&jointVelocity, sot, joints);

  std::cout << "Optimal SNS w/ sm Velocity IK result: " << std::endl
      << jointVelocity.transpose() << std::endl;
  std::cout << "-----------------------------" << std::endl;

  FSNSVelocityIK ikVelSolver_fsns(7, 0.01);
  ikVelSolver_fsns.setJointsCapabilities(-3.0*l, 3.0*l, l, 0.5*l);
  ikVelSolver_fsns.getJointVelocity(&jointVelocity, sot, joints);

  std::cout << "Fast SNS Velocity IK result: " << std::endl
      << jointVelocity.transpose() << std::endl;
  std::cout << "-----------------------------" << std::endl;

  //Definition of a kinematic chain & add segments to the chain
  KDL::Chain chain;
  chain.addSegment(Segment(Joint(Joint::RotZ),Frame(Vector(0.0,0.0,1.020))));
  chain.addSegment(Segment(Joint(Joint::RotX),Frame(Vector(0.0,0.0,0.480))));
  chain.addSegment(Segment(Joint(Joint::RotX),Frame(Vector(0.0,0.0,0.645))));
  chain.addSegment(Segment(Joint(Joint::RotZ)));
  chain.addSegment(Segment(Joint(Joint::RotX),Frame(Vector(0.0,0.0,0.120))));
  chain.addSegment(Segment(Joint(Joint::RotZ)));
  chain.addSegment(Segment(Joint(Joint::RotX),Frame(Vector(0.0,0.0,0.10))));

  KDL::JntArray q_goal(7), q_tSeed(7);
  for (int ii = 0; ii < 7; ++ii) {
    q_goal(ii) = joints(ii);
  }

  KDL::Frame goal;
  KDL::ChainFkSolverPos_recursive positionFK(chain);
  positionFK.JntToCart(q_goal, goal);

  SNSPositionIK positionIK(chain, ikVelSolver);

  KDL::JntArray goalJoints;
  KDL::Twist tolerances;  // not currently used
  int result = positionIK.CartToJnt(q_tSeed, goal, &goalJoints, tolerances);

  // DEBUG
  //goal.p.data[0] += 0.01;
  //goal.M = KDL::Rotation::RotX(0.1) * goal.M;
  //int result = positionIK.CartToJnt(q_goal, goal, &goalJoints, tolerances);

  std::cout << "Position IK result: " << result << std::endl
            << goalJoints.data.transpose() << std::endl;
}
