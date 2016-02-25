/* test for the IK library*/

#include <sns_ik/sns_ik_math_utils.hpp>
#include <sns_ik/sns_velocity_ik.hpp>
#include <sns_ik/osns_velocity_ik.hpp>
#include <sns_ik/osns_sm_velocity_ik.hpp>
#include <sns_ik/sns_position_ik.hpp>

#include <Eigen/Dense>
#include <iostream>
#include <kdl/chain.hpp>

using namespace Eigen;
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

  VectorD l = VectorD::Ones(7);

  SNSVelocityIK ikVelSolver(7, 0.01);
  ikVelSolver.setJointsCapabilities(-3.0*l, 3.0*l, l, 0.5*l);
  ikVelSolver.getJointVelocity(&jointVelocity, sot, joints);

  std::cout << "SNS Velocity IK result: " << std::endl
            << jointVelocity.transpose() << std::endl;
  std::cout << "-----------------------------" << std::endl;

  OSNSVelocityIK ikVelSolver_osns(7, 0.01);
  ikVelSolver_osns.setJointsCapabilities(-3.0*l, 3.0*l, l, 0.5*l);
  ikVelSolver_osns.getJointVelocity(&jointVelocity, sot, joints);

  std::cout << "Optimal SNS Velocity IK result: " << std::endl
      << jointVelocity.transpose() << std::endl;
  std::cout << "-----------------------------" << std::endl;

  OSNSsmVelocityIK ikVelSolver_osns_sm(7, 0.01);
  ikVelSolver_osns_sm.setJointsCapabilities(-3.0*l, 3.0*l, l, 0.5*l);
  ikVelSolver_osns_sm.getJointVelocity(&jointVelocity, sot, joints);

  std::cout << "Optimal SNS w/ sm Velocity IK result: " << std::endl
      << jointVelocity.transpose() << std::endl;
  std::cout << "-----------------------------" << std::endl;

  KDL::Chain chain;
  KDL::JntArray jointSeed(7);
  for (int ii = 0; ii < 7; ++ii) {
    chain.addSegment(KDL::Segment());
    jointSeed(ii) = joints(ii);
  }

  SNSPositionIK positionIK(chain, ikVelSolver);

  KDL::Frame goal;  // TODO: randomize
  KDL::JntArray goalJoints;
  KDL::Twist tolerances;  // not currently used
  positionIK.CartToJnt(jointSeed, goal, &goalJoints, tolerances);

  std::cout << "Position IK result: " << std::endl
             << goalJoints.data.transpose() << std::endl;
}
