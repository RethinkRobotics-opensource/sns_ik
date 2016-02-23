/* test for the IK library*/

#include <sns_ik/sns_ik_math_utils.hpp>
#include <sns_ik/sns_velocity_ik.hpp>
#include <sns_ik/sns_position_ik.hpp>

#include <Eigen/Dense>
#include <iostream>
#include <kdl/chain.hpp>

using namespace Eigen;
using namespace sns_ik;

int main(int argc, char** argv) {
  SNSVelocityIK ikVelSolver(7, 0.01);
  StackOfTasks sot;
  Task task;

  VectorD jointVelocity;

  task.jacobian = MatrixD::Random(3,7);
  task.desired =  MatrixD::Random(3,1);

  std::cout << "desired: " << task.desired.transpose() << std::endl;
  std::cout << "jacobian: " << task.jacobian << std::endl;

  sot.push_back(task);

  VectorD l = VectorD::Ones(7);
  ikVelSolver.setJointsCapabilities(-3.0*l, 3.0*l, l, 0.5*l);

  VectorD joints = VectorD::Random(7);
  ikVelSolver.getJointVelocity(&jointVelocity, sot, joints);

  std::cout << "Velocity IK result: "<< jointVelocity.transpose() << std::endl;

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

  //std::count << "Positin IK result: " << goalJoints.transpose() << std::endl;
}
