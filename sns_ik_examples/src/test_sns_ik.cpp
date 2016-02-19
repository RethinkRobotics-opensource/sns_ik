/* test for the IK library*/

#include <sns_ik/sns_ik_math_utils.hpp>
#include <sns_ik/sns_velocity_ik.hpp>

#include <Eigen/Dense>
#include <iostream>

using namespace Eigen;
using namespace sns_ik;

int main(int argc, char** argv) {
  SNSVelocityIK ikVelSolver(7, 0.01);
  StackOfTasks sot;
  Task task;

  VectorD jointVelocity;
  
  task.jacobian = MatrixD::Random(3,7);
  task.desired =  MatrixD::Random(3,1);

  sot.push_back(task);

  VectorD l = VectorD::Ones(7);
  ikVelSolver.setJointsCapabilities(-3.0*l, 3.0*l, l, 0.5*l);

  VectorD joints = VectorD::Random(7);
  ikVelSolver.getJointVelocity(&jointVelocity, sot, joints);

  std::cout << "result :"<< jointVelocity.transpose() << std::endl;
}
