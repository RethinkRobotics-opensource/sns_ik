/* test for the IK library*/

#include <sns_ikl/sns_ikl_math_utils.hpp>
#include <sns_ikl/sns_velocity_ik.hpp>

#include <Eigen/Dense>
#include <iostream>

using namespace Eigen;
using namespace sns_ikl;

int main(int argc, char** argv) {
  SNSVelocityIK ikVelSolver;
  StackOfTasks sot;
  Task task;

  VectorD jointVelocity;
  
  task.jacobian = MatrixD::Random(3,7);
  task.desired =  MatrixD::Random(3,1);

  sot.push_back(task);

  VectorD joints = VectorD::Random(7);

  ikVelSolver.getJointVelocity(&jointVelocity, sot, joints);

  std::cout << "result :"<< jointVelocity.transpose() << std::endl;
}
