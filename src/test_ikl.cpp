/* test for the IK library*/

#include <sns_ikl/sns_ikl.hpp>

#include <Eigen/Dense>
#include <iostream>

using namespace Eigen;
using namespace IKL_;

int main(int argc, char** argv) {
  
  IKL inverseKinematic;
  StackOfTasks sot;
  Task task;
  
  VectorD jointVelocity;
  
  task.jacobian = MatrixD::Random(3,7);
  task.desired =  MatrixD::Random(3,1);

  sot.push_back(task);
  
  inverseKinematic.getJointVelocity(&jointVelocity,&sot);
  
  cout << "result :"<< endl << jointVelocity << endl;
  

}
