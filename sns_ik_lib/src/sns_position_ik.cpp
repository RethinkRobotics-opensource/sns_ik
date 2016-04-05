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
#include <sns_ik/sns_velocity_ik.hpp>
#include <iostream>

namespace sns_ik{

SNSPositionIK::SNSPositionIK(KDL::Chain chain, std::shared_ptr<SNSVelocityIK> velocity_ik) :
    m_chain(chain),
    m_ikVelSolver(velocity_ik),
    m_positionFK(chain),
    m_jacobianSolver(chain),
    m_linearMaxStepSize(0.2),
    m_angularMaxStepSize(0.2),
    m_maxIterations(100),
    m_dt(0.002),
    dt_step(0.02)
{
}

SNSPositionIK::~SNSPositionIK()
{
}
/*
int SNSPositionIK::CartToJnt(const KDL::JntArray& joint_seed,
                             const KDL::Frame& goal_pose,
                             const KDL::JntArray& joint_ns_bias,
                             const MatrixD& ns_jacobian,
                             const std::vector<int>& ns_indicies,
                             KDL::JntArray* return_joints,
                             const KDL::Twist& tolerances)
{
  // TODO: use tolerance twist
  double linearTolerance = 1e-5;
  double angularTolerance = 1e-5;
  VectorD jl_low = m_ikVelSolver->getJointLimitLow();
  VectorD jl_high = m_ikVelSolver->getJointLimitHigh();
  VectorD maxJointVel = m_ikVelSolver->getJointVelocityMax();

  // initialize variables
  bool solutionFound = false;
  KDL::JntArray q_i = joint_seed;
  KDL::Frame pose_i;
  int n_dof = joint_seed.rows();
  StackOfTasks sot(1);
  sot[0].desired = VectorD::Zero(6);
  double stepScale = 1.0;

  // If there's a nullspace bias, create a secondary task
  if (joint_ns_bias.rows()) {
    Task nsTask;
    nsTask.jacobian = ns_jacobian;
    nsTask.desired = VectorD::Zero(joint_ns_bias.rows());
    // the desired task to apply the NS bias will change with each iteration
    sot.push_back(nsTask);
  }

  double L, theta;
  double lineErr, rotErr;
  VectorD qDot(n_dof);
  KDL::Jacobian jacobian;
  jacobian.resize(q_i.rows());
  KDL::Vector rotAxis, trans;
  KDL::Rotation rot;
  //double sf;

  int ii;
  for (ii = 0; ii < m_maxIterations; ++ii) {
    if (!calcPositionAndError(q_i, goal_pose, &pose_i, &L, &theta, &trans, &rotAxis)) {
      // ERROR
      return -1;
    }
    rotErr = theta;
    lineErr = L;

    //std::cout << ii << ": Cartesian error: " << L << " m, " << theta << " rad" << std::endl;

    if (lineErr <= linearTolerance && rotErr <= angularTolerance) {
      solutionFound = true;
      break;
    }

    if (L > stepScale * m_linearMaxStepSize) {
      trans = (stepScale * m_linearMaxStepSize / L) * trans;
    }

    if (theta > stepScale * m_angularMaxStepSize) {
      theta = stepScale * m_angularMaxStepSize;
    }

    // Calculate the desired Cartesian twist
    //sot[0].desired.head<3>() = (1.0/m_dt) * trans.data;
    //sot[0].desired.tail<3>() = theta/m_dt * rotAxis.data;
    sot[0].desired(0) = trans.data[0] / m_dt;
    sot[0].desired(1) = trans.data[1] / m_dt;
    sot[0].desired(2) = trans.data[2] / m_dt;
    sot[0].desired(3) = theta * rotAxis.data[0] / m_dt;
    sot[0].desired(4) = theta * rotAxis.data[1] / m_dt;
    sot[0].desired(5) = theta * rotAxis.data[2] / m_dt;

    if (m_jacobianSolver.JntToJac(q_i, jacobian) < 0)
    {
      // ERROR
      std::cout << "JntToJac failed" << std::endl;
      return -1;
    }
    sot[0].jacobian = jacobian.data;

    if (joint_ns_bias.rows()) {
      for (size_t jj = 0; jj < joint_ns_bias.rows(); ++jj) {
        // This calculates a "nullspace velocity".
        // There is an arbitrary scale factor which will be set by the max scale factor.
        int indx = ns_indicies[jj];
        double vel = 0.1 * (joint_ns_bias(jj) - q_i(indx)) / m_dt; // TODO: step size needs to be optimized
        // TODO: may want to limit the NS velocity to 50% of max joint velocity
        //vel = std::max(-0.5*maxJointVel(indx), std::min(0.5*maxJointVel(indx), vel));
        sot[1].desired(jj) = vel;
      }

    }

    //sf = m_ikVelSolver->getJointVelocity(&qDot, sot, q_i.data);
    m_ikVelSolver->getJointVelocity(&qDot, sot, q_i.data);

    q_i.data += m_dt * qDot;

    for (int j = 0; j < jl_low.rows(); ++j) {
      q_i.data[j] = std::max(std::min(q_i.data[j], jl_high[j]), jl_low[j]);
    }

    //std::cout << "    q: " << q_i.data.transpose() << std::endl;
    //std::cout << "    cart vel: " << sot[0].desired.transpose() << std::endl;
    //std::cout << "    qDot: " << qDot.transpose() << std::endl;

    if (qDot.norm() < 1e-5) {  // TODO: config param
      //std::cout << "ERROR: Solution stuck, iter: "<<ii<<", error: " << lineErr << " m, " << rotErr << " rad" << std::endl;
      return -2;
    }

//    if (sf < stepScale) {
//      stepScale = std::max(0.8*stepScale, 0.01);
//      //std::cout << "New step scale (" << ii << "): " << stepScale << std::endl;
//    } else if (stepScale < 1.0) {
//      stepScale = std::min(1.1*stepScale, 1.0);
//    }
  }

  if (solutionFound) {
    *return_joints = q_i;
    //std::cout << "Solution Found in "<< ii <<" iterations!" << std::endl;
    return 1;  // TODO: return success/fail code
  } else {
    //std::cout << "Reached max iterations:, error: " << lineErr << " m, " << rotErr << " rad" << std::endl;
    return -1;
  }
}
*/

bool SNSPositionIK::calcPositionAndError(const KDL::JntArray& q,
                                         const KDL::Frame& goal,
                                         KDL::Frame* pose,
                                         double* errL,
                                         double* errR,
                                         KDL::Vector* trans,
                                         KDL::Vector* rotAxis)
{
  if (m_positionFK.JntToCart(q, *pose) < 0)
  {
    // ERROR
    std::cout << "JntToCart failed" << std::endl;
    return false;
  }

  // Calculate the offset transform
  *trans = goal.p - pose->p;
  *errL = trans->Norm();
  KDL::Rotation rot = goal.M * pose->M.Inverse();
  *errR = rot.GetRotAngle(*rotAxis);  // returns [0 ... pi]
  return true;
}

int SNSPositionIK::CartToJnt(const KDL::JntArray& joint_seed,
                             const KDL::Frame& goal_pose,
                             const KDL::JntArray& joint_ns_bias,
                             const MatrixD& ns_jacobian,
                             const std::vector<int>& ns_indicies,
                             KDL::JntArray* return_joints,
                             const KDL::Twist& tolerances)
{
  // TODO: use tolerance twist
  double linearTolerance = 1e-5;
  double angularTolerance = 1e-5;
  VectorD jl_low = m_ikVelSolver->getJointLimitLow();
  VectorD jl_high = m_ikVelSolver->getJointLimitHigh();
  VectorD maxJointVel = m_ikVelSolver->getJointVelocityMax();

  // initialize variables
  bool solutionFound = false;
  KDL::JntArray q_i = joint_seed;
  KDL::JntArray q_next;
  KDL::Frame pose_i;
  int n_dof = joint_seed.rows();
  StackOfTasks sot(1);
  sot[0].desired = VectorD::Zero(6);
  double stepScale = 1.0;
  double Tstep_max = INF;
  double dt_step_self;

  q_0=q_i;
/*
  // If there's a nullspace bias, create a secondary task
  if (joint_ns_bias.rows()) {
    Task nsTask;
    nsTask.jacobian = ns_jacobian;
    nsTask.desired = VectorD::Zero(joint_ns_bias.rows());
    // the desired task to apply the NS bias will change with each iteration
    sot.push_back(nsTask);
  }
  */


/*
  Task nsTask;
  nsTask.jacobian = MatrixD::Identity(n_dof,n_dof);
  nsTask.desired = VectorD::Zero(n_dof);
  // the desired task to apply the NS bias will change with each iteration
  sot.push_back(nsTask);
*/



  //double L, theta;
  double lineErr, rotErr;
  VectorD qDot(n_dof);
  KDL::Jacobian jacobian;
  jacobian.resize(q_i.rows());
  KDL::Vector rotAxis, trans;
  KDL::Rotation rot;


  if (!calcPositionAndError(q_i, goal_pose, &pose_i, &lineErr, &rotErr, &trans, &rotAxis)) {
    // ERROR
    return -1;
  }

  int ii;
  for (ii = 0; ii < m_maxIterations; ++ii) {
    //std::cout << ii << ": Cartesian error: " << L << " m, " << theta << " rad" << std::endl;

    if (lineErr <= linearTolerance && rotErr <= angularTolerance) {
      solutionFound = true;
      break;
    }

    // Calculate the desired Cartesian twist
    sot[0].desired(0) = trans.data[0] / m_dt;
    sot[0].desired(1) = trans.data[1] / m_dt;
    sot[0].desired(2) = trans.data[2] / m_dt;
    sot[0].desired(3) = rotErr * rotAxis.data[0] / m_dt;
    sot[0].desired(4) = rotErr * rotAxis.data[1] / m_dt;
    sot[0].desired(5) = rotErr * rotAxis.data[2] / m_dt;
//    sot[0].desired(0) = 1000 * trans.data[0];
//    sot[0].desired(1) = 1000 * trans.data[1];
//    sot[0].desired(2) = 1000 * trans.data[2];
//    sot[0].desired(3) = 1000 * rotErr * rotAxis.data[0];
//    sot[0].desired(4) = 1000 * rotErr * rotAxis.data[1];
//    sot[0].desired(5) = 1000 * rotErr * rotAxis.data[2];

    if (m_jacobianSolver.JntToJac(q_i, jacobian) < 0)
    {
      // ERROR
      std::cout << "JntToJac failed" << std::endl;

      return -1;
    }
    sot[0].jacobian = jacobian.data;
/*
      for (size_t jj = 0; jj < n_dof; ++jj) {
        // This calculates a "nullspace velocity".
        // There is an arbitrary scale factor which will be set by the max scale factor.
        double vel =  100*(q_0(jj) - q_i(jj)) ; // TODO: step size needs to be optimized
        // TODO: may want to limit the NS velocity to 50% of max joint velocity
        //vel = std::max(-0.5*maxJointVel(indx), std::min(0.5*maxJointVel(indx), vel));
        sot[1].desired(jj) = vel;
      }
*/

    m_ikVelSolver->setLoopPeriod(m_dt);
    m_ikVelSolver->getJointVelocity(&qDot, sot, q_i.data);


    // Max time step to guarantee joint limits
    dt_step = 0.2;//std::min(2.0*dt_step, 0.2);
    Tstep_max = INF;
    for (int j = 0; j < jl_low.rows(); ++j) {
          //q_i.data[j] = std::max(std::min(q_i.data[j], jl_high[j]), jl_low[j]);
      if (qDot(j) > 1e-9) {
        Tstep_max = std::min(Tstep_max, (jl_high[j] - q_i.data[j]) / qDot(j));
      } else if (qDot(j) < -1e-9) {
        Tstep_max = std::min(Tstep_max, (jl_low[j] - q_i.data[j]) / qDot(j));
      }
    }
    dt_step=(dt_step>Tstep_max)?Tstep_max:(dt_step<0.001)?0.001:dt_step;



    // search for a time step that decreases the error
    double currentErr = lineErr + rotErr;
    double stepErr = INF;


    while (stepErr > currentErr && dt_step >= 0.001) {
      q_next.data = q_i.data + dt_step * qDot;
      calcPositionAndError(q_next, goal_pose, &pose_i, &lineErr, &rotErr, &trans, &rotAxis);
      stepErr = lineErr + rotErr;
      dt_step /= 2.0;
    }
    if (stepErr < currentErr) dt_step *= 2.0;



    //remove stuck situation
    if (qDot.norm() < 1e-9  || dt_step < 0.001) {  // TODO: config param

        for (size_t jj = 0; jj < n_dof; ++jj) {
          qDot(jj) =  1*( - q_i(jj)) ;
        }
        dt_step=0.02;
        q_next.data = q_i.data + dt_step * qDot;

    }



     q_i = q_next;
  //  std::cout << ii <<": dt("<< dt_step<<"), stepErr: " << stepErr << "("<<lineErr<<", "<<rotErr<<")" << std::endl;

    //std::cout << "    q: " << q_i.data.transpose() << std::endl;
    //std::cout << "    cart vel: " << sot[0].desired.transpose() << std::endl;
    //std::cout << "    qDot: " << qDot.transpose() << std::endl;

    if (qDot.norm() < 1e-9) {  // TODO: config param
      std::cout << "ERROR: Solution stuck, iter: "<<ii<<", qDot.norm: "<<qDot.norm()<<", error: " << lineErr << " m, " << rotErr << " rad" << std::endl;
      return -2;
    }

  }

  if (solutionFound) {
    *return_joints = q_i;
  //  std::cout << "Solution Found in "<< ii <<" iterations!" << std::endl;
    return 1;  // TODO: return success/fail code
  } else {
    //std::cout << "Reached max iterations:, error: " << lineErr << " m, " << rotErr << " rad" << std::endl;
    return -1;
  }
}

}  // namespace sns_ik
