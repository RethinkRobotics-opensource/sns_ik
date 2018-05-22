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
#include <ros/console.h>

#include "sns_ik_math_utils.hpp"

namespace sns_ik {

SNSPositionIK::SNSPositionIK(KDL::Chain chain, std::shared_ptr<SNSVelocityIK> velocity_ik, double eps) :
    m_chain(chain),
    m_ikVelSolver(velocity_ik),
    m_positionFK(chain),
    m_jacobianSolver(chain),
    m_linearMaxStepSize(0.2),
    m_angularMaxStepSize(0.2),
    m_maxIterations(150),
    m_eps(eps),
    m_dt(0.2),
    m_useBarrierFunction(true),
    m_barrierInitAlpha(0.1),
    m_barrierDecay(0.8)
{
}

SNSPositionIK::~SNSPositionIK()
{
}

bool SNSPositionIK::calcPoseError(const KDL::JntArray& q,
                                  const KDL::Frame& goal,
                                  KDL::Frame* pose,
                                  double* errL,
                                  double* errR,
                                  KDL::Vector* trans,
                                  KDL::Vector* rotAxis)
{
  if (m_positionFK.JntToCart(q, *pose) < 0)
  {
    ROS_ERROR("JntToCart failed");
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
                             const Eigen::MatrixXd& ns_jacobian,
                             const std::vector<int>& ns_indicies,
                             const double ns_gain,
                             KDL::JntArray* return_joints,
                             const KDL::Twist& bounds)
{
  Eigen::VectorXd jl_low = m_ikVelSolver->getJointLimitLow();
  Eigen::VectorXd jl_high = m_ikVelSolver->getJointLimitHigh();
  Eigen::VectorXd maxJointVel = m_ikVelSolver->getJointVelocityMax();

  // initialize variables
  bool solutionFound = false;
  KDL::JntArray q_i = joint_seed;
  KDL::Frame pose_i;
  int n_dof = joint_seed.rows();
  std::vector<Task> sot(1);
  sot[0].desired = Eigen::VectorXd::Zero(6);

  // If there's a nullspace bias, create a secondary task
  if (joint_ns_bias.rows()) {
    Task nsTask;
    nsTask.jacobian = ns_jacobian;
    nsTask.desired = Eigen::VectorXd::Zero(joint_ns_bias.rows());
    // the desired task to apply the NS bias will change with each iteration
    sot.push_back(nsTask);
  }

  double theta;
  double lineErr, rotErr;
  Eigen::VectorXd qDot(n_dof);
  KDL::Jacobian jacobian;
  jacobian.resize(q_i.rows());
  KDL::Vector rotAxis, trans;
  KDL::Rotation rot;
  KDL::Twist delta_twist;

  double barrierAlpha = m_barrierInitAlpha;

  int ii;
  for (ii = 0; ii < m_maxIterations; ++ii) {

    if (!calcPoseError(q_i, goal_pose, &pose_i, &lineErr, &rotErr, &trans, &rotAxis)) {
      ROS_ERROR("Failed to calculate pose error!");
      return -1;
    }

    // Check stopping tolerances
    delta_twist = diffRelative(goal_pose, pose_i);

    if (std::abs(delta_twist.vel.x()) <= std::abs(bounds.vel.x()))
      delta_twist.vel.x(0);
    if (std::abs(delta_twist.vel.y()) <= std::abs(bounds.vel.y()))
      delta_twist.vel.y(0);
    if (std::abs(delta_twist.vel.z()) <= std::abs(bounds.vel.z()))
      delta_twist.vel.z(0);
    if (std::abs(delta_twist.rot.x()) <= std::abs(bounds.rot.x()))
      delta_twist.rot.x(0);
    if (std::abs(delta_twist.rot.y()) <= std::abs(bounds.rot.y()))
      delta_twist.rot.y(0);
    if (std::abs(delta_twist.rot.z()) <= std::abs(bounds.rot.z()))
      delta_twist.rot.z(0);

    if(KDL::Equal(delta_twist, KDL::Twist::Zero(), m_eps)) {
      solutionFound = true;
      break;
    }

    // Enforce max linear and rotational step sizes
    if (lineErr > m_linearMaxStepSize) {
      trans = (m_linearMaxStepSize / lineErr) * trans;
    }

    theta = rotErr;
    if (theta > m_angularMaxStepSize) {
      theta = m_angularMaxStepSize;
    }

    // Calculate the desired Cartesian twist
    sot[0].desired(0) = trans.data[0] / m_dt;
    sot[0].desired(1) = trans.data[1] / m_dt;
    sot[0].desired(2) = trans.data[2] / m_dt;
    sot[0].desired(3) = theta * rotAxis.data[0] / m_dt;
    sot[0].desired(4) = theta * rotAxis.data[1] / m_dt;
    sot[0].desired(5) = theta * rotAxis.data[2] / m_dt;

    if (m_jacobianSolver.JntToJac(q_i, jacobian) < 0)
    {
      ROS_ERROR("JntToJac failed");
      return -1;
    }
    sot[0].jacobian = jacobian.data;

    if (joint_ns_bias.rows()) {
      for (size_t jj = 0; jj < joint_ns_bias.rows(); ++jj) {
        // This calculates a "nullspace velocity".
        // There is an arbitrary scale factor which will be set by the max scale factor.
        int indx = ns_indicies[jj];
        double vel = ns_gain * (joint_ns_bias(jj) - q_i(indx)); // TODO: step size needs to be optimized
        // TODO: may want to limit the NS velocity to 50% of max joint velocity
        //vel = std::max(-0.5*maxJointVel(indx), std::min(0.5*maxJointVel(indx), vel));
        sot[1].desired(jj) = vel;
      }

    }

    m_ikVelSolver->getJointVelocity(&qDot, sot, q_i.data);

    if (qDot.norm() < 1e-6) {  // TODO: config param
      ROS_ERROR("qDot.norm() too small!");
      return -2;
    }

    // Update the joint positions
    q_i.data += m_dt * qDot;

    // Apply a decaying barrier function
    // u = upper limit;  l = lower limit
    // B(x) = -log(u - x) - log(-l + x)
    // -alpha * dB(x)/dx === alpha * (1/(x-l) + 1/(x-u))
    if (m_useBarrierFunction && (lineErr > 0.5 || rotErr > 0.5) ) {
      for (int j = 0; j < jl_low.rows(); ++j) {
        // First force the joint within limits.
        // It can not be exactly at the limit since it will cause a division by zero NaN in the barrier function.
        q_i.data[j] = std::max(std::min(q_i.data[j], jl_high[j] - 1e-7), jl_low[j] + 1e-7);
        q_i.data[j] += barrierAlpha * (1/(q_i.data[j] - jl_low[j]) + 1/(q_i.data[j] - jl_high[j]));
      }
      barrierAlpha *= m_barrierDecay;
    } else {
      for (int j = 0; j < jl_low.rows(); ++j) {
        q_i.data[j] = std::max(std::min(q_i.data[j], jl_high[j]), jl_low[j]);
      }
    }
  }

  if (solutionFound) {
    *return_joints = q_i;
    ROS_DEBUG("Solution Found in %d iterations!", ii);
    return 1;  // TODO: return success/fail code
  } else {
    return -1;
  }
}

}  // namespace sns_ik
