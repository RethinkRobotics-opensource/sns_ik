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
// Author: Ian McMahon

#ifndef SNS_IK_HPP
#define SNS_IK_HPP

#include <string>
#include <memory>
#include <kdl/chain.hpp>
#include <kdl/frames.hpp>
#include <kdl/jntarray.hpp>
#include <kdl/chainfksolverpos_recursive.hpp>
#include <sns_ik/sns_position_ik.hpp>

namespace sns_ik {

  enum VelocitySolveType { SNS,
                           SNS_Optimal,
                           SNS_OptimalScaleMargin,
                           SNS_Fast,
                           SNS_FastOptimal,
                           SNS_Base
                         };


  /*
   * Convert velocity solver type to a string (for logging)
   * @param solverType: solve type to convert
   * @return: string representation of the solver type
   */
  std::string toStr(const sns_ik::VelocitySolveType& type);

  // Forward declare SNS Velocity Base Class
  class SNSVelocityIK;
  class SNS_IK
  {
  public:
    SNS_IK(const std::string& base_link, const std::string& tip_link,
           const std::string& URDF_param="/robot_description",
           double loopPeriod=0.01, double eps=1e-5,
           sns_ik::VelocitySolveType type=sns_ik::SNS);

    SNS_IK(const KDL::Chain& chain,
           const KDL::JntArray& q_min, const KDL::JntArray& q_max,
           const KDL::JntArray& v_max, const KDL::JntArray& a_max,
           const std::vector<std::string>& jointNames,
           double loopPeriod=0.01, double eps=1e-5,
           sns_ik::VelocitySolveType type=sns_ik::SNS);

    ~SNS_IK();

    bool setVelocitySolveType(VelocitySolveType type);

    inline bool getPositionSolver(std::shared_ptr<sns_ik::SNSPositionIK>& positionSolver) {
      positionSolver=m_ik_pos_solver;
      return m_initialized;
    }

    inline bool getVelocitySolver(std::shared_ptr<sns_ik::SNSVelocityIK>& velocitySolver) {
      velocitySolver=m_ik_vel_solver;
      return m_initialized;
    }

    inline bool getKDLChain(KDL::Chain& chain) {
      chain=m_chain;
      return m_initialized;
    }

    inline bool getKDLLimits(KDL::JntArray& lb, KDL::JntArray& ub, KDL::JntArray& vel, KDL::JntArray& accel) {
      lb=m_lower_bounds;
      ub=m_upper_bounds;
      vel=m_velocity;
      accel=m_acceleration;
      return m_initialized;
    }

    inline bool getJointNames(std::vector<std::string>& jointNames) {
      jointNames = m_jointNames;
      return m_initialized;
    }

    bool setMaxJointVelocity(const KDL::JntArray& vel);

    bool setMaxJointAcceleration(const KDL::JntArray& accel);

    int CartToJnt(const KDL::JntArray &q_init,
                  const KDL::Frame &p_in,
                  KDL::JntArray &q_out,
                  const KDL::Twist& bounds=KDL::Twist::Zero())
    { return CartToJnt(q_init, p_in, KDL::JntArray(0), std::vector<std::string>(),
                       q_out, bounds);
    }

    // Assumes the NS bias is for all the joints in the correct order
    int CartToJnt(const KDL::JntArray &q_init, const KDL::Frame &p_in,
                  const KDL::JntArray& q_bias,
                  KDL::JntArray &q_out,
                  const KDL::Twist& bounds=KDL::Twist::Zero())
    { return CartToJnt(q_init, p_in, q_bias, m_jointNames,
                       q_out, bounds);
    }

    int CartToJnt(const KDL::JntArray &q_init, const KDL::Frame &p_in,
                  const KDL::JntArray& q_bias,
                  const std::vector<std::string>& biasNames,
                  KDL::JntArray &q_out,
                  const KDL::Twist& bounds=KDL::Twist::Zero());

    int CartToJntVel(const KDL::JntArray& q_in,
                     const KDL::Twist& v_in,
                     KDL::JntArray& qdot_out)
    { return CartToJntVel(q_in, v_in, KDL::JntArray(0), std::vector<std::string>(),
                          KDL::JntArray(0), qdot_out); }

    // Assumes the NS bias is for all the joints in the correct order
    int CartToJntVel(const KDL::JntArray& q_in,
                     const KDL::Twist& v_in,
                     const KDL::JntArray& q_bias,
                     KDL::JntArray& qdot_out)
    { return CartToJntVel(q_in, v_in, q_bias, m_jointNames, KDL::JntArray(0), qdot_out); }

    int CartToJntVel(const KDL::JntArray& q_in,
                         const KDL::Twist& v_in,
                         const KDL::JntArray& q_bias,
                         const KDL::JntArray& q_vel_bias,
                         KDL::JntArray& qdot_out)
    { return CartToJntVel(q_in, v_in, q_bias, m_jointNames, q_vel_bias, qdot_out); }

    int CartToJntVel(const KDL::JntArray& q_in,
                     const KDL::Twist& v_in,
                     const KDL::JntArray& q_bias,
                     const std::vector<std::string>& biasNames,
                     const KDL::JntArray& q_vel_bias,
                     KDL::JntArray& qdot_out);

    // Nullspace gain should be specified between 0 and 1.0
    double getNullspaceGain() { return m_nullspaceGain; }
    void setNullspaceGain(double gain)
    { m_nullspaceGain = std::max(std::min(gain, 1.0), 0.0); }

    // Set time step in seconds
    void setLoopPeriod(double loopPeriod);
    double getLoopPeriod() { return m_loopPeriod; }

    bool getTaskScaleFactors(std::vector<double>& scaleFactors);

  private:
    bool m_initialized;
    double m_eps;
    double m_loopPeriod;
    double m_nullspaceGain;
    VelocitySolveType m_solvetype;
    KDL::Chain m_chain;
    KDL::JntArray m_lower_bounds, m_upper_bounds, m_velocity, m_acceleration;
    enum JointType { Revolute, Prismatic, Continuous };
    std::vector<JointType> m_types;
    std::vector<std::string> m_jointNames;

    std::vector<KDL::JntArray> m_solutions;
    std::shared_ptr<SNSVelocityIK> m_ik_vel_solver;
    std::shared_ptr<SNSPositionIK> m_ik_pos_solver;
    std::shared_ptr<KDL::ChainJntToJacSolver> m_jacobianSolver;

    void initialize();

    bool nullspaceBiasTask(const KDL::JntArray& q_bias,
                           const std::vector<std::string>& biasNames,
                           Eigen::MatrixXd* jacobian, std::vector<int>* indicies);

  };
}  //namespace
#endif
