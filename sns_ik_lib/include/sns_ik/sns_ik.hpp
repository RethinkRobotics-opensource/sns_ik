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
#include <sns_ik/sns_velocity_ik.hpp>
#include <sns_ik/osns_velocity_ik.hpp>
#include <sns_ik/osns_sm_velocity_ik.hpp>
#include <sns_ik/sns_position_ik.hpp>

namespace sns_ik {

  enum SolveType { Position, Velocity };

  class SNS_IK
  {
  public:
    SNS_IK(const std::string& base_link, const std::string& tip_link,
           const std::string& URDF_param="/robot_description",
           double maxtime=0.005, double eps=1e-5,
           sns_ik::SolveType type=sns_ik::Velocity);

    SNS_IK(const KDL::Chain& chain,
           const KDL::JntArray& q_min, const KDL::JntArray& q_max,
           const KDL::JntArray& v_max, const KDL::JntArray& a_max,
           double maxtime=0.005, double eps=1e-5,
           sns_ik::SolveType type=sns_ik::Velocity);

    ~SNS_IK();

    bool getKDLChain(KDL::Chain& chain) {
      chain=m_chain;
      return m_initialized;
    }

    bool getKDLLimits(KDL::JntArray& lb, KDL::JntArray& ub, KDL::JntArray& vel, KDL::JntArray& accel) {
      lb=m_lower_bounds;
      ub=m_upper_bounds;
      vel=m_velocity;
      accel=m_acceleration;
      return m_initialized;
    }

    int CartToJnt(const KDL::JntArray &q_init, const KDL::Frame &p_in, KDL::JntArray &q_out, const KDL::Twist& bounds=KDL::Twist::Zero());
    int CartToJnt(const KDL::JntArray& q_in, const KDL::Twist& v_in, KDL::JntArray& qdot_out);

    inline void SetSolveType(SolveType type) {
      m_solvetype = type;
    }


  private:
    bool m_initialized;
    double m_eps;
    double m_maxtime;
    SolveType m_solvetype;
    KDL::Chain m_chain;
    KDL::JntArray m_lower_bounds, m_upper_bounds, m_velocity, m_acceleration;
    enum JointType { Revolute, Prismatic, Continuous };
    std::vector<JointType> m_types;
    std::vector<KDL::JntArray> m_solutions;
    std::shared_ptr<SNSVelocityIK> m_ik_vel_solver;
    std::shared_ptr<SNSPositionIK> m_ik_pos_solver;
    std::shared_ptr<KDL::ChainJntToJacSolver> m_jacobianSolver;
    void initialize();

  };
}  //namespace
#endif
