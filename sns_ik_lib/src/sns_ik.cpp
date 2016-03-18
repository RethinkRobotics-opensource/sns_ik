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

#include <sns_ik/sns_ik.hpp>
#include <ros/ros.h>
#include <kdl_parser/kdl_parser.hpp>
#include <urdf/model.h>
#include <sns_ik/sns_velocity_ik.hpp>
#include <sns_ik/osns_velocity_ik.hpp>
#include <sns_ik/osns_sm_velocity_ik.hpp>
#include <sns_ik/fsns_velocity_ik.hpp>
#include <sns_ik/fosns_velocity_ik.hpp>

namespace sns_ik {

  SNS_IK::SNS_IK(const std::string& base_link, const std::string& tip_link,
                 const std::string& URDF_param, double looprate, double eps,
                 sns_ik::VelocitySolveType type) :
    m_initialized(false),
    m_looprate(looprate),
    m_eps(eps),
    m_solvetype(type)
  {
    ros::NodeHandle node_handle("~");
    urdf::Model robot_model;
    std::string xml_string;
    std::string urdf_xml, full_urdf_xml;
    node_handle.param("urdf_param",urdf_xml,URDF_param);
    node_handle.searchParam(urdf_xml,full_urdf_xml);

    ROS_DEBUG_NAMED("sns_ik","Reading xml file from parameter server");
    if (!node_handle.getParam(full_urdf_xml, xml_string)) {
      ROS_FATAL_NAMED("sns_ik","Could not load the xml from parameter server: %s", urdf_xml.c_str());
      return;
    }

    node_handle.param(full_urdf_xml, xml_string, std::string());
    robot_model.initString(xml_string);

    ROS_DEBUG_STREAM_NAMED("sns_ik","Reading joints and links from URDF");
    KDL::Tree tree;
    if (!kdl_parser::treeFromUrdfModel(robot_model, tree)) {
      ROS_FATAL("Failed to extract kdl tree from xml robot description.");
      return;
    }

    if(!tree.getChain(base_link, tip_link, m_chain)) {
      ROS_FATAL("Couldn't find chain %s to %s",base_link.c_str(),tip_link.c_str());
    }

    std::vector<KDL::Segment> chain_segments = m_chain.segments;
    boost::shared_ptr<const urdf::Joint> joint;
    m_lower_bounds.resize(m_chain.getNrOfJoints());
    m_upper_bounds.resize(m_chain.getNrOfJoints());
    m_velocity.resize(m_chain.getNrOfJoints());
    m_acceleration.resize(m_chain.getNrOfJoints());

    unsigned int joint_num=0;
    for(std::size_t i = 0; i < chain_segments.size(); ++i) {
      joint = robot_model.getJoint(chain_segments[i].getJoint().getName());
      if (joint->type != urdf::Joint::UNKNOWN && joint->type != urdf::Joint::FIXED) {
        double lower=0; //TODO Better default values? Error if these arent found?
        double upper=0;
        double velocity=0;
        double acceleration=0;

        if ( joint->type == urdf::Joint::CONTINUOUS ) {
            lower=std::numeric_limits<float>::lowest();
            upper=std::numeric_limits<float>::max();
        } else {
          if(joint->safety) {
            lower = std::max(joint->limits->lower, joint->safety->soft_lower_limit);
            upper = std::min(joint->limits->upper, joint->safety->soft_upper_limit);
          } else {
            lower = joint->limits->lower;
            upper = joint->limits->upper;
          }
          velocity = std::abs(joint->limits->velocity);
        }
        // Checking the Param server for limit modifications
        // and acceleration limits
        std::string prefix = urdf_xml + "_planning/joint_limits/" + joint->name + "/";
        double ul;
        if(node_handle.getParam(prefix + "max_position", ul)){
          upper = std::min(upper, ul);
        }
        double ll;
        if(node_handle.getParam(prefix + "min_position", ll)){
          lower = std::max(lower, ll);
        }
        double vel;
        if(node_handle.getParam(prefix + "max_velocity", vel)){
          velocity = std::min(velocity, std::abs(vel));
        }
        node_handle.getParam(prefix + "max_acceleration", acceleration);

        m_lower_bounds(joint_num)=lower;
        m_upper_bounds(joint_num)=upper;
        m_velocity(joint_num) = velocity;
        m_acceleration(joint_num) = std::abs(acceleration);

        ROS_INFO_STREAM("sns_ik Using joint "<<joint->name<<" lb:"
                        <<m_lower_bounds(joint_num)<<" ub:"<<m_upper_bounds(joint_num)
                        << " v:"<<m_velocity(joint_num)<<" a:"<<m_acceleration(joint_num));
        joint_num++;
      }
    }
    initialize();
  }


  SNS_IK::SNS_IK(const KDL::Chain& chain, const KDL::JntArray& q_min,
                 const KDL::JntArray& q_max, const KDL::JntArray& v_max,
                 const KDL::JntArray& a_max, double looprate, double eps,
                 sns_ik::VelocitySolveType type):
    m_initialized(false),
    m_looprate(looprate),
    m_eps(eps),
    m_solvetype(type),
    m_chain(chain),
    m_lower_bounds(q_min),
    m_upper_bounds(q_max),
    m_velocity(v_max),
    m_acceleration(a_max)
  {
    initialize();
  }

  void SNS_IK::initialize() {

    ROS_ASSERT_MSG(m_chain.getNrOfJoints()==m_lower_bounds.data.size(), \
                "SNS_IK: Number of joint lower bounds does not equal number of joints");
    ROS_ASSERT_MSG(m_chain.getNrOfJoints()==m_upper_bounds.data.size(), \
                "SNS_IK: Number of joint upper bounds does not equal number of joints");
    ROS_ASSERT_MSG(m_chain.getNrOfJoints()==m_velocity.data.size(), \
                "SNS_IK: Number of max joint velocity bounds does not equal number of joints");
    ROS_ASSERT_MSG(m_chain.getNrOfJoints()==m_acceleration.data.size(), \
                "SNS_IK: Number of max joint acceleration bounds does not equal number of joints");
    // Populate a vector cooresponding to the type for every joint
    for (std::size_t i=0; i<m_chain.segments.size(); i++) {
      std::string type = m_chain.segments[i].getJoint().getTypeName();
      if (type.find("Rot")!=std::string::npos) {
        if (m_upper_bounds(m_types.size())>=std::numeric_limits<float>::max() &&
            m_lower_bounds(m_types.size())<=std::numeric_limits<float>::lowest()){
          m_types.push_back(SNS_IK::JointType::Continuous);
        } else {
          m_types.push_back(SNS_IK::JointType::Revolute);
        }
      } else if (type.find("Trans")!=std::string::npos) {
        m_types.push_back(SNS_IK::JointType::Prismatic);
      }
    }
    ROS_ASSERT_MSG(m_types.size()==(unsigned int)m_lower_bounds.data.size(),
                   "SNS_IK: Could not determine joint limits for all non-continuous joints");

    m_jacobianSolver = std::shared_ptr<KDL::ChainJntToJacSolver>(new KDL::ChainJntToJacSolver(m_chain));
    ROS_ASSERT_MSG(setVelocitySolveType(m_solvetype),
                   "SNS_IK: Failed to create a new SNS velocity and position solver."); //TODO make loop rate configurable
  }

bool SNS_IK::setVelocitySolveType(VelocitySolveType type) {
  // If the requested solve type is different or there is no velocity solver
  if(m_solvetype != type || !m_ik_vel_solver){
    switch (type) {
      case sns_ik::SNS_OptimalScaleMargin:
        m_ik_vel_solver = std::shared_ptr<OSNS_sm_VelocityIK>(new OSNS_sm_VelocityIK(m_chain.getNrOfJoints(), m_looprate));
        ROS_INFO("SNS_IK: Set Velocity solver to SNS Optimal Scale Margin solver.");
        break;
      case sns_ik::SNS_Optimal:
        m_ik_vel_solver = std::shared_ptr<OSNSVelocityIK>(new OSNSVelocityIK(m_chain.getNrOfJoints(), m_looprate));
        ROS_INFO("SNS_IK: Set Velocity solver to SNS Optimal solver.");
        break;
      case sns_ik::SNS_Fast:
        m_ik_vel_solver = std::shared_ptr<FSNSVelocityIK>(new FSNSVelocityIK(m_chain.getNrOfJoints(), m_looprate));
        ROS_INFO("SNS_IK: Set Velocity solver to Fast SNS solver.");
        break;
      case sns_ik::SNS_FastOptimal:
        m_ik_vel_solver = std::shared_ptr<FOSNSVelocityIK>(new FOSNSVelocityIK(m_chain.getNrOfJoints(), m_looprate));
        ROS_INFO("SNS_IK: Set Velocity solver to Fast Optimal SNS solver.");
        break;
      case sns_ik::SNS:
        m_ik_vel_solver = std::shared_ptr<SNSVelocityIK>(new SNSVelocityIK(m_chain.getNrOfJoints(), m_looprate));
        ROS_INFO("SNS_IK: Set Velocity solver to Standard SNS solver.");
        break;
      default:
        ROS_ERROR("SNS_IK: Unknow Velocity solver type requested.");
        return false;
    }
    m_ik_vel_solver->setJointsCapabilities(m_lower_bounds.data, m_upper_bounds.data,
                                           m_velocity.data, m_acceleration.data);
    m_ik_pos_solver = std::shared_ptr<SNSPositionIK>(new SNSPositionIK(m_chain, m_ik_vel_solver));
    m_solvetype = type;
    m_initialized = true;
    return true;
 }
 return false;
}

int SNS_IK::CartToJnt(const KDL::JntArray &q_init, const KDL::Frame &p_in, KDL::JntArray &q_out, const KDL::Twist& tolerances) {

  if (!m_initialized) {
    ROS_ERROR("SNS_IK was not properly initialized with a valid chain or limits.");
    return -1;
  }
  return m_ik_pos_solver->CartToJnt(q_init, p_in, &q_out, tolerances);
}

int SNS_IK::CartToJnt(const KDL::JntArray& q_in, const KDL::Twist& v_in, KDL::JntArray& qdot_out) {

  if (!m_initialized) {
    ROS_ERROR("SNS_IK was not properly initialized with a valid chain or limits.");
    return -1;
  }
  KDL::Jacobian jacobian;
  jacobian.resize(q_in.rows());
  if (m_jacobianSolver->JntToJac(q_in, jacobian) < 0)
  {
    std::cout << "JntToJac failed" << std::endl;
    return -1;
  }
  StackOfTasks sot;
  Task task;
  task.jacobian = jacobian.data;
  task.desired = VectorD::Zero(6);
  // twistEigenToKDL
  for(size_t i = 0; i < 6; i++)
      task.desired(i) = v_in[i];
  sot.push_back(task);
  return m_ik_vel_solver->getJointVelocity(&qdot_out.data, sot, q_in.data);
}

SNS_IK::~SNS_IK(){}

} // sns_ik namespace
