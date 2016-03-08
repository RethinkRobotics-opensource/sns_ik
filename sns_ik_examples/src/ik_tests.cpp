/********************************************************************************
Copyright (c) 2016, Rethink Robotics, Inc.
Copyright (c) 2016, TRACLabs, Inc.
All rights reserved.

Redistribution and use in source and binary forms, with or without modification,
 are permitted provided that the following conditions are met:

    1. Redistributions of source code must retain the above copyright notice,
       this list of conditions and the following disclaimer.

    2. Redistributions in binary form must reproduce the above copyright notice,
       this list of conditions and the following disclaimer in the documentation
       and/or other materials provided with the distribution.

    3. Neither the name of the copyright holder nor the names of its contributors
       may be used to endorse or promote products derived from this software
       without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE
OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED
OF THE POSSIBILITY OF SUCH DAMAGE.
********************************************************************************/

#include <boost/date_time.hpp>
#include <trac_ik/trac_ik.hpp>
#include <sns_ik/sns_ik.hpp>
#include <ros/ros.h>
#include <kdl/chainiksolverpos_nr_jl.hpp>
#include <kdl/chainfksolvervel_recursive.hpp>
#include <kdl/framevel.hpp>

double fRand(double min, double max)
{
  double f = (double)rand() / RAND_MAX;
  return min + f * (max - min);
}

bool in_vel_bounds(KDL::JntArray vals, KDL::JntArray vels)
{
  for(size_t i; i < vels.data.size(); i++){
      if(vals(i) < -vels(i) || vals(i) > vels(i)){
          return false;
      }
      /*else{
          std::cout<<vals(i)<<' ';
      }*/
  }
  //std::cout<<std::endl;
  return true;
}

void test(ros::NodeHandle& nh, double num_samples, std::string chain_start, std::string chain_end, double timeout, std::string urdf_param)
{

  double eps = 1e-5;

  // This constructor parses the URDF loaded in rosparm urdf_param into the
  // needed KDL structures.  We then pull these out to compare against the KDL
  // IK solver.
  TRAC_IK::TRAC_IK tracik_solver(chain_start, chain_end, urdf_param, timeout, eps);

  KDL::Chain chain;
  KDL::JntArray ll, ul, vl, al; //lower joint limits, upper joint limits

  bool valid = tracik_solver.getKDLChain(chain);

  if (!valid) {
    ROS_ERROR("There was no valid KDL chain found");
    return;
  }

  valid = tracik_solver.getKDLLimits(ll,ul);

  if (!valid) {
    ROS_ERROR("There were no valid KDL joint limits found");
    return;
  }

  assert(chain.getNrOfJoints() == ll.data.size());
  assert(chain.getNrOfJoints() == ul.data.size());

  ROS_INFO ("Using %d joints",chain.getNrOfJoints());


  // Set up KDL IK
  KDL::ChainFkSolverPos_recursive fk_solver(chain); // Forward kin. solver
  KDL::ChainIkSolverVel_pinv vik_solver(chain); // PseudoInverse vel solver
  KDL::ChainIkSolverPos_NR_JL kdl_solver(chain,ll,ul,fk_solver, vik_solver, 1, eps); // Joint Limit Solver
  // 1 iteration per solve (will wrap in timed loop to compare with TRAC-IK)


  // Create Nominal chain configuration midway between all joint limits
  KDL::JntArray nominal(chain.getNrOfJoints());

  for (uint j=0; j<nominal.data.size(); j++) {
    nominal(j) = (ll(j)+ul(j))/2.0;
  }

  // Create desired number of valid, random joint configurations
  std::vector<KDL::JntArray> JointList;
  KDL::JntArray q(chain.getNrOfJoints());

  for (uint i=0; i < num_samples; i++) {
    for (uint j=0; j<ll.data.size(); j++) {
      q(j)=fRand(ll(j), ul(j));
    }
    JointList.push_back(q);
  }

  boost::posix_time::ptime start_time;
  boost::posix_time::time_duration diff;

  KDL::JntArray result;
  KDL::Frame end_effector_pose;
  int rc;

  double total_time=0;
  uint success=0;

  ROS_INFO_STREAM("*** Testing KDL with "<<num_samples<<" random samples");

  for (uint i=0; i < num_samples; i++) {
    fk_solver.JntToCart(JointList[i],end_effector_pose);
    double elapsed = 0;
    result=nominal; // start with nominal
    start_time = boost::posix_time::microsec_clock::local_time();
    do {
      q=result; // when iterating start with last solution
      rc=kdl_solver.CartToJnt(q,end_effector_pose,result);
      diff = boost::posix_time::microsec_clock::local_time() - start_time;
      elapsed = diff.total_nanoseconds() / 1e9;
    } while (rc < 0 && elapsed < timeout);
    total_time+=elapsed;
    if (rc>=0)
      success++;

    if (int((double)i/num_samples*100)%10 == 0)
      ROS_INFO_STREAM_THROTTLE(1,int((i)/num_samples*100)<<"\% done");
  }

  ROS_INFO_STREAM("KDL found "<<success<<" solutions ("<<100.0*success/num_samples<<"\%) with an average of "<<total_time/num_samples<<" secs per sample");


  total_time=0;
  success=0;

  ROS_INFO_STREAM("*** Testing TRAC-IK with "<<num_samples<<" random samples");

  for (uint i=0; i < num_samples; i++) {
    fk_solver.JntToCart(JointList[i],end_effector_pose);
    double elapsed = 0;
    start_time = boost::posix_time::microsec_clock::local_time();
    rc=tracik_solver.CartToJnt(nominal,end_effector_pose,result);
    diff = boost::posix_time::microsec_clock::local_time() - start_time;
    elapsed = diff.total_nanoseconds() / 1e9;
    total_time+=elapsed;
    if (rc>=0)
      success++;

    if (int((double)i/num_samples*100)%10 == 0)
      ROS_INFO_STREAM_THROTTLE(1,int((i)/num_samples*100)<<"\% done");
  }

  ROS_INFO_STREAM("TRAC-IK found "<<success<<" solutions ("<<100.0*success/num_samples<<"\%) with an average of "<<total_time/num_samples<<" secs per sample");


  sns_ik::SNS_IK snsik_solver(chain_start, chain_end, urdf_param, timeout, eps);
  valid = snsik_solver.getKDLChain(chain);
  if (!valid) {
    ROS_ERROR("SNS: There was no valid KDL chain found");
    return;
  }
  valid = snsik_solver.getKDLLimits(ll,ul,vl,al);
  if (!valid) {
    ROS_ERROR("SNS: There were no valid KDL joint limits found");
    return;
  }
  assert(chain.getNrOfJoints() == ll.data.size());
  assert(chain.getNrOfJoints() == ul.data.size());
  assert(chain.getNrOfJoints() == vl.data.size());
  assert(chain.getNrOfJoints() == al.data.size());
  /* Position Tests
  total_time=0;
  success=0;
  ROS_INFO_STREAM("*** Testing SNS-IK with "<<num_samples<<" random samples");
  for (uint i=0; i < num_samples; i++) {
    fk_solver.JntToCart(JointList[i],end_effector_pose);
    double elapsed = 0;
    start_time = boost::posix_time::microsec_clock::local_time();
    rc=snsik_solver.CartToJnt(nominal,end_effector_pose,result);
    diff = boost::posix_time::microsec_clock::local_time() - start_time;
    elapsed = diff.total_nanoseconds() / 1e9;
    total_time+=elapsed;
    if (rc>=0)
      success++;
    if (int((double)i/num_samples*100)%10 == 0)
      ROS_INFO_STREAM_THROTTLE(1,int((i)/num_samples*100)<<"\% done");
  }
  ROS_INFO_STREAM("SNS-IK found "<<success<<" solutions ("<<100.0*success/num_samples<<"\%) with an average of "<<total_time/num_samples<<" secs per sample");
  */
  // Create random velocities within the limits
  total_time=0;
  success=0;
  KDL::JntArrayVel result_vel_array;
  KDL::JntArray result_vel;
  KDL::FrameVel end_effector_vel;
  KDL::FrameVel result_end_effector_vel;
  KDL::ChainFkSolverVel_recursive vfk_solver(chain); // Foward kin. vel solver
  std::vector<KDL::JntArrayVel> JointVelList;
  KDL::JntArrayVel v(chain.getNrOfJoints());
  // Create Random velocities within the vel limits
  for (uint i=0; i < num_samples; i++) {
    v.q = JointList[i];
    for (uint j=0; j<vl.data.size(); j++) {
      v.qdot(j)=fRand(-vl(j), vl(j));
    }
    JointVelList.push_back(v);
  }
  ROS_INFO_STREAM("*** Testing SNS-IK Velocities with "<<num_samples<<" random samples");
  for (uint i=0; i < num_samples; i++) {
    // add position to my vel
    vfk_solver.JntToCart(JointVelList[i],end_effector_vel);
    double elapsed = 0;
    start_time = boost::posix_time::microsec_clock::local_time();
    rc=snsik_solver.CartToJnt(JointVelList[i].q, end_effector_vel.GetTwist(), result_vel);
    // check to make sure vel is within limit

    diff = boost::posix_time::microsec_clock::local_time() - start_time;
    elapsed = diff.total_nanoseconds() / 1e9;
    total_time+=elapsed;
    result_vel_array.q = JointVelList[i].q;
    result_vel_array.qdot = result_vel;
    vfk_solver.JntToCart(result_vel_array, result_end_effector_vel);
    if (rc>=0 && in_vel_bounds(result_vel, vl) && Equal(end_effector_vel, result_end_effector_vel, 1e-3))
      success++;
    if (int((double)i/num_samples*100)%10 == 0)
      ROS_INFO_STREAM_THROTTLE(1,int((i)/num_samples*100)<<"\% done");
  }
  ROS_INFO_STREAM("SNS-IK Velocities found "<<success<<" solutions ("<<100.0*success/num_samples<<"\%) with an average of "<<total_time/num_samples<<" secs per sample");

  ROS_INFO_STREAM("*** Testing KDL-IK Velocities with "<<num_samples<<" random samples");
  total_time=0;
  success=0;
  for (uint i=0; i < num_samples; i++) {
    // add position to my vel
    vfk_solver.JntToCart(JointVelList[i],end_effector_vel);
    double elapsed = 0;
    start_time = boost::posix_time::microsec_clock::local_time();
    rc=vik_solver.CartToJnt(JointVelList[i].q, end_effector_vel.GetTwist(), result_vel);
    // check to make sure vel is within limit

    diff = boost::posix_time::microsec_clock::local_time() - start_time;
    elapsed = diff.total_nanoseconds() / 1e9;
    total_time+=elapsed;
    result_vel_array.q = JointVelList[i].q;
    result_vel_array.qdot = result_vel;
    vfk_solver.JntToCart(result_vel_array, result_end_effector_vel);
    if (rc>=0 && in_vel_bounds(result_vel, vl) && Equal(end_effector_vel, result_end_effector_vel, 1e-3))
      success++;
    if (int((double)i/num_samples*100)%10 == 0)
      ROS_INFO_STREAM_THROTTLE(1,int((i)/num_samples*100)<<"\% done");
  }
  ROS_INFO_STREAM("KDL-IK Velocities found "<<success<<" solutions ("<<100.0*success/num_samples<<"\%) with an average of "<<total_time/num_samples<<" secs per sample");


}

int main(int argc, char** argv)
{
  srand(1);
  ros::init(argc, argv, "ik_tests");
  ros::NodeHandle nh("~");

  int num_samples;
  std::string chain_start, chain_end, urdf_param;
  double timeout;

  nh.param("num_samples", num_samples, 1000);
  nh.param("chain_start", chain_start, std::string(""));
  nh.param("chain_end", chain_end, std::string(""));

  if (chain_start=="" || chain_end=="") {
    ROS_FATAL("Missing chain info in launch file");
    exit (-1);
  }

  nh.param("timeout", timeout, 0.005);
  nh.param("urdf_param", urdf_param, std::string("/robot_description"));

  if (num_samples < 1)
    num_samples = 1;

  test(nh, num_samples, chain_start, chain_end, timeout, urdf_param);

  // Useful when you make a script that loops over multiple launch files that test different robot chains
  std::vector<char *> commandVector;
  commandVector.push_back((char*)"killall");
  commandVector.push_back((char*)"-9");
  commandVector.push_back((char*)"roslaunch");
  commandVector.push_back(NULL);

  char **command = &commandVector[0];
  execvp(command[0],command);

  return 0;
}
