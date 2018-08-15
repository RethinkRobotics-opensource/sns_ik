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
#include <sns_ik/sns_ik.hpp>
#include <sns_ik/sns_vel_ik_base.hpp>
#include <ros/ros.h>
#include <kdl/chainiksolverpos_nr_jl.hpp>
#include <kdl/chainfksolvervel_recursive.hpp>
#include <kdl/chainiksolvervel_pinv_nso.hpp>
#include <kdl/chainiksolvervel_pinv.hpp>
#include <kdl/framevel.hpp>
#include <time.h>
// Set USE_TRAC to 1 to test against trac_ik
#define USE_TRAC 0
#if USE_TRAC
  #include <trac_ik/trac_ik.hpp>
#endif

double fRand(double min, double max)
{
  double f = (double)rand() / RAND_MAX;
  return min + f * (max - min);
}

double getDeltaWithLimits(double value, double desired_delta,
                          double limit_min, double limit_max)
{
  // Validate delta is not more than half of joint limits
  double input_delta = std::min(desired_delta, (limit_max-limit_min)/2.0);
  double lower_delta = value-input_delta;
  double upper_delta = value+input_delta;
  // Randomize between upper delta and lower delta
  if (rand() % 2) {
    return upper_delta > limit_max ? lower_delta : upper_delta;
  } else {
    return lower_delta < limit_min ? upper_delta : lower_delta;
  }
}

bool in_vel_bounds(const KDL::JntArray& vel_values, const KDL::JntArray& vel_limits)
{
  for (int i = 0; i < vel_limits.data.size(); i++) {
    if (vel_values(i) < -vel_limits(i)-1e-6 || vel_values(i) > vel_limits(i)+1e-6) {
      return false;
    }
  }
  return true;
}

bool in_pos_bounds(const KDL::JntArray& jnt_values, const KDL::JntArray& jnt_lower,
                   const KDL::JntArray& jnt_upper)
{
  for (int i = 0; i < jnt_values.data.size(); i++){
    if (jnt_values(i) < jnt_lower(i)-1e-6 || jnt_values(i) > jnt_upper(i)+1e-6) {
      return false;
    }
  }
  return true;
}

double nullspace_l2_norm_ratio(const KDL::JntArray& jnt_result, const KDL::JntArray& jnt_ns_result,
                           const KDL::JntArray&  jnt_ns_bias)
{
  double error = 0;
  double ns_error = 0;
  for (int i = 0; i < jnt_ns_bias.data.size(); i++) {
    error += std::pow(jnt_result(i) - jnt_ns_bias(i), 2);
    ns_error += std::pow(jnt_ns_result(i) - jnt_ns_bias(i), 2);
  }
  return std::pow(ns_error, 0.5) / std::pow(error, 0.5);
}

// Compares linear and rotational velocities to see if they are they are scaled properly
bool velocityIsScaled(KDL::FrameVel fv1, KDL::FrameVel fv2, double eps, double *scale)
{
  KDL::Vector v1 = fv1.p.v;
  KDL::Vector v2 = fv2.p.v;
  double v1norm = v1.Norm();
  double v2norm = v2.Norm();
  *scale = v2norm / v1norm;

  // calculate inner product of the velocity vectors
  // theta = acos(cosTheta)
  double cosTheta = KDL::dot(v1, v2) / (v1norm * v2norm);

  // small angle approximation: std::cos(eps) ~ 1 - eps^2/2
  if (cosTheta < 1 - eps*eps/2) {
    return false; // linear velocity not scaled correctly
  }

  // compare rotation scaling
  KDL::Vector w1 = fv1.M.w;
  KDL::Vector w2 = fv2.M.w;
  double w1norm = w1.Norm();
  double w2norm = w2.Norm();
  double rotScale = w2norm / w1norm;

  // Check if the rotation scale matches the linear scale
  if (w1norm > eps && std::fabs(*scale-rotScale) > eps) {
    return false;
  }

  // calculate inner product of the roational velocity vectors
  double cosThetaW = KDL::dot(w1, w2) / (w1norm * w2norm);
  if (cosThetaW < 1 - eps*eps/2) {
    return false; // rotational velocity not scaled correctly
  }

  return true;
}

double standardDeviation(std::vector<double> values, double mean){
  double val_sum = 0;
  for(double sample: values){
    val_sum += std::pow((sample - mean), 2);
  }
  return std::pow(val_sum/double(values.size()), 0.5);
}


void test(ros::NodeHandle& nh, double num_samples_pos, double num_samples_vel,
          std::string chain_start, std::string chain_end, double timeout, double loop_period,
          std::string urdf_param, bool use_random_position_seed,
          bool use_delta_position_seed, double delta_position_seed_value,
          bool use_nullspace_bias_task, double nullspace_gain,
          bool nominal_nullspace, bool delta_nullspace, double delta_nullspace_value)
{

  double eps = 1e-5;

  // This constructor parses the URDF loaded in rosparm urdf_param into the
  // needed KDL structures.  We then pull these out to compare against the KDL
  // IK solver.
  KDL::Chain chain;
  KDL::JntArray ll, ul, vl, al; //lower joint limits, upper joint limits
  bool valid = false;
  sns_ik::SNS_IK snsik_solver(chain_start, chain_end, urdf_param, loop_period, eps, sns_ik::SNS);
  valid = snsik_solver.getKDLChain(chain);
  if (!valid) {
    ROS_ERROR("There was no valid KDL chain found");
    return;
  }
  valid = snsik_solver.getKDLLimits(ll,ul,vl,al);
  if (!valid) {
    ROS_ERROR("There were no valid KDL joint limits found");
    return;
  }
  assert(chain.getNrOfJoints() == ll.data.size());
  assert(chain.getNrOfJoints() == ul.data.size());
  assert(chain.getNrOfJoints() == vl.data.size());
  assert(chain.getNrOfJoints() == al.data.size());

  // Create Nominal chain configuration midway between all joint limits
  KDL::JntArray nominal(chain.getNrOfJoints());

  for (uint j=0; j<nominal.data.size(); j++) {
    nominal(j) = (ll(j)+ul(j))/2.0;
  }

  // Set up KDL IK
  KDL::ChainFkSolverPos_recursive fk_solver(chain); // Forward kin. solver
  KDL::ChainIkSolverVel_pinv vik_solver(chain); // PseudoInverse vel solver
  KDL::ChainIkSolverPos_NR_JL kdl_solver(chain, ll, ul, fk_solver, vik_solver, 1, eps); // Joint Limit Solver
  // 1 iteration per solve (will wrap in timed loop to compare with TRAC-IK)

  // PseudoInverse vel solver with nullspace optimization
  KDL::ChainIkSolverVel_pinv_nso vik_nso_solver(chain);
  vik_nso_solver.setOptPos(nominal);
  vik_nso_solver.setAlpha(nullspace_gain);
  KDL::JntArray nsWeights(chain.getNrOfJoints());
  nsWeights.data = Eigen::VectorXd::Ones(chain.getNrOfJoints());
  vik_nso_solver.setWeights(nsWeights);
  KDL::ChainIkSolverPos_NR_JL kdl_nso_solver(chain, ll, ul, fk_solver, vik_nso_solver, 1, eps);

  // Create desired number of valid, random joint configurations
  std::vector<KDL::JntArray> JointList;
  std::vector<KDL::JntArray> JointSeed;
  std::vector<KDL::JntArray> NullSpaceBias;
  KDL::JntArray q(chain.getNrOfJoints());
  KDL::JntArray q_delta(chain.getNrOfJoints());
  KDL::JntArray q_nullspace(chain.getNrOfJoints());

  uint num_joint_pos = std::max(num_samples_pos, num_samples_vel);
  for (uint i=0; i < num_joint_pos; i++) {
    for (uint j=0; j<ll.data.size(); j++) {
      q(j)=fRand(ll(j), ul(j));
    }
    JointList.push_back(q);
    // Set joint seed
    if (use_delta_position_seed){
      for (uint j=0; j<ll.data.size(); j++) {
        q_delta(j)=getDeltaWithLimits(q(j), delta_position_seed_value, ll(j), ul(j));
      }
      JointSeed.push_back(q_delta);
    } else if (i == 0 || !use_random_position_seed){
      JointSeed.push_back(nominal);
    } else { // "random" seed
      JointSeed.push_back(JointList[i-1]);
    }
    // Determine the NullSpace
    if (use_nullspace_bias_task) {
      if(delta_nullspace){
        for (uint j=0; j<ll.data.size(); j++) {
          q_nullspace(j)=getDeltaWithLimits(q(j), delta_nullspace_value, ll(j), ul(j));
        }
      } else if (nominal_nullspace){
        q_nullspace=nominal;
      } else { // random nullspace
        for (uint j=0; j<ll.data.size(); j++) {
          q_nullspace(j)=fRand(ll(j), ul(j));
        }
      }
      NullSpaceBias.push_back(q_nullspace);
    }
  }

  boost::posix_time::ptime start_time;
  boost::posix_time::time_duration diff;

  KDL::JntArray result;
  KDL::JntArray ns_result;
  KDL::Frame end_effector_pose;
  KDL::Frame end_effector_pose_check;
  int rc, ns_rc;
  int both_success_cnt=0;
  uint both_success=0;
  double total_time=0;
  double ns_total_time=0;
  double elapsed = 0;
  uint success=0;
  uint ns_success=0;
  double total_ns_l2_norm_ratio=0.0;

  std::vector<double> kdlPos_indivTime;

  ROS_INFO_STREAM("*** Testing KDL with "<<num_samples_pos<<" random samples");

  for (uint i=0; i < num_samples_pos; i++) {
    elapsed = 0;
    both_success=0;

    // Solve End effector pose for both natural and nullspace
    fk_solver.JntToCart(JointList[i], end_effector_pose);

    result = JointSeed[i];
    int cnt = 0;  // keep track of iteration count to enforce max number of iterations
    start_time = boost::posix_time::microsec_clock::local_time();
    do {
      q = result; // when iterating start with last solution
      rc = kdl_solver.CartToJnt(q, end_effector_pose, result);
      diff = boost::posix_time::microsec_clock::local_time() - start_time;
      elapsed = diff.total_nanoseconds() / 1e9;
    } while (rc < 0 && elapsed < timeout && cnt++ < 100);
    total_time += elapsed;
    kdlPos_indivTime.push_back(elapsed);
    fk_solver.JntToCart(result, end_effector_pose_check);
    bool inPosBounds = in_pos_bounds(result, ll, ul);
    if (rc>=0 && inPosBounds && Equal(end_effector_pose, end_effector_pose_check, 1e-3)) {
      success++;
      both_success++;
    }

    // Compare previous Inverse Kinematics calls against that use the nullspace
    if(use_nullspace_bias_task) {
      cnt = 0;
      elapsed = 0;
      start_time = boost::posix_time::microsec_clock::local_time();
      ns_result = JointSeed[i];
      do {
        q = ns_result; // when iterating start with last solution
        vik_nso_solver.setOptPos(NullSpaceBias[i]);
        ns_rc = kdl_nso_solver.CartToJnt(q, end_effector_pose, ns_result);
        diff = boost::posix_time::microsec_clock::local_time() - start_time;
        elapsed = diff.total_nanoseconds() / 1e9;
      } while (ns_rc < 0 && elapsed < timeout && cnt++ < 100);
      ns_total_time += elapsed;
      fk_solver.JntToCart(ns_result, end_effector_pose_check);
      if (ns_rc>=0 && in_pos_bounds(ns_result, ll, ul)
                && Equal(end_effector_pose, end_effector_pose_check, 1e-3)) {
        ns_success++;
        both_success++;
      }
      if(both_success==2) {
        total_ns_l2_norm_ratio += nullspace_l2_norm_ratio(result, ns_result, NullSpaceBias[i]);
        both_success_cnt++;
      }
    }

    if (int((double)i/num_samples_pos*100)%10 == 0)
      ROS_INFO_STREAM_THROTTLE(1,int((i)/num_samples_pos*100)<<"\% done");
  }
  double kdlPos_ns_avgTime = ns_total_time/num_samples_pos;
  double kdlPos_ns_avgL2Score = total_ns_l2_norm_ratio/both_success_cnt;
  double kdlPos_ns_successRate = ns_success/num_samples_pos;
  double kdlPos_successRate = success/num_samples_pos;
  double kdlPos_avgTime = total_time/num_samples_pos;
  double kdlPos_stdDev = standardDeviation(kdlPos_indivTime, kdlPos_avgTime);
  ROS_INFO_STREAM("KDL found " << success << " solutions (" << 100.0 * kdlPos_successRate
                  <<"\%) with an average of " << kdlPos_avgTime
                  << " secs per sample");
  ROS_INFO_STREAM("KDL nullspace success rate: "
                   << 100.0 * kdlPos_ns_successRate << ", avg l2_norm ratio: "
                   << kdlPos_ns_avgL2Score << ", avg time: " << kdlPos_ns_avgTime);
  #if USE_TRAC
    total_time=0;
    success=0;
    std::vector<double> tracPos_indivTime;
    TRAC_IK::TRAC_IK tracik_solver(chain_start, chain_end, urdf_param, timeout, eps);
    valid = tracik_solver.getKDLChain(chain);
    if (!valid) {
      ROS_ERROR("There was no valid KDL chain found");
      return;
    }

    ROS_INFO_STREAM("*** Testing TRAC-IK with "<<num_samples_pos<<" random samples");
    if (use_nullspace_bias_task)
      ROS_WARN("TRAC-IK does not support secondary tasks.");

    for (uint i=0; i < num_samples_pos; i++) {
      fk_solver.JntToCart(JointList[i],end_effector_pose);
      double elapsed = 0;
      start_time = boost::posix_time::microsec_clock::local_time();
      rc=tracik_solver.CartToJnt(JointSeed[i],end_effector_pose,result);
      diff = boost::posix_time::microsec_clock::local_time() - start_time;
      elapsed = diff.total_nanoseconds() / 1e9;
      total_time+=elapsed;
      tracPos_indivTime.push_back(elapsed);
      fk_solver.JntToCart(result, end_effector_pose_check);
      bool inPosBounds = in_pos_bounds(result, ll, ul);
      if (rc>=0 && inPosBounds && Equal(end_effector_pose, end_effector_pose_check, 1e-3))
        success++;

      if (int((double)i/num_samples_pos*100)%10 == 0)
        ROS_INFO_STREAM_THROTTLE(1,int((i)/num_samples_pos*100)<<"\% done");
    }

    double tracPos_successRate = success/num_samples_pos;
    double tracPos_avgTime = total_time/num_samples_pos;
    double tracPos_stdDev = standardDeviation(tracPos_indivTime, tracPos_avgTime);

    ROS_INFO_STREAM("TRAC-IK found " << success << " solutions (" << 100.0 * tracPos_successRate
                    << "\%) with an average of " << tracPos_avgTime << " secs per sample");
  #endif

  // SNS Position Tests
  snsik_solver.setNullspaceGain(nullspace_gain);
  struct velocitySolverData {
    sns_ik::VelocitySolveType type;
    std::string               name;
    double             successRate;
    double     scaling_successRate;
    double                avg_time;
    double    avg_ns_l2_norm_ratio;
    double          avg_ns_success;
    double             avg_ns_time;
    std::vector<double>  indiv_time;
    std::vector<double>  indiv_ns_time;
  };

  std::vector<velocitySolverData> vel_solver_data;
  velocitySolverData sns = {sns_ik::SNS,"SNS",0.0,0.0,0.0,0.0,0.0,0.0};
  vel_solver_data.push_back(sns);
  velocitySolverData sns_optimalsm = {sns_ik::SNS_OptimalScaleMargin,"SNS Optimal Scale Margin",0.0,0.0,0.0,0.0,0.0,0.0};
  vel_solver_data.push_back(sns_optimalsm);
  velocitySolverData sns_optimal = {sns_ik::SNS_Optimal,"SNS Optimal",0.0,0.0,0.0,0.0,0.0,0.0};
  vel_solver_data.push_back(sns_optimal);
  velocitySolverData sns_fast = {sns_ik::SNS_Fast,"SNS Fast",0.0,0.0,0.0,0.0,0.0,0.0};
  vel_solver_data.push_back(sns_fast);
  velocitySolverData sns_fastoptimal = {sns_ik::SNS_FastOptimal,"SNS Fast Optimal",0.0,0.0,0.0,0.0,0.0,0.0};
  vel_solver_data.push_back(sns_fastoptimal);
  velocitySolverData sns_base = {sns_ik::SNS_Base,"SNS Base",0.0,0.0,0.0,0.0,0.0,0.0};
  vel_solver_data.push_back(sns_base);

  for(auto& vst: vel_solver_data){
    snsik_solver.setVelocitySolveType(vst.type);
    // Initialize Solver Variables
    total_time=0;
    ns_total_time=0;
    success=0;
    ns_success=0;
    both_success_cnt=0;
    total_ns_l2_norm_ratio = 0.0;
    ROS_INFO_STREAM("*** Testing SNS-IK with "<<num_samples_pos<<" random samples");
    for (uint i=0; i < num_samples_pos; i++) {
      // Initialize Iteration Variables
      elapsed = 0;
      both_success = 0;
      //Solve forward kinematics for both
      fk_solver.JntToCart(JointList[i],end_effector_pose);
      start_time = boost::posix_time::microsec_clock::local_time();
      rc = snsik_solver.CartToJnt(JointSeed[i], end_effector_pose, result);
      diff = boost::posix_time::microsec_clock::local_time() - start_time;
      elapsed = diff.total_nanoseconds() / 1e9;
      total_time+=elapsed;
      vst.indiv_time.push_back(elapsed);
      fk_solver.JntToCart(result, end_effector_pose_check);
      if (rc>=0 && in_pos_bounds(result, ll, ul)
                && Equal(end_effector_pose, end_effector_pose_check, 1e-3)){
        success++;
        both_success++;
       }

      if(use_nullspace_bias_task) {
        elapsed = 0;
        start_time = boost::posix_time::microsec_clock::local_time();
        ns_rc = snsik_solver.CartToJnt(JointSeed[i], end_effector_pose, NullSpaceBias[i], ns_result);
        diff = boost::posix_time::microsec_clock::local_time() - start_time;
        elapsed = diff.total_nanoseconds() / 1e9;
        ns_total_time += elapsed;
        vst.indiv_ns_time.push_back(elapsed);
        fk_solver.JntToCart(ns_result, end_effector_pose_check);
        if (ns_rc>=0 && in_pos_bounds(ns_result, ll, ul)
                  && Equal(end_effector_pose, end_effector_pose_check, 1e-3)){
          ns_success++;
          both_success++;
        }
        if(both_success==2) {
          total_ns_l2_norm_ratio += nullspace_l2_norm_ratio(result, ns_result, NullSpaceBias[i]);
          both_success_cnt++;
        }
      }

      if (int((double)i/num_samples_pos*100)%10 == 0)
        ROS_INFO_STREAM_THROTTLE(1,int((i)/num_samples_pos*100)<<"\% done");
    }
    vst.successRate = success/num_samples_pos;
    vst.avg_time = total_time/num_samples_pos;
    ROS_INFO_STREAM(vst.name << " found " << success << " solutions ("
                    << 100*vst.successRate << "\%) with an average of " << vst.avg_time
                    << " secs per sample");
    if(use_nullspace_bias_task) {
        vst.avg_ns_l2_norm_ratio = total_ns_l2_norm_ratio/both_success_cnt;
        vst.avg_ns_success = ns_success/num_samples_pos;
        vst.avg_ns_time = ns_total_time/num_samples_pos;
        ROS_INFO_STREAM(vst.name <<" nullspace success rate: "<<vst.avg_ns_success
                    <<", avg l2_norm ratio: "<<vst.avg_ns_l2_norm_ratio
                    <<", avg time: "<<vst.avg_ns_time);
    }
  }

  ROS_INFO("\n************************************");
  ROS_INFO("Position IK Summary:");
  for(auto& vst: vel_solver_data){
      double std_dev = standardDeviation(vst.indiv_time, vst.avg_time);
      ROS_INFO("%s: %.2f%% success rate with (time: %.2f \u00b1 %.2f ms)",
               vst.name.c_str(), 100*vst.successRate, 1000*vst.avg_time, 1000*std_dev);
  }
  ROS_INFO("KDL: %.2f%% success rate with (time: %.2f \u00b1 %.2f ms)",
           100.*kdlPos_successRate, 1000*kdlPos_avgTime, 1000*kdlPos_stdDev);
  #if USE_TRAC
  ROS_INFO("TRAC: %.2f%% success rate with (time: %.2f \u00b1 %.2f ms)",
           100.*tracPos_successRate, 1000*tracPos_avgTime, 1000*tracPos_stdDev);
  #endif
  ROS_INFO("\n************************************\n");

  // Create random velocities within the limits
  uint successWithScaling;
  KDL::JntArrayVel result_vel_array;
  KDL::JntArray result_vel;
  KDL::FrameVel end_effector_vel;
  KDL::FrameVel result_end_effector_vel;
  KDL::ChainFkSolverVel_recursive vfk_solver(chain); // Foward kin. vel solver
  std::vector<KDL::JntArrayVel> JointVelList;
  KDL::JntArrayVel v(chain.getNrOfJoints());
  // Create Random velocities within the vel limits
  for (uint i=0; i < num_samples_vel; i++) {
    v.q = JointList[i];
    for (uint j=0; j<vl.data.size(); j++) {
      v.qdot(j)=fRand(-vl(j), vl(j));
    }
    JointVelList.push_back(v);
  }

  for(auto& vst: vel_solver_data){
    total_time=0;
    success=0;
    successWithScaling = 0;
    //ROS_INFO_STREAM("*** Testing SNS-IK Velocity: "<<vst.name);
    snsik_solver.setVelocitySolveType(vst.type);
    for (uint i=0; i < num_samples_vel; i++) {
      // add position to my vel
      vfk_solver.JntToCart(JointVelList[i],end_effector_vel);
      double elapsed = 0;
      start_time = boost::posix_time::microsec_clock::local_time();
      if(use_nullspace_bias_task)
        rc = snsik_solver.CartToJntVel(JointVelList[i].q, end_effector_vel.GetTwist(), NullSpaceBias[i], result_vel);
      else
        rc = snsik_solver.CartToJntVel(JointVelList[i].q, end_effector_vel.GetTwist(), result_vel);

      diff = boost::posix_time::microsec_clock::local_time() - start_time;
      elapsed = diff.total_nanoseconds() / 1e9;
      total_time+=elapsed;

      // check to make sure vel is within limit
      result_vel_array.q = JointVelList[i].q;
      result_vel_array.qdot = result_vel;
      vfk_solver.JntToCart(result_vel_array, result_end_effector_vel);
      bool inVelBounds = in_vel_bounds(result_vel, vl);
      if (rc>=0 && inVelBounds && Equal(end_effector_vel, result_end_effector_vel, 1e-3))
        success++;
      double scale;
      if (rc>=0 && velocityIsScaled(end_effector_vel, result_end_effector_vel, 1e-3, &scale) && inVelBounds)
        successWithScaling++;

      /*if (int((double)i/num_samples_vel*100)%10 == 0)
        ROS_INFO_STREAM_THROTTLE(1,int((i)/num_samples_vel*100)<<"\% done");*/
    }
    //ROS_INFO_STREAM("Velocities Solver found "<<success<<" solutions ("<<100.0*success/num_samples_vel<<"\%) with an average of "<<total_time/num_samples_vel<<" secs per sample");
    //ROS_INFO_STREAM("Velocity Scaling found " << successWithScaling << " solutions ("
    //            << 100.0*successWithScaling/num_samples_vel << "\%)");
    vst.successRate = success/num_samples_vel;
    vst.scaling_successRate = successWithScaling/num_samples_vel;
    vst.avg_time = total_time/num_samples_vel;
  }

  ROS_INFO_STREAM("*** Testing KDL-IK Velocities with " << num_samples_vel << " random samples");
  total_time = 0;
  success = 0;
  successWithScaling = 0;
  for (uint i = 0; i < num_samples_vel; i++) {
    // add position to my vel
    vfk_solver.JntToCart(JointVelList[i], end_effector_vel);
    double elapsed = 0;
    start_time = boost::posix_time::microsec_clock::local_time();

    if (use_nullspace_bias_task) {
      vik_nso_solver.setOptPos(NullSpaceBias[i]);
      rc = vik_nso_solver.CartToJnt(JointVelList[i].q, end_effector_vel.GetTwist(), result_vel);
    } else {
      rc = vik_solver.CartToJnt(JointVelList[i].q, end_effector_vel.GetTwist(), result_vel);
    }
    diff = boost::posix_time::microsec_clock::local_time() - start_time;
    elapsed = diff.total_nanoseconds() / 1e9;
    total_time+=elapsed;
    // check to make sure vel is within limit
    result_vel_array.q = JointVelList[i].q;
    result_vel_array.qdot = result_vel;
    vfk_solver.JntToCart(result_vel_array, result_end_effector_vel);
    bool inVelBounds = in_vel_bounds(result_vel, vl);
    if (rc>=0 && inVelBounds && Equal(end_effector_vel, result_end_effector_vel, 1e-3))
      success++;
    double scale;
    if (rc>=0 && velocityIsScaled(end_effector_vel, result_end_effector_vel, 1e-3, &scale) && inVelBounds)
      successWithScaling++;

    if (int((double)i/num_samples_vel*100)%10 == 0)
      ROS_INFO_STREAM_THROTTLE(1,int((i)/num_samples_vel*100)<<"\% done");
  }

  double kdlVel_successRate = success/num_samples_vel;
  double kdlVel_scalingSuccessRate = successWithScaling/num_samples_vel;
  double kdlVel_avgTime = total_time/num_samples_vel;
  ROS_INFO_STREAM("KDL Velocity found " << success << " solutions (" << 100.0 * kdlVel_successRate
                    << "\%) with an average of " << kdlVel_avgTime << " secs per sample");
  ROS_INFO_STREAM("KDL Velocity Scaling Score " << successWithScaling << " solutions ("
                   << 100.0*kdlVel_scalingSuccessRate << "\%)");

  ROS_INFO("\n************************************");
  ROS_INFO("Velocity IK Summary:");
  for(auto& vst: vel_solver_data){
      ROS_INFO("%s: %.2f%% w/o and %.2f%% w/ scaling success rates (%.3f ms)",
               vst.name.c_str(), 100*vst.successRate, 100*vst.scaling_successRate, 1000*vst.avg_time);
  }
  ROS_INFO("KDL Velocity: %.2f%% w/o and %.2f%% w/ scaling success rates (%.3f ms)",
           100*kdlVel_successRate, 100*kdlVel_scalingSuccessRate, 1000*kdlVel_avgTime);
  ROS_INFO("\n************************************");

}

int main(int argc, char** argv)
{
  srand(time(NULL));
  ros::init(argc, argv, "ik_tests");
  ros::NodeHandle nh("~");

  setlocale(LC_CTYPE,"");  // needed to print ±

  int num_samples_pos, num_samples_vel;
  std::string chain_start, chain_end, urdf_param;
  double timeout, loop_period;
  bool use_random_position_seed;
  bool use_delta_position_seed;
  double delta_position_seed_value;
  bool use_nullspace_bias_task;
  double nullspace_gain;
  bool nominal_nullspace, delta_nullspace;
  double delta_nullspace_value;
  nh.param("num_samples_pos", num_samples_pos, 100);
  nh.param("num_samples_vel", num_samples_vel, 1000);
  nh.param("chain_start", chain_start, std::string(""));
  nh.param("chain_end", chain_end, std::string(""));
  nh.param("use_random_position_seed", use_random_position_seed, false);
  nh.param("use_delta_position_seed", use_delta_position_seed, false);
  nh.param("delta_position_seed_value", delta_position_seed_value, 0.2);
  nh.param("use_nullspace_bias_task", use_nullspace_bias_task, false);
  nh.param("nullspace_gain", nullspace_gain, 0.3);
  nh.param("nominal_nullspace", nominal_nullspace, false);
  nh.param("delta_nullspace", delta_nullspace, false);
  nh.param("delta_nullspace_value", delta_nullspace_value, 0.1);

  if (chain_start=="" || chain_end=="") {
    ROS_FATAL("Missing chain info in launch file");
    exit (-1);
  }

  nh.param("timeout", timeout, 0.005);
  nh.param("loop_period", loop_period, 0.005);
  nh.param("urdf_param", urdf_param, std::string("/robot_description"));

  if (num_samples_pos < 1)
    num_samples_pos = 1;

  test(nh, num_samples_pos, num_samples_vel,
       chain_start, chain_end, timeout,loop_period,
       urdf_param, use_random_position_seed,
       use_delta_position_seed, delta_position_seed_value,
       use_nullspace_bias_task, nullspace_gain,
       nominal_nullspace, delta_nullspace, delta_nullspace_value);

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
