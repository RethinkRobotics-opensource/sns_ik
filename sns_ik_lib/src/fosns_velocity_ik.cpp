/*! \file fosns_velocity_ik.cpp
 * \brief Fast Optimal SNS velocity IK solver
 * \author Fabrizio Flacco
 * \author Forrest Rogers-Marcovitz
 */
/*
 *    Copyright 2016 Rethink Robotics
 *
 *    Copyright 2012-2016 Fabrizio Flacco
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

#include <sns_ik/fosns_velocity_ik.hpp>

#include <ros/ros.h>
#include <iostream>

namespace sns_ik {

FOSNSVelocityIK::FOSNSVelocityIK(int dof, double loop_period) :
    FSNSVelocityIK(dof, loop_period),
    scaleMargin(0.98)
{
}

void FOSNSVelocityIK::setNumberOfTasks(int ntasks, int dof)
{
  SNSVelocityIK::setNumberOfTasks(ntasks, dof);

  //for the Fast version
  Eigen::MatrixXd Z = Eigen::MatrixXd::Zero(n_dof, n_dof);
  Eigen::VectorXd zv = Eigen::VectorXd::Zero(n_dof);
  Eigen::VectorXi zvi = Eigen::VectorXi::Zero(n_dof);
  B = Z;

  S.resize(n_tasks, zvi);
  nSat.resize(n_tasks, 0);

  satList.resize(n_tasks);
  lagrangeMu = zv;
  lagrangeMu1 = zv;
  lagrangeMup2w = zv;
}

double FOSNSVelocityIK::getJointVelocity(Eigen::VectorXd *jointVelocity,
    const std::vector<Task> &sot,
    const Eigen::VectorXd &jointConfiguration)
{
  // This will only reset member variables if different from previous values
  setNumberOfTasks(sot.size(), sot[0].jacobian.cols());
  S.resize(n_tasks, Eigen::VectorXi::Zero(n_dof));

  // TODO: check that setJointsCapabilities has been already called

  Eigen::MatrixXd P = Eigen::MatrixXd::Identity(n_dof, n_dof);
  Eigen::MatrixXd PS = Eigen::MatrixXd::Identity(n_dof, n_dof);
  *jointVelocity = Eigen::VectorXd::Zero(n_dof, 1);
  Eigen::VectorXd higherPriorityJointVelocity;
  Eigen::MatrixXd higherPriorityNull;

  shapeJointVelocityBound(jointConfiguration);

  // this is not the best solution... the scale margin should be computed inside FOSNSsingle

  for (int i_task = 0; i_task < n_tasks; i_task++) {  //consider all tasks
    higherPriorityJointVelocity = *jointVelocity;
    higherPriorityNull = P;
    scaleFactors[i_task] = SNSsingle(i_task, higherPriorityJointVelocity, higherPriorityNull,
        sot[i_task].jacobian, sot[i_task].desired, jointVelocity, &PS);

    if (scaleFactors[i_task] > 0.0) {
      if (scaleFactors[i_task] * scaleMargin < 1.0) {
        double taskScale = scaleFactors[i_task] * scaleMargin;
        Eigen::VectorXd scaledTask = sot[i_task].desired * taskScale;
        scaleFactors[i_task] = SNSsingle(i_task, higherPriorityJointVelocity, higherPriorityNull,
            sot[i_task].jacobian, scaledTask, jointVelocity, &P);
        scaleFactors[i_task] = taskScale;
      } else {
        scaleFactors[i_task] = 1.0;
        P = PS;
      }
    }
  }

//  double nSatTot = 0.0;
//  for (int i = 0; i < n_tasks; i++)
//    nSatTot += nSat[i];
//  return nSatTot;
  return 1.0;
}

//#define LOG_ACTIVE

double FOSNSVelocityIK::SNSsingle(int priority,
                                  const Eigen::VectorXd &higherPriorityJointVelocity,
                                  const Eigen::MatrixXd &higherPriorityNull,
                                  const Eigen::MatrixXd &jacobian,
                                  const Eigen::VectorXd &task,
                                  Eigen::VectorXd *jointVelocity,
                                  Eigen::MatrixXd *nullSpaceProjector)
{
  //INITIALIZATION
  Eigen::MatrixXd JPinverse;  //(J_k P_{k-1})^#
  Eigen::ArrayXd a, b;  // used to compute the task scaling factor
  bool limit_excedeed;
  double scalingFactor = 1.0;
  int mostCriticalJoint;
  bool singularTask = false;

  double base_Scale;
  double best_Scale = -1.0;
  Eigen::VectorXd best_dq1;
  Eigen::VectorXd best_dq2;
  Eigen::VectorXd best_dqw;
  //int best_nSat;

  Eigen::MatrixXd tildeZ;
  Eigen::VectorXd dq1, dq2;
  Eigen::VectorXd dqw;
  Eigen::VectorXd dq1_base, dq2_base;
  Eigen::VectorXd dqn = Eigen::VectorXd::Zero(n_dof);
  //Eigen::VectorXd dotQs;

  Eigen::MatrixXd zin;
  Eigen::VectorXd bin, bout;
  double dqw_in;

  //double mu_in;
  double mu_in1;
  double mu_inp2w;
  //double mu_out;
  double mu_out1;
  double mu_outp2w;
  double min_mu;
  int id_min_mu = n_dof + 1;  //just to be not possible
  Eigen::VectorXd scaledMU;

  bool computedScalingFactor = false;

#ifdef LOG_ACTIVE
  int n_in=0,n_out=0;
  std::stringstream log;
#endif

  //compute the base solution
  singularTask = !pinv_QR_Z(jacobian, higherPriorityNull, &JPinverse, &tildeZ);
  *nullSpaceProjector = tildeZ * tildeZ.transpose();
  dq1 = JPinverse * task;
  dq2 = -JPinverse * jacobian * higherPriorityJointVelocity;
  dqw = Eigen::VectorXd::Zero(n_dof);
  dq1_base = dq1;
  dq2_base = dq2;
  dotQ = higherPriorityJointVelocity + dq1 + dq2;
  a = dq1.array();
  b = dotQ.array() - a;
  getTaskScalingFactor(a, b, Eigen::VectorXi::Zero(n_dof), &scalingFactor, &mostCriticalJoint);

#ifdef LOG_ACTIVE
  log<<"task "<<priority<<std::endl;
  log<<"base Z norm "<<higherPriorityNull.norm()<<std::endl;
  log<<"base J*Z norm "<<(jacobian*higherPriorityNull).norm()<<std::endl;
  log<<"scale factor at 0 "<<scalingFactor<<std::endl;
  log<<"base S "<<S[priority].transpose()<<std::endl;
#endif

  if (scalingFactor >= 1.0) {
    // then is clearly the optimum since all joints velocity are computed with the pseudoinverse
    *jointVelocity = dotQ;
    //dotQopt[priority]=dotQ;
    nSat[priority] = 0;
    satList[priority].clear();
    S[priority] = Eigen::VectorXi::Zero(n_dof);
#ifdef LOG_ACTIVE
    log<<"task accomplished without saturations"<<std::endl;
    log<<"scale "<<scalingFactor<<std::endl;
    //log<<"last dotQ "<<higherPriorityJointVelocity->transpose()<<std::endl;
    //log<<"dotQ "<<dotQ.transpose()<<std::endl;
    //log<<"dotQmin "<<dotQmin.transpose()<<std::endl;
    //log<<"dotQmax "<<dotQmax.transpose()<<std::endl;
    if (singularTask){
      log<<"THE TASK IS SINGULAR"<<std::endl;
      ROS_WARN("\n%s\n\n",log.str().c_str());
    }
#endif
    return scalingFactor;
  } else {
    base_Scale = scalingFactor;
    if ((scalingFactor > best_Scale)) {
      //save best solution so far
      best_Scale = scalingFactor;
      best_dq1 = dq1;
      best_dq2 = dq2;
      best_dqw = dqw;
      //best_nSat = nSat[priority];
    }
  }

  if ((singularTask) || (base_Scale < 0)) {
    // the task is singular so return a scaled damped solution (no SNS possible)

    if (scalingFactor > 0.0) {
      *jointVelocity = higherPriorityJointVelocity + scalingFactor * dq1 + dq2;
      //dotQopt[priority]=(*jointVelocity);
      *nullSpaceProjector = tildeZ * tildeZ.transpose();
    } else {
      // the task is not executed
      *jointVelocity = higherPriorityJointVelocity;
      //dotQopt[priority]=(*jointVelocity);
      *nullSpaceProjector = higherPriorityNull;
    }
#ifdef LOG_ACTIVE
    if (singularTask) log<<"the task is singular"<<std::endl;
    if (base_Scale<0) log<<"base scale < 0"<<std::endl;
    log<<"scale "<<scalingFactor<<std::endl;
    //  log<<"dotQ "<<dotQ.transpose();
    //log<<"dotQmin "<<dotQmin.transpose();
    //  log<<"dotQmax "<<dotQmax.transpose();
    ROS_WARN("\n%s\n\n",log.str().c_str());
#endif
    nSat[priority] = 0;
    satList[priority].clear();
    S[priority] = Eigen::VectorXi::Zero(n_dof);
    return scalingFactor;
  }

//_______________________________________________________END of base computation

#define WARM_START
#ifdef  WARM_START
//############### WARM START
//B=Eigen::MatrixXd::Zero(n_dof,n_dof);
  if (nSat[priority]) {
//    log<<"started with "<<nSat[priority]<< " saturated\n";
    Eigen::MatrixXd Zws = Eigen::MatrixXd::Zero(nSat[priority], tildeZ.cols());
    Eigen::MatrixXd invZws = Eigen::MatrixXd::Zero(tildeZ.cols(), nSat[priority]);
    Eigen::MatrixXd forMu = Eigen::MatrixXd::Zero(nSat[priority], nSat[priority]);
    Eigen::MatrixXd Bws = Eigen::MatrixXd::Zero(n_dof, nSat[priority]);
    Eigen::VectorXd dq1_ws = Eigen::VectorXd::Zero(nSat[priority]);
    Eigen::VectorXd dq2_ws = Eigen::VectorXd::Zero(nSat[priority]);
    Eigen::VectorXd dqw_ws = Eigen::VectorXd::Zero(nSat[priority]);
    int idws = 0;
    bool invertibelZws;

    for (it = satList[priority].begin(); it != satList[priority].end();) {
      int id = *it;

      if (tildeZ.row(id).norm() < 1e-10) {
        //it means that the joint has been already saturated by the previous task
        if (it == satList[priority].begin()) {
          it = satList[priority].erase_after(satList[priority].before_begin());
        } else {
          it = satList[priority].erase_after(prev_it);
        }
        S[priority](id) = 0;
        nSat[priority]--;
      } else {

#ifdef LOG_ACTIVE
        n_in++;
        log<<id<< " ";
#endif

        Zws.row(idws) = tildeZ.row(id);
        dq1_ws(idws) = dq1_base(id);
        dq2_ws(idws) = dq2_base(id);
        if (S[priority](id) > 0.0) {
          dqw_ws(idws) = dotQmax(id) - higherPriorityJointVelocity(id);
        } else {
          dqw_ws(idws) = dotQmin(id) - higherPriorityJointVelocity(id);
        }
        dqn(id) = dqw_ws(idws);

        idws++;
        prev_it = it;
        it++;
      }

    }

    if (nSat[priority]) {
      Zws.conservativeResize(nSat[priority], tildeZ.cols());
      invertibelZws = pinv_QR(Zws, &invZws);

      if (!invertibelZws) {
        //  ROS_WARN("Zws is not invertible... what should I do?");
#ifdef LOG_ACTIVE_NO
        //log<<std::endl<<"Z \n"<<tildeZ;
        //log<<std::endl<<"Zws \n"<<Zws;
        log<<"\nS\n"<<S[priority].transpose()<<std::endl<<"norm Zws ";
        for (int i=0;i<nSat[priority];i++) {
          //for (it=satList[priority].begin(); it!=satList[priority].end(); ++it){
          //  int id=*it;
          //  log<<" "<<tildeZ.row(id).norm();
          log<<" "<<Zws.row(i).norm();
        }
        log<<"\n Zws rows "<<Zws.rows();
        ROS_WARN("%s",log.str().c_str());
        exit(1);
#endif
        satList[priority].clear();
        nSat[priority] = 0;
        S[priority] = Eigen::VectorXi::Zero(n_dof);

      } else {

        Bws = tildeZ * invZws;
        forMu = invZws.transpose() * invZws;
        //      forMu=Bws.transpose()*Bws;
        idws = 0;
        for (it = satList[priority].begin(); it != satList[priority].end(); ++it) {
          int id = *it;
          B.col(id) = Bws.col(idws);

          lagrangeMu1(id) = -forMu.row(idws) * dq1_ws;
          lagrangeMup2w(id) = forMu.row(idws) * (dqw_ws - dq2_ws);
          lagrangeMu(id) = lagrangeMu1(id) + lagrangeMup2w(id);

          idws++;
        }
        dq1 = dq1_base - Bws * dq1_ws;
        dq2 = dq2_base - Bws * dq2_ws;
        dqw = Bws * dqw_ws;
        dotQ = higherPriorityJointVelocity + dq1 + dq2 + dqw;

        tildeZ = tildeZ - Bws * Zws;
#ifdef LOG_ACTIVE
        //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        a=dq1.array();
        b=dotQ.array() - a;
        getTaskScalingFactor(a, b, S[priority], &scalingFactor, &mostCriticalJoint);
        scaledMU=scalingFactor*lagrangeMu1 + lagrangeMup2w;
        //find the minimum negative mu
        for (it=satList[priority].begin(); it!=satList[priority].end(); ++it) {
          int id=*it;
          if (dqn(id)>=0) scaledMU(id)=-scaledMU(id);
        }
        log<<"\nstart scale "<< scalingFactor<<std::endl;
        //log<<"\nstart scaled Mu "<< scaledMU.transpose()<<std::endl;
        log<<"\nstart scaled dq "<< (higherPriorityJointVelocity+scalingFactor*dq1+dq2+dqw).transpose()<<std::endl;
        log<<"\nstart S "<<S[priority].transpose();
        //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
#endif
        //find the minimum negative mu
        min_mu = 0;
        id_min_mu = n_dof + 1;
        for (it = satList[priority].begin(); it != satList[priority].end(); ++it) {
          int id = *it;
          if (dqn(id) >= 0)
            lagrangeMu(id) = -lagrangeMu(id);
          if ((lagrangeMu(id) < min_mu) && (abs(dqn(id)) > 1e-12)) {
            min_mu = lagrangeMu(id);
            id_min_mu = id;
          }
        }
#ifdef LOG_ACTIVE
        log<<"\nstart Mu "<< lagrangeMu.transpose()<<std::endl;
#endif
      }
    }
  }
  computedScalingFactor = false;
#else
//######################### else
  satList[priority].clear();
  nSat[priority]=0;
  S[priority]=Eigen::VectorXi::Zero(n_dof);

//  lagrangeMu1=Eigen::VectorXd::Zero(n_dof);
//  lagrangeMup2w=Eigen::VectorXd::Zero(n_dof);
//#########################
#endif

  //SNS
  int count = 0;
  do {
    count++;
    //ROS_INFO("%d",count);
    if (count > 2 * n_dof) {
#ifndef _ONLY_WARNING_ON_ERROR
      ROS_ERROR("Infinite loop on SNS for task (%d)", priority);

      exit(1);
#else
      ROS_WARN("Infinite loop on SNS for task (%d): nSat=%d ",priority,nSat[priority]);
      /*      //##############################
       string s;
       stringstream buffer;
       streambuf * old = std::cout.rdbuf(buffer.rdbuf());
       cout << "nSat\n" << nSat[priority]<<std::endl;
       s = buffer.str();
       ROS_INFO("\n p %d \n %s",priority,s.c_str());
       //#############################
       */
      // the task is not executed
      //nSat[priority]=0;
      satList[priority].clear();
      nSat[priority] = 0;
      S[priority] = Eigen::VectorXi::Zero(n_dof);
      *jointVelocity = higherPriorityJointVelocity;
      //dotQopt[priority]=(*jointVelocity);
      *nullSpaceProjector = higherPriorityNull;
      limit_excedeed=false;
#ifdef LOG_ACTIVE
      //log<<std::endl<<"last dq "<<dotQ.transpose()<< "\ndq1"<<dq1.transpose()<< "\ndq2"<<dq2.transpose()  << "\ndqn"<<dqn.transpose()<< "\ndqw"<<dqw.transpose()<< "\n\nS"<<S[priority].transpose()<<std::endl;
      log<<"last scaling factor "<< scalingFactor;
      ROS_WARN("%s",log.str().c_str());
      exit(1);
#endif
      //continue;
      return -1.0;
#endif
    }
    limit_excedeed = true;

    if (!computedScalingFactor) {
      a = dq1.array();
      b = dotQ.array() - a;
      getTaskScalingFactor(a, b, S[priority], &scalingFactor, &mostCriticalJoint);
    }
    computedScalingFactor = false;

#ifdef LOG_ACTIVE
    //if (n_in>2*n_dof-3){
//    log<<std::endl<<"last dq "<<dotQ.transpose()<< "\ndq1"<<dq1.transpose()<< "\ndq2"<<dq2.transpose()  << "\ndqn"<<dqn.transpose()<< "\ndqw"<<dqw.transpose()<< "\n\nS"<<S[priority].transpose()<<std::endl;
    log<<" last scaling factor "<< scalingFactor;
    log<<"\n n_in "<< n_in<<" n_out "<<n_out<<"  -> S "<<S[priority].transpose();
//      log<<"\n nSat"<< nSat[priority];

    log<<"\n best scale factor "<< best_Scale<<"\n\n";

    scaledMU=scalingFactor*lagrangeMu1 + lagrangeMup2w;
    //find the minimum negative mu
    for (it=satList[priority].begin(); it!=satList[priority].end(); ++it) {
      int id=*it;
      if (dqn(id)>=0) scaledMU(id)=-scaledMU(id);
    }
    //log<<"/nstop scaled Mu "<< scaledMU.transpose()<<std::endl;
    //log<<"stop scaled dq "<< (*higherPriorityJointVelocity)+scalingFactor*dq1+dq2+dqw<<std::endl;
    log<<"error "<<(jacobian*dotQ - task).norm();

    ROS_WARN("\n%s\n\n",log.str().c_str());
    log.str("");
    //    exit(1);
    //  }
#endif

    if ((scalingFactor >= 1.0) || (scalingFactor < 0.0)) {
      //check the optimality of the solution
      if (id_min_mu < n_dof) {
#ifdef LOG_ACTIVE
        n_out++;
        log<<" O"<<id_min_mu;
#endif
        //remove the id_min_mu joint from saturation and update B
        bout = B.col(id_min_mu);
        mu_out1 = lagrangeMu1(id_min_mu);
        mu_outp2w = lagrangeMup2w(id_min_mu);
        for (it = satList[priority].begin(); it != satList[priority].end();) {
          int id = *it;
          if (id == id_min_mu) {
            //remove id-min_mu from list of saturated joints
            if (it == satList[priority].begin()) {
              it = satList[priority].erase_after(satList[priority].before_begin());
            } else {
              it = satList[priority].erase_after(prev_it);
            }
          } else {
            //update B
            double baux = (double) bout.dot(B.col(id)) / bout.squaredNorm();
            B.col(id) -= bout * baux;
            prev_it = it;
            it++;
          }
        }
        //update the solution
        double dq1X = dq1_base(id_min_mu);
        double dq2X = dq2_base(id_min_mu);
        double dqwX = dqn(id_min_mu);
        Eigen::MatrixXd ZX = tildeZ.row(id_min_mu);
        for (it = satList[priority].begin(); it != satList[priority].end(); ++it) {
          int id = *it;
          lagrangeMu1(id) += bout(id) * mu_out1;
          lagrangeMup2w(id) += bout(id) * mu_outp2w;
          lagrangeMu(id) = lagrangeMu1(id) + lagrangeMup2w(id);
          dq1X -= B(id_min_mu, id) * dq1_base(id);
          dq2X -= B(id_min_mu, id) * dq2_base(id);
          dqwX -= B(id_min_mu, id) * dqn(id);
          ZX -= B(id_min_mu, id) * tildeZ.row(id);

        }
        dq1 += bout * dq1X;
        dq2 += bout * dq2X;
        dqw -= bout * dqwX;
        tildeZ += bout * ZX;
        nSat[priority]--;
        S[priority](id_min_mu) = 0;
        dotQ = higherPriorityJointVelocity + dq1 + dq2 + dqw;

        //find the minimum negative mu
        min_mu = 0;
        id_min_mu = n_dof + 1;
        for (it = satList[priority].begin(); it != satList[priority].end(); ++it) {
          int id = *it;
          if (dqn(id) >= 0)
            lagrangeMu(id) = -lagrangeMu(id);
          if ((lagrangeMu(id) < min_mu) && (abs(dqn(id)) > 1e-12)) {
            min_mu = lagrangeMu(id);
            id_min_mu = id;
          }
        }
        continue;
      }
    }
    if (scalingFactor >= 1.0) {
      // task accomplished
      *jointVelocity = dotQ;
      //dotQopt[priority]=(*jointVelocity);
      *nullSpaceProjector = tildeZ * tildeZ.transpose();  //if start net task from previous saturations
//        ROS_INFO("n_in %d n_out %d",n_in,n_out);

      return scalingFactor;

    } else if ((best_Scale > 0.0) || (scalingFactor > 0.0)) {
      if ((scalingFactor > best_Scale)) {
        //save best solution so far
        best_Scale = scalingFactor;
        best_dq1 = dq1;
        best_dq2 = dq2;
        best_dqw = dqw;
        //best_nSat = nSat[priority];
      }
      //check if the scaled solution is the optimum
      scaledMU = scalingFactor * lagrangeMu1 + lagrangeMup2w;
      //find the minimum negative mu
      min_mu = 0;
      id_min_mu = n_dof + 1;
      for (it = satList[priority].begin(); it != satList[priority].end(); ++it) {
        int id = *it;
        if (dqn(id) >= 0)
          scaledMU(id) = -scaledMU(id);
        //if ((scaledMU(id)<0) && (scaledMU(id)>min_mu)){
        if ((scaledMU(id) < min_mu) && (abs(dqn(id)) > 1e-12)) {
          min_mu = scaledMU(id);
          id_min_mu = id;
        }
      }

      if (id_min_mu < n_dof) {
#ifdef LOG_ACTIVE
        n_out++;
        log<<" Os"<<id_min_mu;
#endif
        //remove the id_min_mu joint from saturation and update B
        bout = B.col(id_min_mu);
        mu_out1 = lagrangeMu1(id_min_mu);
        mu_outp2w = lagrangeMup2w(id_min_mu);
        for (it = satList[priority].begin(); it != satList[priority].end();) {
          int id = *it;
          if (id == id_min_mu) {
            //remove id-min_mu from list of saturated joints
            if (it == satList[priority].begin()) {
              it = satList[priority].erase_after(satList[priority].before_begin());
            } else {
              it = satList[priority].erase_after(prev_it);
            }
          } else {
            //update B
            double baux = (double) bout.dot(B.col(id)) / bout.squaredNorm();
            B.col(id) -= bout * baux;
            prev_it = it;
            it++;
          }
        }
        //update the solution
        double dq1X = dq1_base(id_min_mu);
        double dq2X = dq2_base(id_min_mu);
        double dqwX = dqn(id_min_mu);
        Eigen::MatrixXd ZX = tildeZ.row(id_min_mu);
        for (it = satList[priority].begin(); it != satList[priority].end(); ++it) {
          int id = *it;
          lagrangeMu1(id) += bout(id) * mu_out1;
          lagrangeMup2w(id) += bout(id) * mu_outp2w;
          lagrangeMu(id) = lagrangeMu1(id) + lagrangeMup2w(id);
          dq1X -= B(id_min_mu, id) * dq1_base(id);
          dq2X -= B(id_min_mu, id) * dq2_base(id);
          dqwX -= B(id_min_mu, id) * dqn(id);
          ZX -= B(id_min_mu, id) * tildeZ.row(id);

        }
        dq1 += bout * dq1X;
        dq2 += bout * dq2X;
        dqw -= bout * dqwX;
        tildeZ += bout * ZX;
        nSat[priority]--;
        S[priority](id_min_mu) = 0;
        dotQ = higherPriorityJointVelocity + dq1 + dq2 + dqw;

        //find the minimum negative mu
        min_mu = 0;
        id_min_mu = n_dof + 1;
        for (it = satList[priority].begin(); it != satList[priority].end(); ++it) {
          int id = *it;
          if (dqn(id) >= 0)
            lagrangeMu(id) = -lagrangeMu(id);
          if ((lagrangeMu(id) < min_mu) && (abs(dqn(id)) > 1e-12)) {
            min_mu = lagrangeMu(id);
            id_min_mu = id;
          }
        }
        continue;
      }

    }

    int idxW = mostCriticalJoint;

    //saturate the most critical joint
    zin = tildeZ.row(idxW);
    //if we do not use norm(zin) then this part can go first
    bin = tildeZ * (zin.transpose() / zin.squaredNorm());

    dq1 -= bin * dq1(idxW);
    dq2 -= bin * dq2(idxW);

    if (dotQ(idxW) >= 0.0) {
      dqw_in = dotQmax(idxW) - higherPriorityJointVelocity(idxW);
      S[priority](idxW) = +1;
    } else {
      dqw_in = dotQmin(idxW) - higherPriorityJointVelocity(idxW);
      S[priority](idxW) = -1;
    }
    dqw += bin * (dqw_in - dqw(idxW));
    dqn(idxW) = dqw_in;

    dotQ = higherPriorityJointVelocity + dq1 + dq2 + dqw;

    a = dq1.array();
    b = dotQ.array() - a;
    getTaskScalingFactor(a, b, S[priority], &scalingFactor, &mostCriticalJoint);
    computedScalingFactor = true;

#ifdef LOG_ACTIVE
    log<<" I "<<idxW<< "norm zin "<<zin.norm()<<" at "<<dqw_in<<" with scale "<<scalingFactor<<std::endl;
    if (scalingFactor<0) log<<"IT WAS < 0 "<<std::endl;
    if (scalingFactor<1e-12) log<<"IT WAS < eps "<<std::endl;
    n_in++;
    //  ROS_WARN("\n%s\n\n",log.str().c_str());
    //  log.str("");

#endif
    if ((zin.norm() < 1e-8) || (scalingFactor < 1e-12)) {
      S[priority](idxW) = 0;  //put it back to the right value
      if (best_Scale >= 0) {
        //take the best solution
        *jointVelocity = higherPriorityJointVelocity + best_Scale * best_dq1 + best_dq2 + best_dqw;
        //dotQopt[priority]=(*jointVelocity);
        //nSat[priority]=best_nSat;
        *nullSpaceProjector = tildeZ * tildeZ.transpose();  //if start net task from previous saturations
//        ROS_INFO("n_in %d n_out %d",n_in,n_out);
        if (best_Scale == base_Scale) {
          //no saturation was needed to obtain the best scale
          satList[priority].clear();
          nSat[priority] = 0;
          S[priority] = Eigen::VectorXi::Zero(n_dof);

        }
        //if (priority==1) *jointVelocity=(*higherPriorityJointVelocity);
        return best_Scale;
      } else {
        //no solution
        //nSat[priority]=0;
        *jointVelocity = higherPriorityJointVelocity;
        //dotQopt[priority]=(*jointVelocity);
        *nullSpaceProjector = higherPriorityNull;
        limit_excedeed = false;
        //continue;
#ifdef LOG_ACTIVE_UNFEASIBLE
        log<<std::endl<<"last dq "<<dotQ.transpose()<< "\ndq1"<<dq1.transpose()<< "\ndq2"<<dq2.transpose() << "\ndqn"<<dqn.transpose()<< "\ndqw"<<dqw.transpose()<< "\n\nS"<<S[priority].transpose()<<std::endl;
        log<<"last scaling factor "<< scalingFactor;
        ROS_WARN("%s",log.str().c_str());
        exit(1);
#endif
        return -1.0;
      }
    }

    nSat[priority]++;
    tildeZ -= bin * zin;

    //update mu and B
    min_mu = 0;
    id_min_mu = n_dof + 1;
    mu_in1 = bin.dot(dq1);
    mu_inp2w = bin.dot(dotQ - dq1);
    for (it = satList[priority].begin(); it != satList[priority].end(); ++it) {
      int id = *it;
      lagrangeMu1(id) -= B(idxW, id) * mu_in1;
      lagrangeMup2w(id) -= B(idxW, id) * mu_inp2w;
      lagrangeMu(id) = lagrangeMu1(id) + lagrangeMup2w(id);
      B.col(id) -= bin * B(idxW, id);
      if (dqn(id) >= 0)
        lagrangeMu(id) = -lagrangeMu(id);
      if ((lagrangeMu(id) < min_mu) && (abs(dqn(id)) > 1e-12)) {
        min_mu = lagrangeMu(id);
        id_min_mu = id;
      }
    }
    satList[priority].push_front(idxW);
    lagrangeMu1(idxW) = mu_in1;
    lagrangeMup2w(idxW) = mu_inp2w;
    //lagrangeMu(idxW)=lagrangeMu1(idxW)+lagrangeMup2w(idxW);

  } while (limit_excedeed);  //actually in this implementation if we use while(1) it would be the same

  *nullSpaceProjector = tildeZ * tildeZ.transpose();
  (*jointVelocity) = dotQ;
  //dotQopt[priority]=(*jointVelocity);
  return scalingFactor;
}

}  // namespace sns_ik
