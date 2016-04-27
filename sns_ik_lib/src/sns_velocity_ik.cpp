/*! \file sns_velocity_ik.cpp
 * \brief Basic SNS velocity IK solver
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

#include <sns_ik/sns_velocity_ik.hpp>

#include <ros/ros.h>
#include <iostream>

using namespace Eigen;
using namespace sns_ik;

SNSVelocityIK::SNSVelocityIK(int dof, Scalar loop_period) :
  n_dof(0),
  n_tasks(0),
  m_usePositionLimits(true)
{
  setNumberOfDOF(dof);
  setLoopPeriod(loop_period);
}

bool SNSVelocityIK::setJointsCapabilities(VectorD limit_low, VectorD limit_high,
                                          VectorD maxVelocity, VectorD maxAcceleration)
{
  if (limit_low.rows() != n_dof || limit_high.rows() != n_dof
      || maxVelocity.rows() != n_dof || maxAcceleration.rows() != n_dof) {
    return false;
  }

  jointLimit_low = limit_low;
  jointLimit_high = limit_high;
  maxJointVelocity = maxVelocity;
  maxJointAcceleration = maxAcceleration;

  dotQmin = -maxJointVelocity.array();
  dotQmax = maxJointVelocity.array();

  return true;
}

void SNSVelocityIK::setNumberOfDOF(int dof)
{
  if (dof > 0 && dof != n_dof) {
    n_dof = dof;
    I = MatrixD::Identity(n_dof, n_dof);
    dotQ = VectorD::Zero(n_dof);
  }
}

void SNSVelocityIK::setNumberOfTasks(int ntasks, int dof)
{
  setNumberOfDOF(dof);

  if (n_tasks != ntasks) {
    n_tasks = ntasks;
    Scalar scale = 1.0;
    VectorD dq = VectorD::Zero(n_dof);

    W.resize(n_tasks, I);
    scaleFactors.resize(n_tasks, scale);
    dotQopt.resize(n_tasks, dq);
    nSat.resize(n_tasks, 0);
  }
}

Scalar SNSVelocityIK::getJointVelocity_STD(VectorD *jointVelocity,
                                           const std::vector<Task> &sot)
{
  int n_task = sot.size();
  int robotDOF = sot[0].jacobian.cols();

  //P_0=I
  //dq_0=0
  MatrixD P = MatrixD::Identity(robotDOF, robotDOF);
  *jointVelocity = VectorD::Zero(robotDOF, 1);
  MatrixD tmp, invJ;

  for (int i_task = 0; i_task < n_task; i_task++) {  //consider all tasks
    // dq_k = dq_{k-1} - (J_k P_{k-1})^# (dx_k - J_k dq_{k-1})
    // P_k = P_{k-1} - (J_k P_{k-1})^# J_k P_{k-1}
    tmp = sot[i_task].jacobian * P;

    pinv_damped_P(tmp, &invJ, &P);

    *jointVelocity = ((*jointVelocity) + invJ * (sot[i_task].desired - sot[i_task].jacobian * (*jointVelocity)));
  }

  return 1.0;
}

Scalar SNSVelocityIK::getJointVelocity(VectorD *jointVelocity, const std::vector<Task> &sot,
                                       const VectorD &jointConfiguration)
{
  // This will only reset member variables if different from previous values
  setNumberOfTasks(sot.size(), sot[0].jacobian.cols());

  // TODO: check that setJointsCapabilities has been already called

  //P_0=I
  //dq_0=0
  MatrixD P = MatrixD::Identity(n_dof, n_dof);
  *jointVelocity = VectorD::Zero(n_dof, 1);
  VectorD higherPriorityJointVelocity;
  MatrixD higherPriorityNull;

  shapeJointVelocityBound(jointConfiguration);

  for (int i_task = 0; i_task < n_tasks; i_task++) {  //consider all tasks
    higherPriorityJointVelocity = *jointVelocity;
    higherPriorityNull = P;
    scaleFactors[i_task] = SNSsingle(i_task, higherPriorityJointVelocity, higherPriorityNull,
        sot[i_task].jacobian, sot[i_task].desired, jointVelocity, &P);
  }

  //return 1.0;
  return scaleFactors[0];
}

void SNSVelocityIK::shapeJointVelocityBound(const VectorD &actualJointConfiguration, double margin) {

  // it could be written using the Eigen::Array potentiality
  double step, max, stop;

  for (int i = 0; i < n_dof; i++) {
    // for the minimum bound
    max = -maxJointVelocity(i);
    if (m_usePositionLimits) {
      step = (jointLimit_low(i) - actualJointConfiguration(i)) / loop_period;
      stop = -std::sqrt(2 * maxJointAcceleration(i) * (actualJointConfiguration(i) - jointLimit_low(i)));
      // take the maximum
      dotQmin(i) = std::max({step, max, stop});
    } else {
      dotQmin(i) = max;
    }

    // for the maximum bound
    max = maxJointVelocity(i);
    if (m_usePositionLimits) {
      step = (jointLimit_high(i) - actualJointConfiguration(i)) / loop_period;
      stop = std::sqrt(2 * maxJointAcceleration(i) * (jointLimit_high(i) - actualJointConfiguration(i)));
      // take the minimum
      dotQmax(i) = std::min({step, max, stop});
    } else {
      dotQmax(i) = max;
    }
  }

  dotQmin *= margin;
  dotQmax *= margin;
}

Scalar SNSVelocityIK::SNSsingle(int priority,
                                const VectorD &higherPriorityJointVelocity,
                                const MatrixD &higherPriorityNull,
                                const MatrixD &jacobian,
                                const VectorD &task,
                                VectorD *jointVelocity,
                                MatrixD *nullSpaceProjector)
{
  //INITIALIZATION
  VectorD tildeDotQ;
  MatrixD projectorSaturated;  //(((I-W_k)*P_{k-1})^#
  MatrixD JPinverse;  //(J_k P_{k-1})^#
  MatrixD temp;
  bool isW_identity;
  MatrixD barP = higherPriorityNull;
  Array<Scalar, Dynamic, 1> a, b;  // used to compute the task scaling factor
  bool limit_excedeed;
  bool singularTask = false;
  bool reachedSingularity = false;
  Scalar scalingFactor = 1.0;
  int mostCriticalJoint;
  //best solution
  Scalar bestScale = -1.0;
  VectorD bestTildeDotQ;
  MatrixD bestInvJP;
  VectorD bestDotQn;
  VectorD dotQn;  //saturate velocity in the null space

  //INIT
  W[priority] = I;
  isW_identity = true;
  dotQn = VectorD::Zero(n_dof);

  //SNS
  int count = 0;
  do {
    count++;
    //ROS_INFO("%d",count);
    if (count > 2 * n_dof) {
#ifndef _ONLY_WARNING_ON_ERROR
      //ROS_ERROR("Infinite loop on SNS for task (%d)", priority);

      //ROS_INFO("p:%d  scale:%f  mc:%d  sing:%d", priority, scalingFactor, mostCriticalJoint, (int) reachedSingularity);
      //##############################
      string s;
      stringstream buffer;
      streambuf * old = std::cout.rdbuf(buffer.rdbuf());
      cout << "W\n" << W[priority] << endl << "P(k-1)\n" << (*higherPriorityNull) << endl << "Psat\n"
          << projectorSaturated << endl << "dotQ\n" << dotQ << endl << "JPinverse\n" << JPinverse << endl;
      s = buffer.str();
      //ROS_INFO("\n p %d \n %s", priority, s.c_str());
      //#############################

      exit(1);
#else
      ROS_WARN("Infinite loop on SNS for task (%d)", priority);
      ROS_INFO("p:%d  scale:%f  mc:%d  sing:%d", priority, scalingFactor, mostCriticalJoint, (int)reachedSingularity);
      // the task is not executed
      *jointVelocity = higherPriorityJointVelocity;
      *nullSpaceProjector = higherPriorityNull;
      limit_excedeed = false;
      continue;
#endif
    }
    limit_excedeed = false;

    // remember that in the SNS W==I always and only on the first loop
    if (isW_identity) {
      tildeDotQ = higherPriorityJointVelocity;
      //compute (J P)^#
      temp = jacobian * higherPriorityNull;
      singularTask = !pinv_damped_P(temp, &JPinverse, nullSpaceProjector);
    } else {
      //JPinverse is already computed
      tildeDotQ = higherPriorityJointVelocity + projectorSaturated * dotQn;
    }
    dotQ = tildeDotQ + JPinverse * (task - jacobian * tildeDotQ);

    a = (JPinverse * task).array();
    b = dotQ.array() - a;

    getTaskScalingFactor(a, b, W[priority], &scalingFactor, &mostCriticalJoint);

    if (scalingFactor >= 1.0) {
      // task accomplished
    } else {
      limit_excedeed = true;
      if (singularTask) {
        // the task is singular so return a scaled damped solution (no SNS possible)
        //ROS_WARN("task %d is singular, scaling factor: %f",priority,scalingFactor);
        if (scalingFactor >= 0.0) {
          (*jointVelocity) = tildeDotQ + JPinverse * (scalingFactor * task - jacobian * tildeDotQ);
        } else {
          // the task is not executed
          //ROS_INFO("task not executed: J sing");
          //W[priority]=I;
          //dotQn=VectorD::Zero(n_dof);
          *jointVelocity = higherPriorityJointVelocity;
          *nullSpaceProjector = higherPriorityNull;
        }

        return scalingFactor;
      }

      if ((scalingFactor > bestScale)) {
        // save best solution so far
        bestScale = scalingFactor;
        bestTildeDotQ = tildeDotQ;
        bestInvJP = JPinverse;
        bestDotQn = dotQn;
      }

      // saturate the most critical join
      W[priority](mostCriticalJoint, mostCriticalJoint) = 0.0;
      isW_identity = false;
      if (dotQ(mostCriticalJoint) > dotQmax(mostCriticalJoint)) {
        dotQn(mostCriticalJoint) = dotQmax(mostCriticalJoint) - higherPriorityJointVelocity(mostCriticalJoint);
      } else {
        dotQn(mostCriticalJoint) = dotQmin(mostCriticalJoint) - higherPriorityJointVelocity(mostCriticalJoint);
      }

      if (priority == 0) {  //for the primary task higherPriorityNull==I
        barP = W[0];
        projectorSaturated = (I - W[0]);
      } else {
        temp = (I - W[priority]);
        reachedSingularity |= !pinv_forBarP(temp, higherPriorityNull, &projectorSaturated);

        barP = (I - projectorSaturated) * higherPriorityNull;
      }

      temp = jacobian * barP;

      reachedSingularity |= !pinv(temp, &JPinverse);

      if (reachedSingularity) {
        if (bestScale >= 0.0) {
          //ROS_INFO("best solution %f",bestScale);
          dotQn = bestDotQn;
          dotQ = bestTildeDotQ + bestInvJP * (bestScale * task - jacobian * bestTildeDotQ);
          //use the best solution found... no further saturation possible
          (*jointVelocity) = dotQ;
        } else {
          // the task is not executed
          ROS_WARN("task not executed: reached sing");
          *jointVelocity = higherPriorityJointVelocity;
          *nullSpaceProjector = higherPriorityNull;
        }

        return bestScale;
      }

    }

  } while (limit_excedeed);

  (*jointVelocity) = dotQ;
  return 1.0;
}

void SNSVelocityIK::getTaskScalingFactor(const Array<Scalar, Dynamic, 1> &a,
                                         const Array<Scalar, Dynamic, 1> &b,
                                         const MatrixD &W, Scalar *scalingFactor,
                                         int *mostCriticalJoint)
{
  Array<Scalar, Dynamic, 1> Smin, Smax;
  Scalar temp, smax, smin;
  Scalar inf = INF;
  int col;

  Smin = (dotQmin - b) / a;
  Smax = (dotQmax - b) / a;

  for (int i = 0; i < a.rows(); i++) {
    //switch
    if (Smin(i) > Smax(i)) {
      temp = Smin(i);
      Smin(i) = Smax(i);
      Smax(i) = temp;
    }
    //remove saturated
    if ((W(i, i) < 0.2) || (a(i) == 0)) {  // if it is not 1 (safer)
      Smin(i) = -INF;
      Smax(i) = INF;
    }
  }

  smax = Smax.minCoeff(mostCriticalJoint, &col);
  smin = Smin.maxCoeff();

  if ((smin > smax) || (smax < 0.0) || (smin > 1.0) || (smax == inf)) {
    (*scalingFactor) = -1.0;  // the task is not feasible
  } else {
    (*scalingFactor) = smax;
  }

}
