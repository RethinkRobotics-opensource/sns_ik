/*! \file fsns_velocity_ik.cpp
 * \brief Fast SNS velocity IK solver
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

#include <sns_ik/fsns_velocity_ik.hpp>

#include "sns_ik_math_utils.hpp"

namespace sns_ik {

FSNSVelocityIK::FSNSVelocityIK(int dof, double loop_period) :
  SNSVelocityIK(dof, loop_period)
{
}


double FSNSVelocityIK::getJointVelocity(Eigen::VectorXd *jointVelocity,
    const std::vector<Task> &sot,
    const Eigen::VectorXd &jointConfiguration)
{
  // This will only reset member variables if different from previous values
  setNumberOfTasks(sot.size(), sot[0].jacobian.cols());
  S.resize(n_tasks, Eigen::VectorXi::Zero(n_dof));

  // TODO: check that setJointsCapabilities has been already called

  //P_0=I
  //dq_0=0
  Eigen::MatrixXd P = Eigen::MatrixXd::Identity(n_dof, n_dof);
  *jointVelocity = Eigen::VectorXd::Zero(n_dof, 1);
  Eigen::VectorXd higherPriorityJointVelocity;
  Eigen::MatrixXd higherPriorityNull;

  shapeJointVelocityBound(jointConfiguration);

  for (int i_task = 0; i_task < n_tasks; i_task++) {  //consider all tasks
    higherPriorityJointVelocity = *jointVelocity;
    higherPriorityNull = P;
    scaleFactors[i_task] = SNSsingle(i_task, higherPriorityJointVelocity, higherPriorityNull,
        sot[i_task].jacobian, sot[i_task].desired, jointVelocity, &P);

    if (scaleFactors[i_task] > 1)
          scaleFactors[i_task] = 1;
  }

  // TODO: what is being returned here?
  //return 1.0 * nSat[0] + nSat[1];
  return 1.0;
}

double FSNSVelocityIK::SNSsingle(int priority,
                                 const Eigen::VectorXd &higherPriorityJointVelocity,
                                 const Eigen::MatrixXd &higherPriorityNull,
                                 const Eigen::MatrixXd &jacobian,
                                 const Eigen::VectorXd &task,
                                 Eigen::VectorXd *jointVelocity,
                                 Eigen::MatrixXd *nullSpaceProjector)
{
  //FIXME: THERE IS A PROBLEM if we use 3 tasks... to be checked

  //INITIALIZATION
  Eigen::MatrixXd JPinverse;  //(J_k P_{k-1})^#
  Eigen::ArrayXd a, b;  // used to compute the task scaling factor
  bool limit_excedeed;
  double scalingFactor = 1.0;
  int mostCriticalJoint;
  bool singularTask = false;

  double best_Scale = -1.0;
  Eigen::VectorXd best_dq1;
  Eigen::VectorXd best_dq2;
  Eigen::VectorXd best_dqw;
  //int best_nSat;

  Eigen::MatrixXd tildeZ;
  Eigen::VectorXd dq1, dq2, dqw;

  Eigen::MatrixXd bin, zin;
  double dqw_in;

  //initialization
  nSat[priority] = 0;
  S[priority] = Eigen::VectorXi::Zero(n_dof);

  //compute the base solution
  singularTask = !pinv_QR_Z(jacobian, higherPriorityNull, &JPinverse, &tildeZ);
  *nullSpaceProjector = tildeZ * tildeZ.transpose();
  dq1 = JPinverse * task;
  dq2 = -JPinverse * jacobian * higherPriorityJointVelocity;
  dqw = Eigen::VectorXd::Zero(n_dof);

  dotQ = higherPriorityJointVelocity + dq1 + dq2;
  a = dq1.array();
  b = dotQ.array() - a;
  getTaskScalingFactor(a, b, S[priority], &scalingFactor, &mostCriticalJoint);
  //getTaskScalingFactor(a, b, Eigen::MatrixXd::Identity(n_dof, n_dof), &scalingFactor, &mostCriticalJoint);

  if (scalingFactor >= 1.0) {
    // then is clearly the optimum since all joints velocity are computed with the pseudoinverse
    (*jointVelocity) = dotQ;
    nSat[priority] = 0;
    dotQopt[priority] = dotQ;
    return 1.0;
  }

  if (singularTask) {
    // the task is singular so return a scaled damped solution (no SNS possible)
    if (scalingFactor > 0.0) {
      nSat[priority] = 0;
      (*jointVelocity) = higherPriorityJointVelocity + scalingFactor * dq1 + dq2;
      dotQopt[priority] = (*jointVelocity);
    } else {
      // the task is not executed
      nSat[priority] = 0;
      *jointVelocity = higherPriorityJointVelocity;
      dotQopt[priority] = (*jointVelocity);
      *nullSpaceProjector = higherPriorityNull;
    }
    return scalingFactor;
  }

  if (scalingFactor > best_Scale) {
    //save best solution so far
    best_Scale = scalingFactor;
    best_dq1 = dq1;
    best_dq2 = dq2;
    best_dqw = dqw;
    //best_nSat = 0;

  }

//_______________________________________________________END of base computation

  //SNS
  int count = 0;
  do {
    count++;
    if (count > 2 * n_dof) {

      // the task is not executed
      nSat[priority]=0;
      *jointVelocity = higherPriorityJointVelocity;
      dotQopt[priority] = *jointVelocity;
      *nullSpaceProjector = higherPriorityNull;
      limit_excedeed=false;
      //continue;
      return -1.0;
    }
    limit_excedeed = true;

    //saturate the most critical joint
    zin = tildeZ.row(mostCriticalJoint);

    if ((zin.norm() < 1e-8) || (scalingFactor < 1e-6)) {
      if (best_Scale >= 0) {
        //take the best solution
        *jointVelocity = higherPriorityJointVelocity + best_Scale * best_dq1 + best_dq2 + best_dqw;
        dotQopt[priority] = (*jointVelocity);
        //nSat[priority]=best_nSat;
        *nullSpaceProjector = tildeZ * tildeZ.transpose();  //if start net task from previous saturations
        return best_Scale;
      } else {
        //no solution
        //nSat[priority]=0;
        *jointVelocity = higherPriorityJointVelocity;
        dotQopt[priority] = (*jointVelocity);
        *nullSpaceProjector = higherPriorityNull;
        limit_excedeed = false;
        //continue;
        return -1.0;
      }
    }

    //if we do not use norm(zin) then this part can go first
    bin = tildeZ * (zin.transpose() / zin.squaredNorm());

    dq1 -= bin * dq1(mostCriticalJoint);
    dq2 -= bin * dq2(mostCriticalJoint);

    if (dotQ(mostCriticalJoint) >= 0.0) {
      dqw_in = dotQmax(mostCriticalJoint) - higherPriorityJointVelocity(mostCriticalJoint);
    } else {
      dqw_in = dotQmin(mostCriticalJoint) - higherPriorityJointVelocity(mostCriticalJoint);
    }
    dqw += bin * (dqw_in - dqw(mostCriticalJoint));

    dotQ = higherPriorityJointVelocity + dq1 + dq2 + dqw;

    nSat[priority]++;
    S[priority](mostCriticalJoint) = nSat[priority];
    tildeZ -= bin * zin;

    a = dq1.array();
    b = dotQ.array() - a;
    getTaskScalingFactor(a, b, S[priority], &scalingFactor, &mostCriticalJoint);
    if (scalingFactor >= 1.0) {
      // task accomplished
      *jointVelocity = dotQ;
      dotQopt[priority] = (*jointVelocity);
      *nullSpaceProjector = tildeZ * tildeZ.transpose();  //if start net task from previous saturations
      return 1.0;
    } else {
      if ((scalingFactor > best_Scale)) {
        //save best solution so far
        best_Scale = scalingFactor;
        best_dq1 = dq1;
        best_dq2 = dq2;
        best_dqw = dqw;
        //best_nSat = nSat[priority];
      }
    }

  } while (limit_excedeed);  //actually in this implementation if we use while(1) it would be the same

  *nullSpaceProjector = tildeZ * tildeZ.transpose();
  *jointVelocity = dotQ;
  return 1.0;
}

void FSNSVelocityIK::getTaskScalingFactor(const Eigen::ArrayXd &a,
                  const Eigen::ArrayXd &b,
                  const Eigen::VectorXi &S, double *scalingFactor,
                  int *mostCriticalJoint)
{
  Eigen::ArrayXd Smin, Smax;
  double temp, smax, smin;
  double inf = INF;
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
    if (S(i) || a(i) == 0) {  // if it is not 0 (safer)
      Smin(i) = -inf;
      Smax(i) = inf;

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

}  // namespace sns_ik
