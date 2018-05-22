/*! \file osns_velocity_ik.cpp
 * \brief Optimal SNS velocity IK solver
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

#include <sns_ik/osns_velocity_ik.hpp>

#include "sns_ik_math_utils.hpp"

namespace sns_ik {

OSNSVelocityIK::OSNSVelocityIK(int dof, double loop_period) :
  SNSVelocityIK(dof, loop_period)
{
}


double OSNSVelocityIK::getJointVelocity(Eigen::VectorXd *jointVelocity,
    const std::vector<Task> &sot,
    const Eigen::VectorXd &jointConfiguration)
{
  // This will only reset member variables if different from previous values
  setNumberOfTasks(sot.size(), sot[0].jacobian.cols());

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
  }

  // TODO: verify what is being set here
  //return 1.0 * nSat[0] + nSat[1] + nSat[2] + nSat[3] + nSat[4];
  return 1.0;
}


double OSNSVelocityIK::SNSsingle(int priority,
                                const Eigen::VectorXd &higherPriorityJointVelocity,
                                const Eigen::MatrixXd &higherPriorityNull,
                                const Eigen::MatrixXd &jacobian,
                                const Eigen::VectorXd &task,
                                Eigen::VectorXd *jointVelocity,
                                Eigen::MatrixXd *nullSpaceProjector)
{
  //INITIALIZATION
  //Eigen::VectorXd tildeDotQ;
  Eigen::MatrixXd projectorSaturated;  //(((I-W_k)*P_{k-1})^#
  Eigen::MatrixXd JPinverse;  //(J_k P_{k-1})^#
  Eigen::MatrixXd temp;
  bool isW_identity;
  Eigen::MatrixXd barP = higherPriorityNull;  // remove the pointer arguments advantage... but I need to modify it
  Eigen::ArrayXd a, b;  // used to compute the task scaling factor
  bool limit_excedeed;
  bool singularTask = false;
  bool reachedSingularity = false;
  double scalingFactor = 1.0;
  int mostCriticalJoint;

  //best solution
  double bestScale = -0.1;
  Eigen::MatrixXd bestW;  //(only in OSNS)
  Eigen::MatrixXd bestInvJP;
  Eigen::VectorXd bestDotQn;
  Eigen::MatrixXd bestTildeP;

  Eigen::VectorXd dotQn;  //saturate velocity in the null space
  Eigen::VectorXd dotQs;
  Eigen::MatrixXd tildeP;  // used in the  OSNS

  //these two are needed to consider also W=I in case of a non feasible task
  bool invJPcomputed = false;

  //nSat[priority] = 0;

  //Compute the solution with W=I it is needed anyway to obtain nullSpaceProjector
  //compute (J P)^#
  temp = jacobian * higherPriorityNull;
  singularTask = !pinv_damped_P(temp, &JPinverse, nullSpaceProjector);

  tildeP = Eigen::MatrixXd::Zero(n_dof, n_dof);
  dotQs = higherPriorityJointVelocity + JPinverse * (task - jacobian * higherPriorityJointVelocity);
  a = (JPinverse * task).array();
  b = dotQs.array() - a;
  getTaskScalingFactor(a, b, I, &scalingFactor, &mostCriticalJoint);

  //double scalingI=scalingFactor;
  if (scalingFactor >= 1.0) {
    // this is clearly the optimum since all joints velocity are computed with the pseudoinverse
    (*jointVelocity) = dotQs;
    W[priority] = I;
    dotQopt[priority] = dotQs;
    return scalingFactor;
  }

  if (singularTask) {
    // the task is singular so return a scaled damped solution (no SNS possible)
    if (scalingFactor >= 0.0) {
      W[priority] = I;
      (*jointVelocity) = higherPriorityJointVelocity + scalingFactor * JPinverse * task
          + tildeP * higherPriorityJointVelocity;
      dotQopt[priority] = *jointVelocity;
    } else {
      // the task is not executed
      W[priority] = I;
      *jointVelocity = higherPriorityJointVelocity;
      dotQopt[priority] = *jointVelocity;
      *nullSpaceProjector = higherPriorityNull;
    }
    return scalingFactor;
  }

  if (scalingFactor > bestScale) {
    //save best solution so far
    bestScale = scalingFactor;
    //bestTildeDotQ=tildeDotQ;
    bestInvJP = JPinverse;
    bestW = I;
    bestTildeP = tildeP;
    //bestPS=projectorSaturated;
    bestDotQn = Eigen::VectorXd::Zero(n_dof);
  }

  //W[priority] = I;  //test: do not use the warm start
//----------------------------------------------------------------------- END W=I

  //INIT
  dotQn = Eigen::VectorXd::Zero(n_dof);
  if (isIdentity (W[priority])) {
    isW_identity = true;
    dotQopt[priority] = dotQs;  // use the one computed above
  } else {
    isW_identity = false;
    for (int i = 0; i < n_dof; i++) {
      if ((W[priority])(i, i) < 0.1) {  // equal to 0.0 (safer)
        //nSat[priority]++;
        //I'm not considering a different dotQ for each task... this could be a little improvement
        if (dotQopt[priority](i) >= 0.0) {
          dotQn(i) = dotQmax(i) - higherPriorityJointVelocity(i);
        } else {
          dotQn(i) = dotQmin(i) - higherPriorityJointVelocity(i);
        }
      } else
        dotQn(i) = 0.0;
    }
  }

  //SNS
  int count = 0;
  do {
    reachedSingularity = false;
    count++;
    if (count > 2 * n_dof) {
      // the task is not executed
      if (bestScale >= 0.0) {
        W[priority]=bestW;
        dotQopt[priority] = higherPriorityJointVelocity
            + bestInvJP *( bestScale * task - jacobian * higherPriorityJointVelocity)
            + bestTildeP * bestDotQn;
      } else {
        W[priority] = I;
        dotQopt[priority] = higherPriorityJointVelocity;
      }

      *jointVelocity = dotQopt[priority];
      return bestScale;
    }
    limit_excedeed = false;
    // If W=I everything is done --> go to the saturation phase

    if (!isW_identity && !invJPcomputed) {
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

      invJPcomputed = false;  //it needs to be computed for the next step
    }
    if (!isW_identity) {

      tildeP = (I - JPinverse * jacobian) * projectorSaturated;

      //compute the joint velocity
      dotQopt[priority] = higherPriorityJointVelocity
          + JPinverse * (task - jacobian * higherPriorityJointVelocity)
          + tildeP * dotQn;

      //compute the scaling factor
      a = (JPinverse * task).array();
      b = dotQopt[priority].array() - a;
      getTaskScalingFactor(a, b, W[priority], &scalingFactor, &mostCriticalJoint);

      if ((scalingFactor >= 1.0) || (scalingFactor < 0)) {

        //check optimality
        if (!isOptimal(priority, dotQopt[priority], tildeP, &W[priority], &dotQn)) {
          //modified W and dotQn
          limit_excedeed = true;
          //ROS_INFO("non OPT");
          continue;
        }
      }
      if (scalingFactor >= 1.0) {  //solution found

        //here limit_excedeed=false
        continue;
      }

      //if no solution found
      limit_excedeed = true;

      // is scaled an optimum
      if (scalingFactor >= 0) {
        dotQs = higherPriorityJointVelocity
            + JPinverse * (scalingFactor * task - jacobian * higherPriorityJointVelocity) + tildeP * dotQn;
        if (!isOptimal(priority, dotQs, tildeP, &W[priority], &dotQn)) {
          //ROS_INFO("non OPT s");
          //modified W and dotQn
          limit_excedeed = true;
          continue;
        }
      }

      if (scalingFactor > bestScale) {

        //save best solution so far
        bestScale = scalingFactor;
        bestInvJP = JPinverse;
        bestW = W[priority];
        bestTildeP = tildeP;
        bestDotQn = dotQn;
      }

    } else
      limit_excedeed = true;  //if it was W=I, to be here the limit is exceeded

    //nSat[priority]++;
    W[priority](mostCriticalJoint, mostCriticalJoint) = 0.0;
    isW_identity = false;
    if (dotQopt[priority](mostCriticalJoint) > dotQmax(mostCriticalJoint)) {
      dotQn(mostCriticalJoint) = dotQmax(mostCriticalJoint) - higherPriorityJointVelocity(mostCriticalJoint);
    } else {
      dotQn(mostCriticalJoint) = dotQmin(mostCriticalJoint) - higherPriorityJointVelocity(mostCriticalJoint);
    }

    //compute JPinverse
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

    invJPcomputed = true;
    //if reachedSingularity then take the best solution

    if ((reachedSingularity) || (scalingFactor < 1e-12)) {
      if (bestScale >= 0.0) {
        W[priority] = bestW;
        dotQopt[priority] = higherPriorityJointVelocity
            + bestInvJP * (bestScale * task - jacobian * higherPriorityJointVelocity)
            + bestTildeP * bestDotQn;
      } else {
        dotQopt[priority] = higherPriorityJointVelocity;
      }

      *jointVelocity = dotQopt[priority];
      return bestScale;
    }

  } while (limit_excedeed);
  *jointVelocity = dotQopt[priority];
  return scalingFactor;
}

bool OSNSVelocityIK::isOptimal(int priority, const Eigen::VectorXd& dotQ,
                               const Eigen::MatrixXd& tildeP, Eigen::MatrixXd* W,
                               Eigen::VectorXd* dotQn, double eps) {

  Eigen::VectorXd barMu;
  bool isOptimal = true;

  barMu = tildeP.transpose() * dotQ;

  for (int i = 0; i < n_dof; i++) {
    if ((*W)(i, i) < 0.1) {  //equal to 0.0 (safer)

      if (abs(dotQ(i) - dotQmax(i)) < eps)
        barMu(i) = -barMu(i);

      if (barMu(i) < 0.0) {
        (*W)(i, i) = 1.0;
        (*dotQn)(i) = 0.0;
        isOptimal = false;
      }
    }

  }
  return isOptimal;
}

}  // namespace sns_ik
