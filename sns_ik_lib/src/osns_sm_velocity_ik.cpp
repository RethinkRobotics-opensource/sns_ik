/*! \file osns_sm_velocity_ik.hpp
 * \brief Optimal SNS, with scale margin, velocity IK solver
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

#include <sns_ik/osns_sm_velocity_ik.hpp>

using namespace Eigen;
using namespace sns_ik;

OSNS_sm_VelocityIK::OSNS_sm_VelocityIK(int dof, Scalar loop_period) :
  OSNSVelocityIK(dof, loop_period),
  m_scaleMargin(0.9)
{
}


Scalar OSNS_sm_VelocityIK::getJointVelocity(VectorD *jointVelocity,
    const StackOfTasks &sot,
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

  shapeJointVelocityBound(jointConfiguration, m_scaleMargin);

  MatrixD PS = MatrixD::Identity(n_dof, n_dof);

  for (int i_task = 0; i_task < n_tasks; i_task++) {  //consider all tasks
    higherPriorityJointVelocity = *jointVelocity;
    higherPriorityNull = P;

    scaleFactors[i_task] = SNSsingle(i_task, higherPriorityJointVelocity, higherPriorityNull,
        sot[i_task].jacobian, sot[i_task].desired, jointVelocity, &PS);

    if (scaleFactors[i_task] < 0) {
      //second chance
      W[i_task] = I;
      PS = higherPriorityNull;
      scaleFactors[i_task] = SNSsingle(i_task, higherPriorityJointVelocity, higherPriorityNull,
          sot[i_task].jacobian, sot[i_task].desired, jointVelocity, &PS);

    }

    if (scaleFactors[i_task] > 0.0) {
      if (scaleFactors[i_task] * m_scaleMargin < (1.0)) {
        double taskScale = scaleFactors[i_task] * m_scaleMargin;
        VectorD scaledTask = sot[i_task].desired * taskScale;
        scaleFactors[i_task] = SNSsingle(i_task, higherPriorityJointVelocity, higherPriorityNull,
            sot[i_task].jacobian, scaledTask, jointVelocity, &P);
        scaleFactors[i_task] = taskScale;

      } else {
        scaleFactors[i_task] = 1.0;
        P = PS;
      }
    }
  }

  return 1.0;
}
