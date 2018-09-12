/*! \file sns_velocity_base_ik.cpp
 * \brief Basic SNS velocity IK solver (Rethink version)
 * \author Fabrizio Flacco
 * \author Forrest Rogers-Marcovitz
 * \author Andy Park
 */
/*
 *    Copyright 2018 Rethink Robotics
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

#include <sns_ik/sns_vel_ik_base_interface.hpp>
#include <ros/console.h>
#include "sns_ik_math_utils.hpp"

namespace sns_ik {

SNSVelIKBaseInterface::SNSVelIKBaseInterface(int dof, double loop_period) :
  SNSVelocityIK(dof, loop_period)
{
  baseIkSolver = SnsVelIkBase::create(n_dof);
}

double SNSVelIKBaseInterface::getJointVelocity(Eigen::VectorXd *jointVelocity,
    const std::vector<Task> &sot,
    const Eigen::VectorXd &jointConfiguration)
{
  // This will only reset member variables if different from previous values
  setNumberOfTasks(sot.size(), sot[0].jacobian.cols());

  // calculate box constraints
  shapeJointVelocityBound(jointConfiguration);

  // solve using SNS base IK solver (andy)
  dqLow = dotQmin;
  dqUpp = dotQmax;
  J = sot[0].jacobian;
  dx = sot[0].desired;
  dqSol.resize(n_dof);
  taskScale = 1.0;

  if (n_tasks > 1) {
    // if the tasks include a nullspace bias task
    taskScaleCS = 1.0;
    dqCS = sot[1].desired;
  }

  // set box constraints
  baseIkSolver->setBounds(dqLow, dqUpp);

  if (n_tasks == 1) {
    exitCode = baseIkSolver->solve(J, dx, &dqSol, &taskScale);
    scaleFactors[0] = taskScale;
  }
  else {
    exitCode = baseIkSolver->solve(J, dx, dqCS, &dqSol, &taskScale, &taskScaleCS);
    scaleFactors[0] = taskScale;
    scaleFactors[1] = taskScaleCS;
  }

  // store solution and scale factor
  *jointVelocity = dqSol;

  // return -1.0 when IK was not successful
  if (exitCode != SnsIkBase::ExitCode::Success)
    return -1.0;

  return 1.0;
}



}  // namespace sns_ik