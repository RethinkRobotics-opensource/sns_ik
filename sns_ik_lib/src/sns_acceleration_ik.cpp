/*! \file sns_acceleration_ik.cpp
 * \brief Basic SNS acceleration IK solver
 * \author Andy Park
 */
/*
 *    Copyright 2018 Rethink Robotics
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

#include <sns_ik/sns_acceleration_ik.hpp>

#include <ros/console.h>

#include "sns_ik_math_utils.hpp"

namespace sns_ik {

static const Eigen::IOFormat EigArrFmt(4, 0, ", ", "\n", "[", "]");

SNSAccelerationIK::SNSAccelerationIK(int dof, double loop_period) :
  n_dof(0),
  n_tasks(0),
  m_usePositionLimits(true)
{
  setNumberOfDOF(dof);
  setLoopPeriod(loop_period);

  // initialize the sns acc ik solver
  baseIkSolver = SnsAccIkBase::create(n_dof);
}

void SNSAccelerationIK::setNumberOfDOF(int dof)
{
  if (dof > 0 && dof != n_dof) {
    n_dof = dof;
  }
}

void SNSAccelerationIK::setNumberOfTasks(int ntasks, int dof)
{
  setNumberOfDOF(dof);

  if (n_tasks != ntasks) {
    n_tasks = ntasks;
    Eigen::VectorXd ddq = Eigen::VectorXd::Zero(n_dof);

    double scale = 1.0;
    scaleFactors.resize(n_tasks, scale);
  }
}


bool SNSAccelerationIK::setJointsCapabilities(const Eigen::VectorXd limit_low, const Eigen::VectorXd limit_high,
                                          const Eigen::VectorXd maxVelocity, const Eigen::VectorXd maxAcceleration)
{
  if (limit_low.rows() != n_dof || limit_high.rows() != n_dof
      || maxVelocity.rows() != n_dof || maxAcceleration.rows() != n_dof) {
    return false;
  }

  jointLimit_low = limit_low;
  jointLimit_high = limit_high;
  maxJointVelocity = maxVelocity;
  maxJointAcceleration = maxAcceleration;

  ddotQmin = -maxJointAcceleration.array();
  ddotQmax = maxJointAcceleration.array();

  return true;
}

bool SNSAccelerationIK::setMaxJointAcceleration(const Eigen::VectorXd maxAcceleration)
{
  if (maxAcceleration.rows() != n_dof) {
    return false;
  }

  maxJointAcceleration = maxAcceleration;
  ddotQmin = -maxAcceleration.array();
  ddotQmax = maxAcceleration.array();

  return true;
}

// The box constraints (ddotQmin and ddotQmax) in acceleration level are calculated 
// from joint position, velocity and acceleration limits based on the equations 
// in the following paper: 
// Flacco, Fabrizio, Alessandro De Luca, and Oussama Khatib. 
// "Motion control of redundant robots under joint constraints: 
// Saturation in the null space." In IEEE International Conference on Robotics 
// and Automation (ICRA), pp. 285-292, 2012.
void SNSAccelerationIK::shapeJointAccelerationBound(const Eigen::VectorXd &actualJointConfiguration,
                                        const Eigen::VectorXd &actualJointVelocities, double margin) {

  // TODO: rewrite this using the Eigen::Array
  double step, max, stop;

  for (int i = 0; i < n_dof; i++) {
    // for the minimum bound
    max = -maxJointAcceleration(i);
    if (m_usePositionLimits) {
      step = 2 * (jointLimit_low(i) - actualJointConfiguration(i) - actualJointVelocities(i)*loop_period) / (loop_period * loop_period);
      stop = (-maxJointVelocity(i) - actualJointVelocities(i)) / loop_period;
      // take the maximum
      ddotQmin(i) = std::max({step, max, stop});
    } else {
      ddotQmin(i) = max;
    }

    // for the maximum bound
    max = maxJointAcceleration(i);
    if (m_usePositionLimits) {
      step = 2 * (jointLimit_high(i) - actualJointConfiguration(i) - actualJointVelocities(i)*loop_period) / (loop_period * loop_period);
      stop = (maxJointVelocity(i) - actualJointVelocities(i)) / loop_period;
      // take the minimum
      ddotQmax(i) = std::min({step, max, stop});
    } else {
      ddotQmax(i) = max;
    }

    // check whether the acceleration limits are correctly shaped
    if (ddotQmin(i) > ddotQmax(i)) {
      if (ddotQmin(i) > 0) {
      ROS_ERROR_STREAM("Lower limit has been reversed for J"<< i <<"!\n"<< "jointLimit_Low: " << jointLimit_low(i) 
      << ",  actualJointConfiguration: " << actualJointConfiguration(i) <<"\nminJointVelocity: " << -maxJointVelocity(i)
      << ",  actualJointVelocities: " << actualJointVelocities(i));
      }
      if (ddotQmax(i) < 0) {
        ROS_ERROR_STREAM("Upper limit has been reversed for J"<< i <<"!\n"<< "jointLimit_high: " << jointLimit_high(i)
        << ",  actualJointConfiguration: " << actualJointConfiguration(i) << "\nmaxJointVelocity: " << maxJointVelocity(i)
        << ",  actualJointVelocities: " << actualJointVelocities(i));
      }
    }
  }

  ddotQmin *= margin;
  ddotQmax *= margin;
}

int SNSAccelerationIK::getJointAcceleration(Eigen::VectorXd *jointAcceleration,
    const std::vector<TaskAcc> &sot, const Eigen::VectorXd &jointConfiguration, const Eigen::VectorXd &jointVelocities)
{
  // check whether tasks are defined or not
  if (sot.size()==0) {
    ROS_ERROR("Tasks have not been correctly set!");
    return -1;
  } 

  // This will only reset member variables if different from previous values
  setNumberOfTasks(sot.size(), sot[0].jacobian.cols());

  // calculate box constraints
  shapeJointAccelerationBound(jointConfiguration, jointVelocities);

  // define variables
  Eigen::ArrayXd ddqLow = ddotQmin;
  Eigen::ArrayXd ddqUpp = ddotQmax;
  Eigen::MatrixXd J = sot[0].jacobian;
  Eigen::VectorXd dJdq = sot[0].dJdq;
  Eigen::VectorXd ddx = sot[0].desired;
  Eigen::VectorXd ddqSol;
  Eigen::VectorXd ddqCS;
  double taskScale, taskScaleCS;
  SnsIkBase::ExitCode exitCode;

  // set default values 
  ddqSol.resize(n_dof);
  taskScale = 1.0;
  taskScaleCS = 1.0;

  if (n_tasks > 1) {
    // if the tasks include a nullspace bias task
    ddqCS = sot[1].desired;
  }

  // set box constraints
  baseIkSolver->setBounds(ddqLow, ddqUpp);

  if (n_tasks == 1) {
    exitCode = baseIkSolver->solve(J, dJdq, ddx, &ddqSol, &taskScale);
    scaleFactors[0] = taskScale;
  }
  else {
    exitCode = baseIkSolver->solve(J, dJdq, ddx, ddqCS, &ddqSol, &taskScale, &taskScaleCS);
    scaleFactors[0] = taskScale;
    scaleFactors[1] = taskScaleCS;
  }

  // store solution and scale factor
  *jointAcceleration = ddqSol;

  if (exitCode != SnsIkBase::ExitCode::Success)
    return -1;

  return 1;
}

}  // namespace sns_ik
