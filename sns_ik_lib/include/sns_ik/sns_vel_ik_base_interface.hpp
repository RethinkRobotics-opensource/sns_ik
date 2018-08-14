/*! \file sns_velocity_base_ik.hpp
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

#ifndef SNS_IK_VELOCITY_BASE_IK
#define SNS_IK_VELOCITY_BASE_IK

#include <Eigen/Dense>
#include "sns_ik/sns_velocity_ik.hpp"
#include "sns_vel_ik_base.hpp"

namespace sns_ik {

class SNSVelIKBaseInterface : public SNSVelocityIK {
  public:
    SNSVelIKBaseInterface(int dof, double loop_period);
    virtual ~SNSVelIKBaseInterface() {};

    // Optimal SNS Velocity IK
    virtual double getJointVelocity(Eigen::VectorXd *jointVelocity, const std::vector<Task> &sot,
                  const Eigen::VectorXd &jointConfiguration);

  protected:
   Eigen::ArrayXd dqLow;
   Eigen::ArrayXd dqUpp;
   Eigen::MatrixXd J;
   Eigen::VectorXd dx;
   Eigen::VectorXd dqCS;
   Eigen::VectorXd dqSol;

   double taskScale, taskScaleCS;

   SnsVelIkBase::uPtr baseIkSolver;

   SnsIkBase::ExitCode exitCode;
};

}  // namespace sns_ik

#endif
