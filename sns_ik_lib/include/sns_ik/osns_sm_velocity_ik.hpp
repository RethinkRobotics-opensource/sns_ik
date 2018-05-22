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

#ifndef OSNS_SM_IK_VELOCITY_IK
#define OSNS_SM_IK_VELOCITY_IK

#include <Eigen/Dense>

#include "sns_ik/osns_velocity_ik.hpp"

namespace sns_ik {

class OSNS_sm_VelocityIK : public OSNSVelocityIK {
  public:
    OSNS_sm_VelocityIK(int dof, double loop_period);
    virtual ~OSNS_sm_VelocityIK() {};

    // Optimal SNS Velocity IK
    virtual double getJointVelocity(Eigen::VectorXd *jointVelocity, const std::vector<Task> &sot,
                            const Eigen::VectorXd &jointConfiguration);

    void setScaleMargin(double scale)
      { m_scaleMargin = scale; }
    double getScaleMargin()
      { return m_scaleMargin; }

  protected:
    double m_scaleMargin;
};

}  // namespace sns_ik

#endif
