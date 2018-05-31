/*! \file fsns_velocity_ik.hpp
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

#ifndef FSNS_IK_VELOCITY_IK
#define FSNS_IK_VELOCITY_IK

#include <Eigen/Dense>

#include "sns_ik/sns_velocity_ik.hpp"

namespace sns_ik {

class FSNSVelocityIK : public SNSVelocityIK {
  public:
    FSNSVelocityIK(int dof, double loop_period);
    virtual ~FSNSVelocityIK() {};

    // Optimal SNS Velocity IK
    virtual double getJointVelocity(Eigen::VectorXd *jointVelocity, const std::vector<Task> &sot,
                  const Eigen::VectorXd &jointConfiguration);

  protected:
    // Perform the SNS for a single task
    virtual double SNSsingle(int priority, const Eigen::VectorXd &higherPriorityJointVelocity,
                  const Eigen::MatrixXd &higherPriorityNull, const Eigen::MatrixXd &jacobian,
                  const Eigen::VectorXd &task, Eigen::VectorXd *jointVelocity, Eigen::MatrixXd *nullSpaceProjector);

    void getTaskScalingFactor(const Eigen::ArrayXd &a,
                  const Eigen::ArrayXd &b,
                  const Eigen::VectorXi &S, double *scalingFactor,
                  int *mostCriticalJoint);

    // TODO: Does this need to be a member variable?
    std::vector<Eigen::VectorXi> S;  //the i-th element is zero if the i-th joint is not saturate, otherwise contains the position in B
};

}  // namespace sns_ik

#endif
