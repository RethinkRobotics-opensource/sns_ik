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

#include "sns_ik/sns_ik_math_utils.hpp"
#include "sns_ik/sns_velocity_ik.hpp"

//#define LOG_ACTIVE

using namespace Eigen;

namespace sns_ik {

class FSNSVelocityIK : public SNSVelocityIK {
  public:
    FSNSVelocityIK(int dof, Scalar loop_period);
    virtual ~FSNSVelocityIK() {};

    // Optimal SNS Velocity IK
    virtual Scalar getJointVelocity(VectorD *jointVelocity, const std::vector<Task> &sot,
                  const VectorD &jointConfiguration);

  protected:
    // Perform the SNS for a single task
    virtual Scalar SNSsingle(int priority, const VectorD &higherPriorityJointVelocity,
                  const MatrixD &higherPriorityNull, const MatrixD &jacobian,
                  const VectorD &task, VectorD *jointVelocity, MatrixD *nullSpaceProjector);

    void getTaskScalingFactor(const Array<Scalar, Dynamic, 1> &a,
                  const Array<Scalar, Dynamic, 1> &b,
                  const VectorXi &S, Scalar *scalingFactor,
                  int *mostCriticalJoint);

    // TODO: Does this need to be a member variable?
    std::vector<VectorXi> S;  //the i-th element is zero if the i-th joint is not saturate, otherwise contains the position in B
};

}  // namespace sns_ikl

#endif
