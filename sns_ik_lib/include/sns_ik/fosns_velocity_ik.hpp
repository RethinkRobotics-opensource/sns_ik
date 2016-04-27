/*! \file fosns_velocity_ik.hpp
 * \brief Fast Optimal SNS velocity IK solver
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

#ifndef FOSNS_IK_VELOCITY_IK
#define FOSNS_IK_VELOCITY_IK

#include <Eigen/Dense>
#include <forward_list>

#include "sns_ik/sns_ik_math_utils.hpp"
#include "sns_ik/sns_velocity_ik.hpp"
#include "sns_ik/fsns_velocity_ik.hpp"

using namespace Eigen;

namespace sns_ik {

class FOSNSVelocityIK : public FSNSVelocityIK {
  public:
    FOSNSVelocityIK(int dof, Scalar loop_period);
    virtual ~FOSNSVelocityIK() {};

    // Optimal SNS Velocity IK
    virtual Scalar getJointVelocity(VectorD *jointVelocity, const StackOfTasks &sot,
                  const VectorD &jointConfiguration);

    virtual void setNumberOfTasks(int ntasks, int dof);

    void setScaleMargin(double scale)
      { scaleMargin = scale; }
    double getScaleMargin()
      { return scaleMargin; }

  protected:
    double scaleMargin;

    // For the FastOpt version of the SNS
    // TODO: should these be member variables?
    MatrixD B;  //update matrix
    std::vector<std::forward_list<int>> satList;
    VectorD lagrangeMu;
    VectorD lagrangeMu1;
    VectorD lagrangeMup2w;
    std::forward_list<int>::iterator it, prev_it;

    // Perform the SNS for a single task
    virtual Scalar SNSsingle(int priority, const VectorD &higherPriorityJointVelocity,
                   const MatrixD &higherPriorityNull, const MatrixD &jacobian,
                   const VectorD &task, VectorD *jointVelocity, MatrixD *nullSpaceProjector);
};

}  // namespace sns_ikl

#endif
