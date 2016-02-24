/*! \file osns_velocity_ik.hpp
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

#ifndef OSNS_IK_VELOCITY_IK
#define OSNS_IK_VELOCITY_IK

#include <Eigen/Dense>

#include "sns_ik/sns_ik_math_utils.hpp"
#include "sns_ik/sns_velocity_ik.hpp"

using namespace Eigen;

namespace sns_ik {

class OSNSVelocityIK : public SNSVelocityIK {
  public:
    OSNSVelocityIK(int dof, Scalar loop_period);
    
    // Optimal SNS Velocity IK
    Scalar getJointVelocity(VectorD *jointVelocity, const StackOfTasks &sot,
                            const VectorD &jointConfiguration);
    
  protected:
    // Perform the SNS for a single task
    Scalar SNSsingle(int priority, const VectorD &higherPriorityJointVelocity,
                     const MatrixD &higherPriorityNull, const MatrixD &jacobian,
                     const VectorD &task, VectorD *jointVelocity, MatrixD *nullSpaceProjector);
    
    bool isOptimal(int priority, const VectorD& dotQ, 
                   const MatrixD& tildeP, MatrixD* W, 
                   VectorD* dotQn, double eps = 1e-8);
};

}  // namespace sns_ikl

#endif
