/*! \file sns_velocity_ik.hpp
 * \brief Basic SNS velocity IK solver
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

#ifndef SNS_IKL_VELOCITY_IK
#define SNS_IKL_VELOCITY_IK

#include <Eigen/Dense>

#include "sns_ik/sns_ik_math_utils.hpp"

using namespace Eigen;

namespace sns_ik {

/*! \struct Task
 *  A desired robot task
 *  \todo for now we are not able to work in acceleration becouse we need the derivative of the Jacobian, but we should consider also the
 *  pre-Jacobian, post-Jacobian and augmented_Jacobian in the derivate... how can we do that?
 */

struct Task {
    MatrixD jacobian;  //!< the task Jacobian
    MatrixD jacobian_derivate;  //!< time derivartive of the Jacobian, needed only at acceleration control level
    VectorD desired;  //!< desired velocity/acceleration
    Scalar h;
    VectorD q0;
    Scalar h0;
    Scalar min;  //!< min for the constraints
    Scalar max;  //!< max for the constraints
    Scalar value;  //!< value for the constraints
    int is_constraints;
};

/*! \typedef StackOfTasks
 *  the stack of desired tasks
 */
typedef std::vector<Task> StackOfTasks;

#define _ONLY_WARNING_ON_ERROR
/*! \def _ONLY_WARNING_ON_ERROR
 * with this is possible to avoid to stop the code if some non critical errors appears.
 */

#define SHAPE_MARGIN 0.98

class SNSVelocityIK {
  public:
    SNSVelocityIK(int dof, Scalar loop_period);
    
    bool setJointsCapabilities(VectorD limit_low, VectorD limit_high,
                               VectorD maxVelocity, VectorD maxAcceleration);
    void setNumberOfTasks(int ntasks, int dof = -1);
    void setNumberOfDOF(int dof);

    // The control loop period in seconds
    void setLoopPeriod(double period) { loop_period = period; }
    
    // Standard straight inverse jacobian
    Scalar getJointVelocity(VectorD *jointVelocity, const StackOfTasks &sot,
                            const VectorD &jointConfiguration);

    // Standard straight inverse jacobian
    Scalar getJointVelocity_STD(VectorD *jointVelocity, const StackOfTasks &sot);

    // The standard velocity IK solver doesn't need the joint configuration, but it's here for consistancy
    Scalar getJointVelocity_STD(VectorD *jointVelocity, const StackOfTasks &sot,
                            const VectorD &jointConfiguration)
        { return getJointVelocity_STD(jointVelocity, sot); }

    std::vector<Scalar> getTasksScaleFactor()
        { return scaleFactors; }

  private:

    // Shape the joint velocity bound dotQmin and dotQmax
    void shapeJointVelocityBound(const VectorD &actualJointConfiguration, double margin = SHAPE_MARGIN);

    // Perform the SNS for a single task
    Scalar SNSsingle(int priority, const VectorD &higherPriorityJointVelocity,
                     const MatrixD &higherPriorityNull, const MatrixD &jacobian,
                     const VectorD &task, VectorD *jointVelocity, MatrixD *nullSpaceProjector);

    void getTaskScalingFactor(const Array<Scalar, Dynamic, 1> &a,
                              const Array<Scalar, Dynamic, 1> &b,
                              const MatrixD &W, Scalar *scalingFactor,
                              int *mostCriticalJoint);

    int n_dof;  //manipulator degree of freedom
    int n_tasks;  //number of tasks
    Scalar loop_period;  //needed to compute the bounds

    VectorD jointLimit_low;  // low joint limits
    VectorD jointLimit_high;  // high joint limit
    VectorD maxJointVelocity;  // maximum joint velocity
    VectorD maxJointAcceleration;  // maximum joint acceleration

    Array<Scalar, Dynamic, 1> dotQmin;  // lower joint velocity bound
    Array<Scalar, Dynamic, 1> dotQmax;  // higher joint velocity bound

    // TODO: are these needed here???
    VectorD dotQ;  // next solution (bar{\dotqv} in the paper)
    std::vector<MatrixD> W;  //selection matrices  [here to permit a warm start]
    std::vector<VectorD> dotQopt;  // next solution (bar{\dotqv} in the paper)
    MatrixD I;  // identity matrix
    std::vector<Scalar> scaleFactors;
};

}  // namespace sns_ikl

#endif
