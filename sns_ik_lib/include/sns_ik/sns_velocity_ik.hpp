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

#ifndef SNS_IK_VELOCITY_IK
#define SNS_IK_VELOCITY_IK

#include <Eigen/Dense>

#include "sns_ik/sns_ik_math_utils.hpp"

using namespace Eigen;

namespace sns_ik {

/*! \struct Task
 *  A desired robot task
 */

struct Task {
    MatrixD jacobian;  //!< the task Jacobian
    VectorD desired;   //!< desired velocity in task space
};

#define _ONLY_WARNING_ON_ERROR
/*! \def _ONLY_WARNING_ON_ERROR
 * with this is possible to avoid to stop the code if some non critical errors appears.
 */

#define SHAPE_MARGIN 0.98

class SNSVelocityIK {
  public:
    SNSVelocityIK(int dof, Scalar loop_period);
    virtual ~SNSVelocityIK() {};

    bool setJointsCapabilities(const VectorD limit_low, const VectorD limit_high,
                               const VectorD maxVelocity, const VectorD maxAcceleration);
    bool setMaxJointVelocity(const VectorD maxVelocity);
    bool setMaxJointAcceleration(const VectorD maxAcceleration);
    virtual void setNumberOfTasks(int ntasks, int dof = -1);
    virtual void setNumberOfDOF(int dof);

    // The control loop period in seconds
    void setLoopPeriod(double period) { loop_period = period; }

    // SNS Velocity IK
    virtual Scalar getJointVelocity(VectorD *jointVelocity, const std::vector<Task> &sot,
                                    const VectorD &jointConfiguration);

    // Standard straight inverse jacobian
    Scalar getJointVelocity_STD(VectorD *jointVelocity, const std::vector<Task> &sot);

    // The standard velocity IK solver doesn't need the joint configuration, but it's here for consistancy
    Scalar getJointVelocity_STD(VectorD *jointVelocity, const std::vector<Task> &sot,
                                const VectorD &jointConfiguration)
        { return getJointVelocity_STD(jointVelocity, sot); }

    std::vector<Scalar> getTasksScaleFactor()
        { return scaleFactors; }

    VectorD getJointLimitLow() { return jointLimit_low; }
    VectorD getJointLimitHigh() { return jointLimit_high; }
    VectorD getJointVelocityMax() { return maxJointVelocity; }

    void usePositionLimits(bool use) { m_usePositionLimits = use; }

  protected:

    // Shape the joint velocity bound dotQmin and dotQmax
    void shapeJointVelocityBound(const VectorD &actualJointConfiguration, double margin = SHAPE_MARGIN);

    // Perform the SNS for a single task
    virtual Scalar SNSsingle(int priority, const VectorD &higherPriorityJointVelocity,
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
    bool m_usePositionLimits;

    Array<Scalar, Dynamic, 1> dotQmin;  // lower joint velocity bound
    Array<Scalar, Dynamic, 1> dotQmax;  // higher joint velocity bound

    // TODO: are these needed here???
    VectorD dotQ;  // next solution (bar{\dotqv} in the paper)
    std::vector<MatrixD> W;  //selection matrices  [here to permit a warm start]
    std::vector<VectorD> dotQopt;  // next solution (bar{\dotqv} in the paper)
    MatrixD I;  // identity matrix
    std::vector<Scalar> scaleFactors;

    std::vector<int> nSat;  //number of saturated joint
};

}  // namespace sns_ikl

#endif
