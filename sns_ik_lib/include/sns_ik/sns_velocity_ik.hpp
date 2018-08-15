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
#include <vector>
#include <sns_ik/sns_vel_ik_base.hpp>

namespace sns_ik {

/*! \struct Task
 *  A desired robot task
 */

struct Task {
    Eigen::MatrixXd jacobian;  //!< the task Jacobian
    Eigen::VectorXd desired;   //!< desired velocity in task space
};

static const double SHAPE_MARGIN = 0.98;

class SNSVelocityIK {
  public:
    SNSVelocityIK(int dof, double loop_period);
    virtual ~SNSVelocityIK() {};

    bool setJointsCapabilities(const Eigen::VectorXd limit_low, const Eigen::VectorXd limit_high,
                               const Eigen::VectorXd maxVelocity, const Eigen::VectorXd maxAcceleration);
    bool setMaxJointVelocity(const Eigen::VectorXd maxVelocity);
    bool setMaxJointAcceleration(const Eigen::VectorXd maxAcceleration);
    virtual void setNumberOfTasks(int ntasks, int dof = -1);
    virtual void setNumberOfDOF(int dof);

    // The control loop period in seconds
    void setLoopPeriod(double period) { loop_period = period; }

    // SNS Velocity IK
    virtual double getJointVelocity(Eigen::VectorXd *jointVelocity, const std::vector<Task> &sot,
                                    const Eigen::VectorXd &jointConfiguration);

    // Standard straight inverse jacobian
    double getJointVelocity_STD(Eigen::VectorXd *jointVelocity, const std::vector<Task> &sot);

    // The standard velocity IK solver doesn't need the joint configuration, but it's here for consistancy
    double getJointVelocity_STD(Eigen::VectorXd *jointVelocity, const std::vector<Task> &sot,
                                const Eigen::VectorXd &jointConfiguration)
        { return getJointVelocity_STD(jointVelocity, sot); }

    std::vector<double> getTasksScaleFactor()
        { return scaleFactors; }

    Eigen::VectorXd getJointLimitLow() { return jointLimit_low; }
    Eigen::VectorXd getJointLimitHigh() { return jointLimit_high; }
    Eigen::VectorXd getJointVelocityMax() { return maxJointVelocity; }

    void usePositionLimits(bool use) { m_usePositionLimits = use; }

  protected:

    // Shape the joint velocity bound dotQmin and dotQmax
    void shapeJointVelocityBound(const Eigen::VectorXd &actualJointConfiguration, double margin = SHAPE_MARGIN);

    // Perform the SNS for a single task
    virtual double SNSsingle(int priority, const Eigen::VectorXd &higherPriorityJointVelocity,
                     const Eigen::MatrixXd &higherPriorityNull, const Eigen::MatrixXd &jacobian,
                     const Eigen::VectorXd &task, Eigen::VectorXd *jointVelocity, Eigen::MatrixXd *nullSpaceProjector);

    void getTaskScalingFactor(const Eigen::ArrayXd &a,
                              const Eigen::ArrayXd &b,
                              const Eigen::MatrixXd &W, double *scalingFactor,
                              int *mostCriticalJoint);

    int n_dof;  //manipulator degree of freedom
    int n_tasks;  //number of tasks
    double loop_period;  //needed to compute the bounds

    Eigen::VectorXd jointLimit_low;  // low joint limits
    Eigen::VectorXd jointLimit_high;  // high joint limit
    Eigen::VectorXd maxJointVelocity;  // maximum joint velocity
    Eigen::VectorXd maxJointAcceleration;  // maximum joint acceleration
    bool m_usePositionLimits;

    Eigen::ArrayXd dotQmin;  // lower joint velocity bound
    Eigen::ArrayXd dotQmax;  // higher joint velocity bound

    // TODO: are these needed here???
    Eigen::VectorXd dotQ;  // next solution (bar{\dotqv} in the paper)
    std::vector<Eigen::MatrixXd> W;  //selection matrices  [here to permit a warm start]
    std::vector<Eigen::VectorXd> dotQopt;  // next solution (bar{\dotqv} in the paper)
    Eigen::MatrixXd I;  // identity matrix
    std::vector<double> scaleFactors;

    std::vector<int> nSat;  //number of saturated joint
};

}  // namespace sns_ik

#endif
