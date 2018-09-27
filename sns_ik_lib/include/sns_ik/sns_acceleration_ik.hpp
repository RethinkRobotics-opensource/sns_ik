/*! \file sns_acceleration_ik.hpp
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

#ifndef SNS_IK_ACCELERATION_IK
#define SNS_IK_ACCELERATION_IK

#include <Eigen/Dense>
#include <vector>
#include <sns_ik/sns_acc_ik_base.hpp>

namespace sns_ik {

/*! \struct Task
 *  A desired robot task
 */

struct TaskAcc {
    Eigen::MatrixXd jacobian;  //!< the task Jacobian
    Eigen::VectorXd dJdq;  //!< the product of task Jacobian dot and joint velocity
    Eigen::VectorXd desired;   //!< desired velocity in task space
};

class SNSAccelerationIK {
  public:
    SNSAccelerationIK(int dof, double loop_period);
    virtual ~SNSAccelerationIK() {};

    bool setJointsCapabilities(const Eigen::VectorXd limit_low, const Eigen::VectorXd limit_high,
                               const Eigen::VectorXd maxVelocity, const Eigen::VectorXd maxAcceleration);
    bool setMaxJointAcceleration(const Eigen::VectorXd maxAcceleration);
    virtual void setNumberOfTasks(int ntasks, int dof = -1);
    virtual void setNumberOfDOF(int dof);

    // The control loop period in seconds
    void setLoopPeriod(double period) { loop_period = period; }

    // SNS Acceleration IK
    int getJointAcceleration(Eigen::VectorXd *jointAcceleration, const std::vector<TaskAcc> &sot,
                                    const Eigen::VectorXd &jointConfiguration, const Eigen::VectorXd &jointVelocities);

    std::vector<double> getTasksScaleFactor()
        { return scaleFactors; }

    Eigen::VectorXd getJointLimitLow() { return jointLimit_low; }
    Eigen::VectorXd getJointLimitHigh() { return jointLimit_high; }
    Eigen::VectorXd getJointVelocityMax() { return maxJointVelocity; }

    void usePositionLimits(bool use) { m_usePositionLimits = use; }

  protected:

    // Shape the joint acceleration bound ddotQmin and ddotQmax
    void shapeJointAccelerationBound(const Eigen::VectorXd &actualJointConfiguration,
                        const Eigen::VectorXd &actualJointVelocities, double margin = 0.98);

    int n_dof;  //manipulator degree of freedom
    int n_tasks;  //number of tasks
    double loop_period;  //needed to compute the bounds

    Eigen::VectorXd jointLimit_low;  // low joint limits
    Eigen::VectorXd jointLimit_high;  // high joint limit
    Eigen::VectorXd maxJointVelocity;  // maximum joint velocity
    Eigen::VectorXd maxJointAcceleration;  // maximum joint acceleration
    bool m_usePositionLimits;

    Eigen::ArrayXd ddotQmin;  // lower joint velocity bound
    Eigen::ArrayXd ddotQmax;  // higher joint velocity bound
    std::vector<double> scaleFactors;

    // variables for base acc ik solver
    Eigen::ArrayXd ddqLow;
    Eigen::ArrayXd ddqUpp;
    Eigen::MatrixXd J;
    Eigen::VectorXd dJdq;
    Eigen::VectorXd ddx;
    Eigen::VectorXd ddqCS;
    Eigen::VectorXd ddqSol;

    double taskScale, taskScaleCS;

    SnsAccIkBase::uPtr baseIkSolver;

    SnsIkBase::ExitCode exitCode;
};

}  // namespace sns_ik

#endif
