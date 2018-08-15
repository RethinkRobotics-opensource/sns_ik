/*! \file sns_position_ik.hpp
 * \brief Basic SNS Position IK solver
 * \author Forrest Rogers-Marcovitz
 */
/*
 *    Copyright 2016 Rethink Robotics
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

#ifndef SNS_IK_POSITION_IK
#define SNS_IK_POSITION_IK

#include <memory>
#include <Eigen/Dense>
#include <kdl/chain.hpp>
#include <kdl/chainjnttojacsolver.hpp>
#include <kdl/chainfksolverpos_recursive.hpp>
#include <sns_ik/sns_vel_ik_base.hpp>

namespace sns_ik {

// Forward declare SNS Velocity Base Class
class SNSVelocityIK;
class SNSPositionIK {
  public:
    SNSPositionIK(KDL::Chain chain, std::shared_ptr<SNSVelocityIK> velocity_ik, double eps=1e-5);
    ~SNSPositionIK();

    int CartToJnt(const KDL::JntArray& joint_seed,
                  const KDL::Frame& goal_pose,
                  KDL::JntArray* return_joints,
                  const KDL::Twist& bounds=KDL::Twist::Zero())
    { return CartToJnt(joint_seed, goal_pose, KDL::JntArray(0), Eigen::MatrixXd(0,0),
                       std::vector<int>(0), 0.0, return_joints, bounds); }

    int CartToJnt(const KDL::JntArray& joint_seed,
                  const KDL::Frame& goal_pose,
                  const KDL::JntArray& joint_ns_bias,
                  const Eigen::MatrixXd& ns_jacobian,
                  const std::vector<int>& ns_indicies,
                  const double ns_gain,
                  KDL::JntArray* return_joints,
                  const KDL::Twist& bounds=KDL::Twist::Zero());

    // TODO: looks like this would require the KDL solvers to be wrapped in smart pointers
    //void setChain(const KDL::Chain chain);
    KDL::Chain getChain() { return m_chain; }

    inline bool getVelocityIK(std::shared_ptr<SNSVelocityIK>& velocitySolver) {
      velocitySolver = m_ikVelSolver;
      return m_ikVelSolver != NULL;
    }

    inline void setStepSize(double linearMaxStepSize, double angularMaxStepSize){
      m_linearMaxStepSize = linearMaxStepSize;
      m_angularMaxStepSize = angularMaxStepSize;
    }

    inline void setMaxIterations(double maxIterations) {
      m_maxIterations = maxIterations;
    }

    void setDeltaTime(double dt) {
      m_dt = dt;
    }

    void setUseBarrierFunction(bool use) {
      m_useBarrierFunction = use;
    }
    void setBarrierInitAlpha(double alpha) {
      m_barrierInitAlpha = alpha;
    }
    bool setBarrierDecay(double decay) {
      if (decay > 0 && decay <= 1.0) {
        m_barrierDecay = decay;
        return true;
      } else {
        return false;
      }
    }

  private:
    KDL::Chain m_chain;
    std::shared_ptr<SNSVelocityIK> m_ikVelSolver;
    KDL::ChainFkSolverPos_recursive m_positionFK;
    KDL::ChainJntToJacSolver m_jacobianSolver;
    double m_linearMaxStepSize;
    double m_angularMaxStepSize;
    double m_maxIterations;
    double m_eps;
    double m_dt;
    bool m_useBarrierFunction;
    double m_barrierInitAlpha;
    double m_barrierDecay;

    /**
     * @brief Calculate the position and rotation errors in base frame
     * @param q - joints input
     * @param goal - desired goal frame
     * @param pose - pose based of FK of q
     * @param errL - translation error magnitude (== trans.Norm())
     * @param errR - rotational error magnitude (angle-axis representation)
     * @param trans - translation vector
     * @param rotAxis - unit rotation vector
     */
    bool calcPoseError(const KDL::JntArray& q,
                       const KDL::Frame& goal,
                       KDL::Frame* pose,
                       double* errL,
                       double* errR,
                       KDL::Vector* trans,
                       KDL::Vector* rotAxis);

    /**
     * determines the rotation axis necessary to rotate from frame b1 to the
     * orientation of frame b2 and the vector necessary to translate the origin
     * of b1 to the origin of b2, and stores the result in a Twist
     * datastructure.  The result is w.r.t. frame b1.
     * \param F_a_b1 frame b1 expressed with respect to some frame a.
     * \param F_a_b2 frame b2 expressed with respect to some frame a.
     * \warning The result is not a real Twist!
     * \warning In contrast to standard KDL diff methods, the result of
     * diffRelative is w.r.t. frame b1 instead of frame a.
     */
    IMETHOD KDL::Twist diffRelative(const KDL::Frame & F_a_b1, const KDL::Frame & F_a_b2, double dt = 1)
    {
        return KDL::Twist(F_a_b1.M.Inverse() * diff(F_a_b1.p, F_a_b2.p, dt),
                     F_a_b1.M.Inverse() * diff(F_a_b1.M, F_a_b2.M, dt));
    }

};

}  // namespace sns_ik

#endif
