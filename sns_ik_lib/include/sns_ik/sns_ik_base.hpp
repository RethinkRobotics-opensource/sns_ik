/** @file sns_ik_base.hpp
 *
 * @brief This file provides the common functions shared by SNS IK base solvers
 *
 * @author Matthew Kelly
 * @author Andy Park
 *
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

#ifndef SNS_IK_LIB__SNS_IK_BASE_H_
#define SNS_IK_LIB__SNS_IK_BASE_H_

#include <Eigen/Dense>
#include <memory>
#include <vector>

#include "sns_linear_solver.hpp"
#include "sns_ik_math_utils.hpp"

namespace sns_ik {

/*
 * This class is an abstract base class that is used by all of the SNS-IK solvers.
 *
 * Note: throughout this class documentation we use the variable q to represent the configuration
 *       space variables. For example, a position solver will use q for position, a velocity solver
 *       will use q for velocity, and an acceleration solver will use q for acceleration.
 */
class SnsIkBase {

public:

  // Smart pointer typedefs. Note: all derived classes MUST override these smart pointers.
  typedef std::shared_ptr<SnsIkBase> Ptr;
  typedef std::unique_ptr<SnsIkBase> uPtr;

  enum class ExitCode {
    Success,  // successfully solver the optimization problem
    BadUserInput,  // user input failed basic checks (eg. inconsistent matrix size)
    InfeasibleTask,  // there is no feasible solution to the primary task
    InternalError  // there was an internal error in the solver (should never happen...)
  };

  // Make sure that class is cleaned-up correctly
  virtual ~SnsIkBase() {};

  /**
   * Set the bounds on the configuration space variable. These limits match the order of
   * the solver, for example, on a velocity solver this will set the velocity limits.
   * Set a bound to infinity (or an arbitrarily large value) to disable it.
   * This method will update the number of joints in the solver based on the size of qLow and qUpp
   * Requirements: inputs must be the same size and qLow < qUpp
   * @param qLow: lower bound on the velocity/acceleration of each joint
   * @param qUpp: upper bound on the velocity/acceleration of each joint
   * @return: true iff successful
   */
  bool setBounds(const Eigen::ArrayXd& qLow, const Eigen::ArrayXd& qUpp);

  /*
   * @return: reference to the lower bound in configuration space (pos / vel / acc)
   */
  const Eigen::ArrayXd& getLowerBounds() const { return qLow_; };

  /*
   * @return: reference to the lower bound in configuration space (pos / vel / acc)
   */
  const Eigen::ArrayXd& getUpperBounds() const { return qUpp_; };

  /*
   * @return: number of joints
   */
  size_t getNrOfJoints() const { return nJnt_; }

protected:

  /*
   * The code of the SNS-IK solver relies on a linear solver. If the linear system is infeasible,
   * then the solver will return the minimum-norm solution with a non-zero (positive) residual.
   * This corresponds to a situation where there is no feasible solution to the task objective, for
   * example, a one-link robot arm where the task velocity is along the axis of the joint. This
   * threshold is used to determine if the task is infeasible.
   */
  static const double LIN_SOLVE_RESIDUAL_TOL;

  // A tolerance used for SVD-based pseudo-inverse operation
  static const double PINV_TOL;

  /*
   * If all goes well, then the solver should always terminate the main loop. If something crazy
   * happens, then we want to prevent an infinite loop. The maximum iteration count is given by
   * this factor multiplied by the number of joints in the task jacobian.
   */
  static const int MAXIMUM_SOLVER_ITERATION_FACTOR;

  // Short names for useful constants
  static const double POS_INF;
  static const double NEG_INF;

  // Any number smaller than this is considered to be zero for numerical purposes
  static const double MINIMUM_FINITE_SCALE_FACTOR;

  // Any number larger than this is considered to be infinite for numerical purposes
  static const double MAXIMUM_FINITE_SCALE_FACTOR;

  // Tolerance for checks on the velocity/acceleration limits
  static const double BOUND_TOLERANCE;

  // Nice formatting option from Eigen
  static const Eigen::IOFormat EigArrFmt;

  /*
   * protected constructor: require factory method to create an object.
   */
  SnsIkBase(int nJnt) : nJnt_(nJnt), qLow_(nJnt), qUpp_(nJnt) {};

  /*
   * Check that qLow_ <= q <= qUpp_
   * @param qUpp: joint acceleration to test
   * @return: true iff qLow <= q <= qUpp
   *          if q.empty() return false
   */
  bool checkBounds(const Eigen::VectorXd& q);

  /*
   * This method sets and solves the decomposition of the matrix that is used by the linear solver.
   * In general, this will be the J*W matrix, which describes the jacobian of the active joints.
   * @param JW: matrix to set in the linear solver.
   * @return: Success if the decomposition was successful
   */
  SnsIkBase::ExitCode setLinearSolver(const Eigen::MatrixXd& JW);

  /*
   * Solve a specific linear system and compute the residual error. Linear system:
   * JW * q = rhs
   * @param rhs: "right hand side" of the linear system.
   * @param[out] q: solution to the linear system
   * @param[out] resErr: residual error in the linear system
   * @return: Success if the solve was successful
   */
  SnsIkBase::ExitCode solveLinearSystem(const Eigen::MatrixXd& rhs, Eigen::VectorXd* q, double* resErr);

  /*
   * @return: rank of the matrix that is currently set in the linear solver
   */
  unsigned int getLinSolverRank() const { return linSolver_.rank(); }

  /*
   * Solve the following equation for the variable qUpp:
   *    J * W * (qUpp - dqNull) = dx - J*dqNull
   * This equation appears throughout the SNS-IK papers. One example is in block "D" of Algorithm 1:
   *  "Control of Redundant Robots Under Hard Joint Constraint: Saturation in the Null Space"
   *   by: Fabrizio Flacco, Alessandro De Luca, Oussama Khatib
   *
   * PRECONDITION:
   * --> Assumes that setLinearSolver(J*W) has been successfully called
   * --> This allows multiple projections to be solved for a single matrix decomposition
   *
   * @param J: Jacobian matrix, mapping from joint to task space. Size = [nTask, nJoint]
   * @param dqNull: null-space joint acceleration. Size = nJoint
   * @param dx: task acceleration vector. Length = nTask
   * @param[out] q: joint velocity/acceleration solution. Length = nJoint
   * @param[out, opt] resErr: residual error (norm-squared)
   * @return: ExitCode::Success: the algorithm worked correctly and satisfied the problem statement
   *          otherwise: something went wrong, exit code specifics the type of problem
   */
  SnsIkBase::ExitCode solveProjectionEquation(const Eigen::MatrixXd& J,
                               const Eigen::VectorXd& dqNull, const Eigen::VectorXd& dx,
                               Eigen::VectorXd* q, double* resErr);

  /*
   * Solve the following equation for the variable ddqUpp:
   *    J * W * (ddq - dqNull) = ddx - dJdq - J*ddqNull
   * This equation appears throughout the SNS-IK papers. One example is in block "D" of Algorithm 1:
   *  "Control of Redundant Robots Under Hard Joint Constraint: Saturation in the Null Space"
   *   by: Fabrizio Flacco, Alessandro De Luca, Oussama Khatib
   *
   * PRECONDITION:
   * --> Assumes that linSolver_ has been initialized with J*W
   *
   * @param J: Jacobian matrix, mapping from joint to task space. Size = [nTask, nJoint]
   * @param dJdq: the product of Jacobian derivative and joint velocity. Length = nTask
   * @param ddqNull: null-space joint acceleration. Size = nJoint
   * @param ddx: task acceleration vector. Length = nTask
   * @param[out] ddq: joint acceleration solution. Length = nJoint
   * @param[out, opt] resErr: residual error (norm-squared)
   * @return: ExitCode::Success: the algorithm worked correctly and satisfied the problem statement
   *          otherwise: something went wrong, exit code specifics the type of problem
   */
  SnsIkBase::ExitCode solveProjectionEquation(const Eigen::MatrixXd& J, const Eigen::VectorXd& dJdq,
                               const Eigen::VectorXd& ddqNull, const Eigen::VectorXd& ddx,
                               Eigen::VectorXd* ddq, double* resErr);

  /*
   * This method implements Algorithm 2 (and a bit of Algorithm 1) from the paper:
   *  "Control of Redundant Robots Under Hard Joint Constraint: Saturation in the Null Space"
   *   by: Fabrizio Flacco, Alessandro De Luca, Oussama Khatib
   *
   * PRECONDITION:
   * --> Assumes that linSolver_ has been initialized with J*W
   *
   * @param J: Jacobian matrix, mapping from joint to task space. Size = [nTask, nJoint]
   * @param desiredTask: task velocity/acceleration vector. Length = nTask
   * @param jointOut: joint velocity/acceleration. Length = nJoint
   * @param jntIsFree: which joints are free to saturate? Length = nJoint
   * @param[out] taskScale: task scale factor
   * @param[out] jntIdx: index corresponding to the most critical joint that is free
   * @param[out] resErr: residual error (norm-squared) in the linear solve
   * @return: ExitCode::Success: the algorithm worked correctly and satisfied the problem statement
   *          otherwise: something went wrong, exit code specifics the type of problem
   */
  ExitCode computeTaskScalingFactor(const Eigen::MatrixXd& J, const Eigen::VectorXd& desiredTask, 
                                const Eigen::VectorXd& jointOut, const std::vector<bool>& jntIsFree,
                                 double* taskScale, int* jntIdx, double* resErr);

  /*
   * This algorithm computes the scale factor that is associated with a given joint, but considering
   * both the sensativity of the joint (a) and the distance to the upper and lower limits.
   * @param low: lower margin
   * @param upp: upper margin
   * @param a: margin scale factor
   * @return: joint scale factor
   */
  static double findScaleFactor(double low, double upp, double a);


private:

  int nJnt_; //!< number of joints

  Eigen::ArrayXd qLow_;  //!< lower bound on joint velocity/acceleration
  Eigen::ArrayXd qUpp_;  //!< upper bound on joint velocity/acceleration

  SnsLinearSolver linSolver_;  //!< linear solver for the core SNS-IK algorithm

  Eigen::MatrixXd JW_;  //!< the matrix that is currently set in the linear solver

};  // class SnsIkBase

}  // namespace sns_ik

#endif  // SNS_IK_LIB__SNS_IK_BASE_H_
