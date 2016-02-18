/*! \file sns_ikl.hpp
 * \brief The SNS Inverse Kinematic Library (\b IKL)
 * \author Fabrizio Flacco
 */
/*
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

#ifndef SNS_IKL
#define SNS_IKL

#include <Eigen/Dense>
#include <vector>
#include <map>
#include <string>
#include <forward_list>

using namespace std;
using namespace Eigen;

namespace sns_ikl {

#define _USE_DOUBLE_  
/*! \def _USE_DOUBLE_  
 * use values with double precision.
 * This is used to have the possibility to switch (at compile time) between float and double precision in order to be faster or more accurate.
 * Inside the \b IKL we will use the follow definition:
 * - \b MatrixD  a dynamic sized Matrix of real numbers
 * - \b VectorD a column vector of real numbers with dynamic length 
 * - \b Scalar a real number
 */

#define _ONLY_WARNING_ON_ERROR
/*! \def _ONLY_WARNING_ON_ERROR
 * with this is possible to avoid to stop the code if some non critical errors appears.
 */

#ifdef _USE_DOUBLE_
/*! \typedef MatrixD
 *  dynamic sized matrix of \em float values
 *  (\em double if \b _USE_DOUBLE_ is defined)
 */
typedef Eigen::Matrix<double, Dynamic, Dynamic> MatrixD;
/*! \typedef VectorD
 *  dynamic sized column vector of \em float values
 *  (\em double if \b _USE_DOUBLE_ is defined)
 */
typedef Eigen::Matrix<double, Dynamic, 1> VectorD;
/*! \typedef Scalar
 *  a \em float value
 *  (\em double if \b _USE_DOUBLE_ is defined)
 */
typedef double Scalar;
#define EPS 0.15
#define LAMBDA_MAX 0.3
#define EPSQ 1e-10
#else
typedef Eigen::Matrix<float,Dynamic,Dynamic> MatrixD;
typedef Eigen::Matrix<float,Dynamic,1> VectorD;
typedef float Scalar;
#define EPS 0.15
#define LAMBDA_MAX 0.3
#define EPSQ 1e-15
#endif

#define INF 1e20;
#define SHAPE_MARGIN 0.98

// TODO: need to define each solver
enum inv_solvers {
  STD, SCALE, SNS, OSNS, FSNS, OSNSsm, FOSNS, QP, CHIAVERINI, RP, STD_MIN_ACC, ACC, RP_ST
};

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
typedef vector<Task> StackOfTasks;

class IKL {
  public:
    IKL(inv_solvers solver = STD);

    // set the method used to compute the inverse kinematic
    void setSolver(inv_solvers solver = STD);
    // set the method used to compute the inverse kinematic
    void setSolver(const string solver);
    void setJointsCapabilities(VectorD limit_low, VectorD limit_high, VectorD maxVelocity, VectorD maxAcceleration);
    void setNumberOfTasks(int ntasks, int dof = -1);
    void setNumberOfDOF(int dof);
    void setloopPeriod(double period) {
      loop_period = period;
    }
    void setScaleMargin(double scale = 0.9) {
      scaleMargin = scale;
    }
    void setNullSpaceDampingFactor(double factor = 0.95) {
      nullSpaceDampingFactor = factor;
    }

    Scalar getJointVelocity(VectorD *jointVelocity, StackOfTasks *sot, VectorD *jointConfiguration = NULL);  //inverse kinematic resolution
    vector<Scalar> getTasksScaleFactor() {
      return scaleFactors;
    }

  protected:
    Scalar getJointVelocity_STD(VectorD *jointVelocity, StackOfTasks *sot, VectorD *jointConfiguration = NULL);
    Scalar getJointVelocity_Chiaverini(VectorD *jointVelocity, StackOfTasks *sot, VectorD *jointConfiguration = NULL);
    Scalar getJointVelocity_Scale(VectorD *jointVelocity, StackOfTasks *sot, VectorD *jointConfiguration = NULL);
    Scalar getJointVelocity_QP(VectorD *jointVelocity, StackOfTasks *sot, VectorD *jointConfiguration = NULL);
    Scalar getJointVelocity_SNS(VectorD *jointVelocity, StackOfTasks *sot, VectorD *jointConfiguration = NULL);
    Scalar getJointVelocity_OSNS(VectorD *jointVelocity, StackOfTasks *sot, VectorD *jointConfiguration = NULL);
    Scalar getJointVelocity_FSNS(VectorD *jointVelocity, StackOfTasks *sot, VectorD *jointConfiguration = NULL);
    Scalar getJointVelocity_FOSNS(VectorD *jointVelocity, StackOfTasks *sot, VectorD *jointConfiguration = NULL);
    Scalar getJointVelocity_OSNSsm(VectorD *jointVelocity, StackOfTasks *sot, VectorD *jointConfiguration = NULL);
    Scalar getJointVelocity_RP(VectorD *jointVelocity, StackOfTasks *sot, VectorD *jointConfiguration = NULL);
    Scalar getJointVelocity_RP_ST(VectorD *jointVelocity, StackOfTasks *sot, VectorD *jointConfiguration = NULL);
    Scalar getJointVelocity_STD_MIN_ACC(VectorD *jointVelocity, StackOfTasks *sot, VectorD *jointConfiguration = NULL);
    Scalar getJointVelocity_ACC(VectorD *jointVelocity, StackOfTasks *sot, VectorD *jointConfiguration = NULL);

    // Used in getJointVelocity_RP_ST
    double activation(double position, double velocity, double min, double max, double bufferPos, double bufferVel);

    MatrixD invJ;

  private:
    inv_solvers inv_solver;

    // robot capabilities
    int n_dof;  //manipulator degree of freedom
    int n_tasks;  //number of tasks
    Scalar loop_period;  //needed to compute the bounds
    Scalar scaleMargin;  //used in the scale Margin approach to recover discontinuities
    Scalar nullSpaceDampingFactor;  //used to damp null space motion in RP_ST
    VectorD jointLimit_low;  //low joint limits
    VectorD jointLimit_high;  //high joint limit
    VectorD maxJointVelocity;  //maximum joint velocity
    VectorD maxJointAcceleration;  //maximum joint acceleration
    VectorD prev_JointVelocity;

    // Pointer to the getJointVelocity function used from outside.
    // The chosen method is associated at run time
    typedef Scalar (IKL::*getJointVelocityPtr)(VectorD *jointVelocity, StackOfTasks *sot, VectorD *jointConfiguration);
    getJointVelocityPtr getJointVelocityP;  //the pointer to the IK algorithm

    JacobiSVD<MatrixD> svd_Jt;
    JacobiSVD<MatrixD> svd_A;  //used to compute the pinv

    vector<MatrixD> W;  //selection matrices  [here to permit a warm start]
    VectorD dotQ;  // next solution (bar{\dotqv} in the paper)
    vector<VectorD> dotQopt;  // next solution (bar{\dotqv} in the paper)
    vector<Scalar> scaleFactors;
    MatrixD I;  //identity matrix
    Array<Scalar, Dynamic, 1> dotQmin;  //lower joint velocity bound
    Array<Scalar, Dynamic, 1> dotQmax;  //higher joint velocity bound

    // For the Fast versions of SNSs
    MatrixD B;  //update matrix
    VectorD dqw;  //saturations
    vector<VectorXi> S;  //the i-th element is zero if the i-th joint is not saturate, otherwise contains the position in B
    vector<VectorXi> sSat;  //contains the saturated joint in saturation order
    vector<int> nSat;  //number of saturated joint

    MatrixD J_prev;  //used to comute dot J
    VectorD dx_prev;  //used to compute ddot x

    // For the FastOpt version of the SNS
    vector<forward_list<int>> satList;
    VectorD lagrangeMu;
    VectorD lagrangeMu1;
    VectorD lagrangeMup2w;
    forward_list<int>::iterator it, prev_it;

    Scalar previousLambda;

    // Shape the joint velocity bound dotQmin and dotQmax
    void shapeJointVelocityBound(VectorD *actualJointConfiguration, double margin = SHAPE_MARGIN);
    // Perform the SNS for a single task
    Scalar SNSsingle(int priority, VectorD *higherPriorityJointVelocity, MatrixD *higherPriorityNull, MatrixD *jacobian,
        VectorD *task, VectorD *jointVelocity, MatrixD *nullSpaceProjector, double start_scale = 1);
    Scalar OSNSsingle(int priority, VectorD *higherPriorityJointVelocity, MatrixD *higherPriorityNull,
        MatrixD *jacobian, VectorD *task, VectorD *jointVelocity, MatrixD *nullSpaceProjector, double start_scale = 1);
    Scalar FSNSsingle(int priority, VectorD *higherPriorityJointVelocity, MatrixD *higherPriorityNull,
        MatrixD *jacobian, VectorD *task, VectorD *jointVelocity, MatrixD *nullSpaceProjector, double start_scale = 1);
    Scalar FOSNSsingle(int priority, VectorD *higherPriorityJointVelocity, MatrixD *higherPriorityNull,
        MatrixD *jacobian, VectorD *task, VectorD *jointVelocity, MatrixD *nullSpaceProjector, double start_scale = 1);
    Scalar QPsingle(int priority, VectorD *higherPriorityJointVelocity, MatrixD *higherPriorityAugmentedJacobian,
        MatrixD *jacobian, VectorD *task, VectorD *jointVelocity, MatrixD *augmentedJacobian, double start_scale = 1);

    // compute the pseudoinverse of A: it return 0 if A is (row) rank deficient
    bool pinv(MatrixD *A, MatrixD *invA, Scalar eps = EPS);
    // compute the pseudoinverse of A: it return 0 if A is (row) rank deficient
    // return also the null space projector P=(P-pinv(A)A)
    bool pinv_P(MatrixD *A, MatrixD *invA, MatrixD *P, Scalar eps = EPS);
    // compute the pseudoinverse of A: it return 0 if A is (row) rank deficient
    bool pinv_damped(MatrixD *A, MatrixD *invA, Scalar lambda_max = LAMBDA_MAX, Scalar eps = EPS);
    bool pinv_damped_P(MatrixD *A, MatrixD *invA, MatrixD *P, Scalar lambda_max = LAMBDA_MAX, Scalar eps = EPS);
    bool pinv_forBarP(MatrixD *W, MatrixD *P, MatrixD *inv);
    bool pinv_QR_Z(MatrixD *A, MatrixD *Z0, MatrixD *invA, MatrixD *Z, Scalar lambda_max = LAMBDA_MAX,
        Scalar eps = EPS);
    bool pinv_QR(MatrixD *A, MatrixD *invA, Scalar eps = EPS);
    void getTaskScalingFactor(Array<Scalar, Dynamic, 1> *a, Array<Scalar, Dynamic, 1> *b, MatrixD *W,
        Scalar *scalingFactor, int *mostCriticalJoint, Scalar maxScalingFactor = 1.0);
    void getTaskScalingFactor_fast(Array<Scalar, Dynamic, 1> *a, Array<Scalar, Dynamic, 1> *b, VectorXi *S,
        Scalar *scalingFactor, int *mostCriticalJoint, Scalar maxScalingFactor = 1.0);
    void getTaskScalingFactor_basic(Array<Scalar, Dynamic, 1> *a, Array<Scalar, Dynamic, 1> *b, Scalar *scalingFactor,
        int *mostCriticalJoint, Scalar maxScalingFactor = 1.0);

    map<string, inv_solvers> s_mapNameIKsolver;  //!<  mnemonic map jacobian name to jacobian value
    void initializeMap();  //!<  utility to initilize the mnemonic map

    bool isIdentity(MatrixD *A);
    bool isOptimal(int priority, VectorD *dotQ, MatrixD *tildeP, MatrixD *W, VectorD *dotQn, double eps = 1e-8);

};

}  //namespace

#endif
