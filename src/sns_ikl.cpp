/*! \file IKL.cpp
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

#include <sns_ikl/sns_ikl.hpp>

#include <iostream>
#include <math.h>
#include <Eigen/Dense>


//only for debug in ROS
#include <ros/ros.h>
#include <string>

using namespace std;
using namespace Eigen;
using namespace IKL_;

IKL::IKL(inv_solvers solver, inv_methods method, level in) {
  initializeMap();

  setSolver(solver);
  inv_method = method;//not used ?
  input_level = in; //not used ?

  saturate = true; //not used ?
  damp = false;//not used ?

  prev_JointVelocity = VectorD::Zero(1,1);
  J_prev = MatrixD::Zero(1,1);
  dx_prev = VectorD::Zero(1,1);
  nullSpaceDampingFactor = 0.95;

}

void IKL::setSolver(inv_solvers solver) {
  inv_solver=solver;  MatrixD J_prev;				//used to comute dot J

  switch(solver){
	  case STD:
      getJointVelocityP=&IKL::getJointVelocity_STD;
      ROS_INFO("solver set to: STD");
      break;
	  case CHIAVERINI:
      getJointVelocityP=&IKL::getJointVelocity_Chiaverini;
      ROS_INFO("solver set to: CHIAVERINI");
      break;
    case SCALE:
  		getJointVelocityP=&IKL::getJointVelocity_Scale;
  		ROS_INFO("solver set to: SCALE");
  		break;
  	case SNS:
  		getJointVelocityP=&IKL::getJointVelocity_SNS;
		  ROS_INFO("solver set to: SNS");
  		break;
  	case OSNS:
  		getJointVelocityP=&IKL::getJointVelocity_OSNS;
		  ROS_INFO("solver set to: OSNS");
		  break;
  	case OSNSsm:
  		getJointVelocityP=&IKL::getJointVelocity_OSNSsm;
		  ROS_INFO("solver set to: OSNSsm");
  		break;
  	case FSNS:
  		getJointVelocityP=&IKL::getJointVelocity_FSNS;
		  ROS_INFO("solver set to: FSNS");
  		break;
  	case FOSNS:
  		getJointVelocityP=&IKL::getJointVelocity_FOSNS;
		  ROS_INFO("solver set to: FOSNS");
  		break;
  	case QP:
  		getJointVelocityP=&IKL::getJointVelocity_QP;
		  ROS_INFO("solver set to: QP");
  		break;
  	case RP:
  		getJointVelocityP=&IKL::getJointVelocity_RP;
		  ROS_INFO("solver set to: RP");
  		break;
  	case STD_MIN_ACC:
  		getJointVelocityP=&IKL::getJointVelocity_STD_MIN_ACC;
		  ROS_INFO("solver set to: STD_MIN_ACC");
  		break;
  	case ACC:
  		getJointVelocityP=&IKL::getJointVelocity_ACC;
		  ROS_INFO("solver set to: ACC");
  		break;
  	case RP_ST:
  		getJointVelocityP=&IKL::getJointVelocity_RP_ST;
		  ROS_INFO("solver set to: RP_ST");
  		break;
  	default:
      cout <<"this IK solver has not be implemented (yet)" <<endl;
      break;
  }
}

void IKL::setSolver(const string solver) {
	setSolver(s_mapNameIKsolver[solver]);

}


void IKL::initializeMap() {
	s_mapNameIKsolver["STD"] = STD;
	s_mapNameIKsolver["CHIAVERINI"] = CHIAVERINI;
	s_mapNameIKsolver["SCALE"] = SCALE;
	s_mapNameIKsolver["SNS"] = SNS;
	s_mapNameIKsolver["OSNS"] = OSNS;
	s_mapNameIKsolver["FSNS"] = FSNS;
	s_mapNameIKsolver["OSNSsm"] = OSNSsm;
	s_mapNameIKsolver["FOSNS"] = FOSNS;
	s_mapNameIKsolver["QP"] = QP;
	s_mapNameIKsolver["RP"] = RP;
	s_mapNameIKsolver["STD_MIN_ACC"] = STD_MIN_ACC;
	s_mapNameIKsolver["ACC"] = ACC;
	s_mapNameIKsolver["RP_ST"] = RP_ST;

}

bool IKL::pinv(MatrixD *A, MatrixD *invA, Scalar eps) {

	//A (m x n) usually comes from a redundant task jacobian, therfore we consider m<n
	int m = A->rows() - 1;
	VectorD sigma;		//vector of singular values

  svd_A.compute(A->transpose(), ComputeThinU | ComputeThinV);
  sigma=svd_A.singularValues();
  if ( ((m > 0) && (sigma(m) > eps)) || ((m == 0) && (A->array().abs() > eps).any()) ) {
    for (int i = 0; i <= m; i++) {
      sigma(i) = 1.0 / sigma(i);
    }
    (*invA) = svd_A.matrixU() * sigma.asDiagonal() * svd_A.matrixV().transpose();
    return true;
  } else {
    return false;
  }
}

bool IKL::pinv_P(MatrixD *A, MatrixD *invA, MatrixD *P, Scalar eps) {

	//A (m x n) usually comes from a redundant task jacobian, therfore we consider m<n
	int m = A->rows() - 1;
	VectorD sigma;	 //vector of singular values

  svd_A.compute(A->transpose(), ComputeThinU | ComputeThinV);
  sigma = svd_A.singularValues();
  if ( ((m > 0)&& (sigma(m) > eps)) || ((m == 0) && (A->array().abs() > eps).any()) ) {
    for (int i = 0; i <= m; i++){
      sigma(i)= 1.0 / sigma(i);
    }
    (*invA) = svd_A.matrixU() * sigma.asDiagonal() * svd_A.matrixV().transpose();
    (*P) = ((*P) - svd_A.matrixU() * svd_A.matrixU().transpose()).eval();
    return true;
  }else{
    return false;
  }

}

bool IKL::pinv_damped(MatrixD *A, MatrixD *invA,Scalar lambda_max, Scalar eps) {

	//A (m x n) usually comes from a redundant task jacobian, therfore we consider m<n
	int m = A->rows() -1 ;
	VectorD sigma;		//vector of singular values
	Scalar lambda2;
	int r = 0;

  svd_A.compute(A->transpose(), ComputeThinU | ComputeThinV);
  sigma = svd_A.singularValues();
  if ( (( m > 0) && (sigma(m) > eps)) || ((m == 0) && (A->array().abs() > eps).any()) ) {
    for (int i = 0; i <= m; i++){
      sigma(i) = 1.0 / sigma(i);
    }
    (*invA) = svd_A.matrixU() * sigma.asDiagonal() * svd_A.matrixV().transpose();
    return true;
  }else{
    lambda2 = (1 - (sigma(m)/eps) * (sigma(m)/eps)) * lambda_max*lambda_max;
    for (int i = 0; i <= m; i++){
      if (sigma(i) > EPSQ) r++;
      sigma(i) = (sigma(i) / (sigma(i)*sigma(i)+lambda2));
    }
    //only U till the rank
    MatrixD subU = svd_A.matrixU().block(0,0,A->cols(),r);
    MatrixD subV = svd_A.matrixV().block(0,0,A->rows(),r);

    (*invA) = subU * sigma.asDiagonal() * subV.transpose();
    return false;
  }

}

bool IKL::pinv_damped_P(MatrixD *A, MatrixD *invA, MatrixD *P, Scalar lambda_max, Scalar eps) {

	//A (m x n) usually comes from a redundant task jacobian, therfore we consider m<n
	int m = A->rows() - 1;
	int r = 0; 				//rank
	VectorD sigma;		//vector of singular values
	Scalar lambda2;

  svd_A.compute(A->transpose(), ComputeThinU | ComputeThinV);
  sigma = svd_A.singularValues();
  if ( ((m > 0) && (sigma(m) > eps)) || ((m == 0) && (A->array().abs() > eps).any()) ) {
    for (int i = 0; i <= m; i++){
      sigma(i)= 1.0 / sigma(i);
    }
    (*invA) = svd_A.matrixU() * sigma.asDiagonal() * svd_A.matrixV().transpose();
    (*P) = ((*P) - svd_A.matrixU() * svd_A.matrixU().transpose()).eval();
    return true;
  }else{
    lambda2 = (1 - (sigma(m)/eps) * (sigma(m)/eps)) * lambda_max*lambda_max;
    for (int i = 0; i <= m; i++){
      if (sigma(i) > EPSQ) r++;
      sigma(i) = (sigma(i) / (sigma(i)*sigma(i)+lambda2));
    }

    //only U till the rank
    MatrixD subU = svd_A.matrixU().block(0,0,A->cols(),r);
    MatrixD subV = svd_A.matrixV().block(0,0,A->rows(),r);
    (*P) = ((*P) - subU*subU.transpose()).eval();

    (*invA) = subU * sigma.asDiagonal() * subV.transpose();
    return false;
  }

}

bool IKL::pinv_QR(MatrixD *A, MatrixD *invA, Scalar eps) {
	MatrixD At = A->transpose();
	HouseholderQR<MatrixD> qr = At.householderQr();
	int m = A->rows();
	int n = A->cols();

	MatrixD Rt = MatrixD::Zero(m,m);
	bool invertible;

	MatrixD hR = (MatrixD) qr.matrixQR();
	MatrixD Y = ((MatrixD) qr.householderQ()).leftCols(m);

	//take the useful part of R
	for (int i = 0; i < m; i++) {
		int j = 0;
		while ( j <= i) {
			Rt(i,j) = hR(j,i);
			j++;
		}
	}
	FullPivLU<MatrixD> invRt(Rt);

	invertible = abs(invRt.determinant()) > eps;

	if (invertible) {
		*invA = Y * invRt.inverse();
		return true;
	}else{
		return false;
	}

}


bool IKL::pinv_QR_Z(MatrixD *A, MatrixD *Z0, MatrixD *invA, MatrixD *Z,Scalar lambda_max, Scalar eps) {
	VectorD sigma;		//vector of singular values
	Scalar lambda2;

	MatrixD AZ0t = ((*A) * (*Z0)).transpose();
	HouseholderQR<MatrixD> qr = AZ0t.householderQr();
	MatrixD J_prev;				//used to comute dot J

	int m = A->rows();
	int p = Z0->cols();

	MatrixD Rt = MatrixD::Zero(m,m);
	bool invertible;
	MatrixD hR = (MatrixD) qr.matrixQR();
	MatrixD Y = ((MatrixD) qr.householderQ()).leftCols(m);

	//take the useful part of R
	for (int i=0; i < m; i++) {
		int j = 0;
		while (j <= i) {
			Rt(i,j) = hR(j,i);
			j++;
		}
	}

	FullPivLU<MatrixD> invRt(Rt);
	invertible = abs(invRt.determinant()) > eps;

	if (invertible) {
		*invA = (*Z0) * Y * invRt.inverse();
		*Z = (*Z0) * (((MatrixD) qr.householderQ()).rightCols(p-m));
		return true;
	}else{
		MatrixD R=MatrixD::Zero(m,m);
		//take the useful part of R
		for ( int i = 0; i < m; i++) {
			int j = i;
			while (j < m) {
				R(i,j) = hR(i,j);
				j++;
			}
		}

		//perform the SVD of R
    svd_A.compute(R, ComputeThinU | ComputeThinV);
    sigma = svd_A.singularValues();
    lambda2 = (1 - (sigma(m-1)/eps) * (sigma(m-1)/eps)) * lambda_max*lambda_max;
    for (int i = 0; i < m; i++) {
      sigma(i) = sigma(i) / (sigma(i)*sigma(i)+lambda2);
    }
    (*invA) = (*Z0) * Y * svd_A.matrixU() * sigma.asDiagonal() * svd_A.matrixV().transpose();

		*Z = (*Z0) * (((MatrixD) qr.householderQ()).rightCols(p-m));
		return false;
	}

}


bool IKL::pinv_forBarP(MatrixD *W,MatrixD *P, MatrixD *inv) {

	MatrixD barW;
	int rowsBarW = 0;

	MatrixD tmp;
	bool invertible;

	for (int i = 0; i < W->rows(); i++) {
		if ((*W)(i,i) > 0.99){ //equal to 1 (safer)
			rowsBarW++;
			barW = (MatrixD(rowsBarW,W->cols()) << barW, W->row(i) ).finished();
		}
	}

	tmp = barW * (*P) * barW.transpose();
	FullPivLU<MatrixD> inversePbar(tmp);

	invertible=inversePbar.isInvertible();

	if (invertible){
		(*inv) = (*P) * barW.transpose() * inversePbar.inverse() * barW;
		return true;
	}else {
		(*inv) = MatrixD::Zero(W->rows(), W->rows());
		return false;
	}
}



bool IKL::isIdentity(MatrixD *A) {

	bool isIdentity = true;
	int n = A->rows();
	int i = 0;
	do {
		isIdentity &= ((*A)(i,i) > 0.9); // equal to 1.0 (safer)
		i++;
	} while (isIdentity && i < n);

	return isIdentity;
}


Scalar IKL::getJointVelocity(VectorD *jointVelocity, StackOfTasks *sot, VectorD *jointConfiguration){
  //! \todo check for saturation if saturate is active 
  return (this->*getJointVelocityP)(jointVelocity,sot,jointConfiguration);
  
}


Scalar IKL::getJointVelocity_STD(VectorD *jointVelocity, StackOfTasks *sot, VectorD *jointConfiguration){
  
  int n_task=sot->size();
  int robotDOF=(*sot)[0].jacobian.cols();


  //P_0=I
  //dq_0=0
  MatrixD P = MatrixD::Identity(robotDOF,robotDOF);
  *jointVelocity = VectorD::Zero(robotDOF,1);
  VectorD higherPriorityJointVelocity;
  MatrixD higherPriorityNull;
  MatrixD tmp;
  
  for(int i_task=0; i_task<n_task; i_task++){ //consider all tasks
	// dq_k = dq_{k-1} - (J_k P_{k-1})^# (dx_k - J_k dq_{k-1})
	// P_k = P_{k-1} - (J_k P_{k-1})^# J_k P_{k-1}
	tmp=  (*sot)[i_task].jacobian * P;

    //if (i_task==0) pinv_damped_P(&tmp, &invJ, &P);
    //else pinv_P(&tmp, &invJ, &P);

    pinv_damped_P(&tmp, &invJ, &P);
    // svd_Jt.compute(((*sot)[i_task].jacobian * P).transpose(), ComputeThinU | ComputeThinV);
    //invJ=svd_Jt.matrixU()*svd_Jt.singularValues().asDiagonal().inverse()*svd_Jt.matrixV().transpose();
    //P=(P-svd_Jt.matrixU()*svd_Jt.matrixU().transpose()).eval();
      //P=(P-invJ*tmp).eval();


    *jointVelocity = ( (*jointVelocity) + invJ * ((*sot)[i_task].desired - (*sot)[i_task].jacobian*(*jointVelocity)));

  }
  return 1.0;
}

 Scalar IKL::getJointVelocity_STD_MIN_ACC(VectorD *jointVelocity, StackOfTasks *sot, VectorD *jointConfiguration){

    int n_task=sot->size();
    int robotDOF=(*sot)[0].jacobian.cols();

    if (prev_JointVelocity.rows()!=robotDOF){
    	prev_JointVelocity=VectorD::Zero(robotDOF,1);
    }

    //P_0=I
    //dq_0=0
    MatrixD P = MatrixD::Identity(robotDOF,robotDOF);
    *jointVelocity = VectorD::Zero(robotDOF,1);
    VectorD higherPriorityJointVelocity;
    MatrixD higherPriorityNull;
    MatrixD tmp;

    for(int i_task=0; i_task<n_task; i_task++){ //consider all tasks
  	// dq_k = dq_{k-1} - (J_k P_{k-1})^# (dx_k - J_k dq_{k-1})
  	// P_k = P_{k-1} - (J_k P_{k-1})^# J_k P_{k-1}
  	tmp=  (*sot)[i_task].jacobian * P;

      //if (i_task==0) pinv_damped_P(&tmp, &invJ, &P);
      //else pinv_P(&tmp, &invJ, &P);

      pinv_damped_P(&tmp, &invJ, &P);
      // svd_Jt.compute(((*sot)[i_task].jacobian * P).transpose(), ComputeThinU | ComputeThinV);
      //invJ=svd_Jt.matrixU()*svd_Jt.singularValues().asDiagonal().inverse()*svd_Jt.matrixV().transpose();
      //P=(P-svd_Jt.matrixU()*svd_Jt.matrixU().transpose()).eval();
        //P=(P-invJ*tmp).eval();


      *jointVelocity = ( (*jointVelocity) + invJ * ((*sot)[i_task].desired - (*sot)[i_task].jacobian*(*jointVelocity)));

    }

    double lambda=0.95;
    *jointVelocity +=  lambda*P*prev_JointVelocity;

    prev_JointVelocity=*jointVelocity;
  
  return 1.0;
}

 Scalar IKL::getJointVelocity_ACC(VectorD *jointVelocity, StackOfTasks *sot, VectorD *jointConfiguration){

     int n_task=sot->size();
     int robotDOF=(*sot)[0].jacobian.cols();
     int m=(*sot)[0].jacobian.rows();

     //THIS implementation considers a single task and is used only to test the algorithm

     if (prev_JointVelocity.rows()!=robotDOF){
    	 J_prev=(*sot)[0].jacobian;
    	 dx_prev=VectorD::Zero(m,1);
       	 prev_JointVelocity=VectorD::Zero(robotDOF,1);

     }

     //P_0=I
     //dq_0=0
     MatrixD P = MatrixD::Identity(robotDOF,robotDOF);
     *jointVelocity = VectorD::Zero(robotDOF,1);
     MatrixD jointAcceleration;
     MatrixD tmp;

     int i_task=0;
     tmp=  (*sot)[i_task].jacobian;

     MatrixD dotJ= (1/loop_period)*(tmp-J_prev);
     VectorD ddx=(1/loop_period)*((*sot)[i_task].desired - dx_prev);

     pinv_damped_P(&tmp, &invJ, &P);

     double lambda=0.95;
     double K=(1-lambda)/loop_period;


     jointAcceleration=  invJ * (ddx - dotJ*prev_JointVelocity) -K*P*prev_JointVelocity;
     *jointVelocity=  prev_JointVelocity + jointAcceleration*loop_period;


     prev_JointVelocity=*jointVelocity;
     J_prev= (*sot)[i_task].jacobian;
     dx_prev=(*sot)[i_task].desired;

   return 1.0;
 }


Scalar IKL::getJointVelocity_Chiaverini(VectorD *jointVelocity, StackOfTasks *sot, VectorD *jointConfiguration){

  int n_task=sot->size();
  int robotDOF=(*sot)[0].jacobian.cols();


  //P_0=I
  //dq_0=0
  MatrixD P = MatrixD::Identity(robotDOF,robotDOF);
  MatrixD Pnew = MatrixD::Identity(robotDOF,robotDOF);
  *jointVelocity = VectorD::Zero(robotDOF,1);
  VectorD higherPriorityJointVelocity;
  MatrixD higherPriorityNull;
  MatrixD tmp,invNonProjectedJ,tmp2;

  for(int i_task=0; i_task<n_task; i_task++){ //consider all tasks
	// dq_k = dq_{k-1} - (J_k P_{k-1})^# (dx_k - J_k dq_{k-1})
	// P_k = P_{k-1} - (J_k P_{k-1})^# J_k P_{k-1}
	tmp=  (*sot)[i_task].jacobian * P;

    //if (i_task==0) pinv_damped_P(&tmp, &invJ, &P);
    //else pinv_P(&tmp, &invJ, &P);

    pinv_damped_P(&tmp, &invJ, &Pnew);

    tmp2=  (*sot)[i_task].jacobian;
    pinv_damped(&tmp2, &invNonProjectedJ);

    // svd_Jt.compute(((*sot)[i_task].jacobian * P).transpose(), ComputeThinU | ComputeThinV);
    //invJ=svd_Jt.matrixU()*svd_Jt.singularValues().asDiagonal().inverse()*svd_Jt.matrixV().transpose();
    //P=(P-svd_Jt.matrixU()*svd_Jt.matrixU().transpose()).eval();
      //P=(P-invJ*tmp).eval();


    *jointVelocity = ( (*jointVelocity) + P*invNonProjectedJ * ((*sot)[i_task].desired - (*sot)[i_task].jacobian*(*jointVelocity)));

    P=Pnew;
  }

  return 1.0;
}


Scalar IKL::getJointVelocity_RP(VectorD *jointVelocity, StackOfTasks *sot, VectorD *jointConfiguration){

  int n_task=sot->size();
  int robotDOF=(*sot)[0].jacobian.cols();
  int n=robotDOF;
  int m,mtot;


  //P_0=I
  //dq_0=0
  MatrixD P = MatrixD::Identity(robotDOF,robotDOF);
  MatrixD I = MatrixD::Identity(robotDOF,robotDOF);
  *jointVelocity = VectorD::Zero(robotDOF,1);
  VectorD higherPriorityJointVelocity;
  MatrixD higherPriorityNull;
  MatrixD tmp,E,K,T,part,invK,part2,part3,J;
  int total_tasks_rows=0;

  for(int k_task=0; k_task<n_task; k_task++){ //consider all tasks
	  total_tasks_rows+=(*sot)[k_task].jacobian.rows();
  }


  MatrixD pJ_prev = MatrixD::Zero(robotDOF,total_tasks_rows+1);


  //fist the task with the lowest priority
  int i_task=n_task-1;

  J=  (*sot)[i_task].jacobian;
  if ((J-I).norm() < 1e-20){
	  *jointVelocity=((*sot)[i_task].desired );
	  i_task--;
  }else printf("?\n");

  m=(*sot)[i_task].jacobian.rows();
  J=  (*sot)[i_task].jacobian ;
  pinv_damped_P(&J, &invJ, &P);
  *jointVelocity = (*jointVelocity) + ( invJ * ((*sot)[i_task].desired - J*(*jointVelocity)));
  pJ_prev.block(0,0,n,m)=invJ;
  mtot=m;

  i_task-=1;
  for(; i_task>=0; --i_task){ //consider all tasks
	J=(*sot)[i_task].jacobian;

	m=J.rows();
	MatrixD Im = MatrixD::Identity(m,m);
	MatrixD pE = MatrixD::Zero(n,m);
	MatrixD pinvT = MatrixD::Identity(m,m);


	E=  J * P;
    pinv_damped_P(&E, &pE, &P);

	part3=pJ_prev.block(0,0,n,mtot);
    part=(part3*part3.transpose())*J.transpose();
    K=(Im + (Im-E*pE)*  J * part*(Im-E*pE));
    invK=K.inverse();
    //pinv_damped(&K, &invK);

    T=pE + (I-pE*J)* part * invK*(Im-E*pE);


    part2=pJ_prev.block(0,0,n,mtot);
    pJ_prev.block(0,0,n,mtot)=part2 - T*J*part2;
    pJ_prev.block(0,mtot+1,n,m)=T;
    mtot+=m;

    tmp=J*T;

    pinv_damped(&tmp, &invJ);
    //invJ=tmp.inverse();



    *jointVelocity = ( (*jointVelocity) + T*invJ * ((*sot)[i_task].desired - J*(*jointVelocity)));


  }

  return 1.0;
}


double IKL::activation(double position, double velocity, double min, double max, double bufferPos,double bufferVel){

	double cx,hx,vx,cx1,hx1,vx1,h;

	cx=(position-(min+bufferPos))/(-bufferPos);
    cx=(cx<0)?0:(cx>1)?1:cx;
    hx=6*cx*cx*cx*cx*cx-15*cx*cx*cx*cx+10*cx*cx*cx;
    vx=(velocity-bufferVel)/(-bufferVel);
    vx=(vx<0)?0:(vx>1)?1:vx;
    vx=6*vx*vx*vx*vx*vx-15*vx*vx*vx*vx+10*vx*vx*vx;

	cx1=(position-(max-bufferPos))/(bufferPos);
    cx1=(cx1<0)?0:(cx1>1)?1:cx1;
    hx1=6*cx1*cx1*cx1*cx1*cx1-15*cx1*cx1*cx1*cx1+10*cx1*cx1*cx1;
    vx1=(velocity+bufferVel)/(bufferVel);
    vx1=(vx1<0)?0:(vx1>1)?1:vx1;
    vx1=6*vx1*vx1*vx1*vx1*vx1-15*vx1*vx1*vx1*vx1+10*vx1*vx1*vx1;

    h=hx*vx + hx1*vx1;

   // ROS_INFO("pos %.2f, vel %.2f, min %.2f,max %.2f, hmin %.2f, hmax %.2f, act %d",position,velocity,min,max,hx*vx,hx1*vx1,(h<1e-20)?0:1);

    return h;

}


Scalar IKL::getJointVelocity_RP_ST(VectorD *jointVelocity, StackOfTasks *sot, VectorD *jointConfiguration){

  int n_task=sot->size();
  int robotDOF=(*sot)[0].jacobian.cols();
  int n=robotDOF;
  int m,mtot=0;
  double h;
  int activeConstraints=0;


/*
  if (prev_JointVelocity.rows()!=robotDOF){
  	prev_JointVelocity=VectorD::Zero(robotDOF,1);
  }
*/
  prev_JointVelocity=(*jointVelocity);

  //P_0=I
  //dq_0=0
  MatrixD P = MatrixD::Identity(robotDOF,robotDOF);
  MatrixD P1 = MatrixD::Identity(robotDOF,robotDOF);
  MatrixD I = MatrixD::Identity(robotDOF,robotDOF);
  *jointVelocity = VectorD::Zero(robotDOF,1);
  VectorD higherPriorityJointVelocity;
  MatrixD higherPriorityNull;
  MatrixD tmp,E,K,T,part,invK,part2,part3,J;
  double lambda;
/*
  //remove smooth transition
  for(int k_task=0; k_task<n_task; k_task++){ //consider all tasks
	  if ((*sot)[k_task].h>1e-20){
		(*sot)[k_task].h=1;
	  }else{
		(*sot)[k_task].h=0;
	  }
  }
*/
/*
  bool inTransition=false;
  for(int k_task=0; k_task<n_task; k_task++){ //consider all tasks
	  if (((*sot)[k_task].h>0)&&((*sot)[k_task].h<1)){
		  inTransition=true;
	  }
  }
  double delta_lambda=0.0005;
  if (inTransition) {
	  lambda=previousLambda + delta_lambda;
  }else{
	  lambda=previousLambda - delta_lambda;
  }
  lambda=(lambda>1)?1:lambda;
  lambda=(lambda<nullSpaceDampingFactor)?nullSpaceDampingFactor:lambda;
 */
  lambda=0.95;//nullSpaceDampingFactor;
  previousLambda=lambda;
  double lambda2=1;//0.2
  *jointVelocity =  lambda*prev_JointVelocity;

  int total_tasks_rows=0;

   for(int k_task=0; k_task<n_task; k_task++){ //consider all tasks
 	  total_tasks_rows+=(*sot)[k_task].jacobian.rows();
   }


   MatrixD pJ_prev = MatrixD::Zero(robotDOF,total_tasks_rows+1);


  //fist the task with the lowest priority
  int i_task=n_task-1;
 // while (((*sot)[i_task].h<1e-20) && (( (*sot)[i_task].h0<1e-20 )||((*sot)[i_task].h0>0.9999999999999999999) ) && (i_task>=0)){
  while (((*sot)[i_task].h<1e-20) && (i_task>=0)){
	  i_task--;
  }
  if (i_task<0) return 0;

  J=  (*sot)[i_task].jacobian ;
  //PG task
  if ((J-I).norm() < 1e-20){
//	  *jointVelocity =  (*jointVelocity) + lambda2*((*sot)[i_task].h*(*sot)[i_task].desired + (*sot)[i_task].h0*((*sot)[i_task].q0));
	  *jointVelocity =  (*jointVelocity) + (1-lambda)*((*sot)[i_task].h*(*sot)[i_task].desired);
	  i_task--;
	  while (((*sot)[i_task].h<1e-20) && (i_task>=0)){
		  i_task--;
	  }
	  if (i_task<0) return 0.1;
  }

  m=(*sot)[i_task].jacobian.rows();
  J=  (*sot)[i_task].jacobian ;
  pinv_damped_P(&J, &invJ, &P);

  if (!(*sot)[i_task].is_constraints){
	  h=(*sot)[i_task].h;
  }else{
	  double velocity=(J*(*jointVelocity)).norm();
	  h=activation((*sot)[i_task].value, velocity, (*sot)[i_task].min, (*sot)[i_task].max, 0.01,0.03);
	  if (h>1e-20) activeConstraints++;
	  while ((h<1e-20) && (i_task>=0)){
		  i_task--;
	  }
	  if (i_task<0) return 0.1;

  }


  *jointVelocity =  (*jointVelocity) + h*invJ * ( ((*sot)[i_task].desired - J*(*jointVelocity)) )  ;
  pJ_prev.block(0,0,n,m)=invJ;
  mtot+=m;

  for(i_task-=1; i_task>=0; i_task--){ //consider all tasks

	  //	if ( ((*sot)[i_task].h<1e-20) && (( (*sot)[i_task].h0<1e-20 )||((*sot)[i_task].h0>0.9999999999999999999) ) ) continue;
		if (((*sot)[i_task].h<1e-20)) continue;

	J=(*sot)[i_task].jacobian;

	  if (!(*sot)[i_task].is_constraints){
		  h=(*sot)[i_task].h;
	  }else{
		  MatrixXd tmp= J*(*jointVelocity); // here I suppose directly that the costraints has dimension 1
		  double velocity=tmp(0);
		  h=activation((*sot)[i_task].value, velocity, (*sot)[i_task].min, (*sot)[i_task].max, 0.01,0.03);
		  if (h>1e-20) activeConstraints++;
//		  if (i_task==0)
//			  ROS_INFO("pos %.2f, vel %.2f, min %.2f,max %.2f, h %.2f, act %d",(*sot)[i_task].value,velocity,(*sot)[i_task].min,(*sot)[i_task].max,h,(h<1e-20)?0:1);

	  }
	if (h<1e-20) continue;

	m=J.rows();
	MatrixD Im = MatrixD::Identity(m,m);
	MatrixD pE = MatrixD::Identity(n,m);
	MatrixD pinvT = MatrixD::Identity(m,m);


	E=  J * P;
	//if (E.norm()>1e-12){
	pinv_damped_P(&E, &pE, &P);


		part3=pJ_prev.block(0,0,n,mtot);
		part=(part3*part3.transpose())*J.transpose();
		K=(Im + (Im-E*pE)*  J * part*(Im-E*pE));
		invK=K.inverse();
		//pinv_damped(&K, &invK);
		T=pE + (I-pE*J)* part * invK*(Im-E*pE);

		part2=pJ_prev.block(0,0,n,mtot);
		pJ_prev.block(0,0,n,mtot)=part2 - T*J*part2;
		pJ_prev.block(0,mtot+1,n,m)=T;
		mtot+=m;

		tmp=J*T;
		//if (!pinv_damped(&tmp, &invJ)) ROS_WARN("damped JT task %d",i_task);
		invJ=tmp.inverse();

		//double h=(*sot)[i_task].h;
		//printf("t %d h %f\n",i_task,h);

	   *jointVelocity =  (*jointVelocity) + h*T*invJ * (((*sot)[i_task].desired - J*(*jointVelocity)) )  ;
//----------------------TEST
	  // pJ_prev*=1e-2;
//-------------------------

	   //   *jointVelocity =  (*jointVelocity) + T*invJ * ( h*(*sot)[i_task].desired + (1-h)*J*(*sot)[i_task].q0 - J*(*jointVelocity)  );
	   // if ((*sot)[i_task].h>1e-10) P=P1;
/*	}else{
		//just to try

		T=pJ_prev.block(0,0,n,m);
		tmp=J*T;
		invJ=tmp.inverse();
	   *jointVelocity =  (*jointVelocity) + T*invJ * ( (*sot)[i_task].h*((*sot)[i_task].desired - J*(*jointVelocity)) + (*sot)[i_task].h0* J *((*sot)[i_task].q0 - (*jointVelocity))  )  ;

		pinv_damped_P(&J, &invJ, &P1);
		*jointVelocity =  (*jointVelocity) + invJ * ( (*sot)[i_task].h*((*sot)[i_task].desired - J*(*jointVelocity)) + (*sot)[i_task].h0* J *((*sot)[i_task].q0 - (*jointVelocity))  )  ;


	}
*/
	   /*
	   if (i_task==1){
		   MatrixXd tmpa= (*sot)[i_task].jacobian*(*jointVelocity); // here I suppose directly that the costraints has dimension 1
		   vel0=tmpa(0);
	   }*/
  }
  /*if (printPlease){
	  MatrixXd tmpa= (*sot)[1].jacobian*(*jointVelocity); // here I suppose directly that the costraints has dimension 1
	  double velocity1=tmpa(0);
	  //ROS_INFO("task 1 vel before %f, after %f",velocity,velocity1);
	  ROS_INFO(" vel %f ",velocity1);
	  //if ((H_X_PER_PROVA>0.9999)) ROS_WARN("error constraints vel was %f to %f h %f",vel0,velocity1,H_X_PER_PROVA);
  }*/
  return activeConstraints;
}


Scalar IKL::getJointVelocity_Scale(VectorD *jointVelocity, StackOfTasks *sot, VectorD *jointConfiguration){

  int n_task=sot->size();
  int robotDOF=(*sot)[0].jacobian.cols();
  Array<Scalar,Dynamic,1> a,b;					   // used to compute the task scaling factor
  double scalingFactor;
  int mostCriticalJoint; // not used here

  //P_0=I
  //dq_0=0
  MatrixD P = MatrixD::Identity(robotDOF,robotDOF);
  *jointVelocity = VectorD::Zero(robotDOF,1);
  VectorD jointVelocityA,jointVelocityB;
  VectorD higherPriorityJointVelocity;
  MatrixD higherPriorityNull;
  MatrixD tmp;

  for(int i=0; i<n_dof; i++){
  		dotQmin(i)=-maxJointVelocity(i);
  		dotQmax(i)=maxJointVelocity(i);
  	}


  double s[2];

  for(int i_task=0; i_task<n_task; i_task++){ //consider all tasks
	// dq_k = dq_{k-1} - (J_k P_{k-1})^# (dx_k - J_k dq_{k-1})
	// P_k = P_{k-1} - (J_k P_{k-1})^# J_k P_{k-1}
/*
	for(int j=0; j<n_dof; j++){
		dotQmin(j)=- maxJointVelocity(j) - (*jointVelocity)(j);
		dotQmax(j)= maxJointVelocity(j) - (*jointVelocity)(j);
	}
*/
	tmp=  (*sot)[i_task].jacobian * P;


    pinv_damped_P(&tmp, &invJ, &P);


    jointVelocityA= invJ *(*sot)[i_task].desired;
    jointVelocityB= ((*jointVelocity) + invJ * (- (*sot)[i_task].jacobian*(*jointVelocity)));

    //*jointVelocity = ( (*jointVelocity) + invJ * ((*sot)[i_task].desired - (*sot)[i_task].jacobian*(*jointVelocity)));
	a=jointVelocityA.array();
	b=jointVelocityB.array();

	getTaskScalingFactor_basic(&a, &b, &scalingFactor, &mostCriticalJoint);
	if (scalingFactor>1) scalingFactor=1;

	if (scalingFactor>0){
		*jointVelocity =scalingFactor*jointVelocityA + jointVelocityB;
	}
	s[i_task]=scalingFactor;

  }
  //ROS_INFO("scale %f  %f",s[0],s[1]);

  return 1.0;
}


void IKL::setNumberOfDOF(int dof){
	n_dof=dof;
	I=MatrixD::Identity(n_dof,n_dof);
    dotQ=VectorD::Zero(n_dof);
}

void IKL::setNumberOfTasks(int ntasks,int dof){

	n_tasks=ntasks;
	if (dof!=-1){
		n_dof=dof;
		I=MatrixD::Identity(n_dof,n_dof);
	}

	Scalar scale=1.0;
	VectorD dq=VectorD::Zero(n_dof);

	for( int i=0; i<n_tasks; i++){
		W.push_back(I);
		scaleFactors.push_back(scale);
		dotQopt.push_back(dq);
	}

	//for the Fast version
	MatrixD Z=MatrixD::Zero(n_dof,n_dof);
	VectorD zv=VectorD::Zero(n_dof);
	VectorXi zvi=VectorXi::Zero(n_dof);
	B=Z;
	dqw=zv;

	for( int i=0; i<n_tasks; i++){
	//	Jw.push_back(Z);
		S.push_back(zvi);
		sSat.push_back(zvi);
		nSat.push_back(0);
	}
	satList.resize(n_tasks);
	lagrangeMu=zv;
	lagrangeMu1=zv;
	lagrangeMup2w=zv;
}

void IKL::setJointsCapabilities(VectorXd limit_low, VectorXd limit_high, VectorXd maxVelocity, VectorXd maxAcceleration){
	jointLimit_low=limit_low;
	jointLimit_high=limit_high;
	maxJointVelocity=maxVelocity;
	maxJointAcceleration=maxAcceleration;

	dotQmin=-maxJointVelocity.array();
	dotQmax=maxJointVelocity.array();
}

void IKL::shapeJointVelocityBound(VectorD *actualJointConfiguration, double margin){

	//it could be written using the Eigen::Array potentiality
	double step,max,stop;

	for (int i=0; i<n_dof; i++) {
		//for the minimum bound
		step=(jointLimit_low(i)-(*actualJointConfiguration)(i))/loop_period;
		max=-maxJointVelocity(i);
		stop=- sqrt(2*maxJointAcceleration(i)*((*actualJointConfiguration)(i)-jointLimit_low(i)));
		dotQmin(i)=(step>max)?((step>stop)?step: stop):((max>stop)?max:stop); //take the maximum

		//for the maximum bound
		step=(jointLimit_high(i)-(*actualJointConfiguration)(i))/loop_period;
		max=maxJointVelocity(i);
		stop=sqrt(2*maxJointAcceleration(i)*(jointLimit_high(i)-(*actualJointConfiguration)(i)));
		dotQmax(i)=(step<max)?((step<stop)?step: stop):((max<stop)?max:stop); //take the minimum
	}

	dotQmin*=margin;
	dotQmax*=margin;
}

void IKL::getTaskScalingFactor(Array<Scalar,Dynamic,1> *a, Array<Scalar,Dynamic,1> *b,MatrixD *W,Scalar *scalingFactor,int *mostCriticalJoint,Scalar maxScalingFactor){

	Array<Scalar,Dynamic,1> Smin,Smax;
	Scalar temp, smax, smin;
	Scalar inf = INF;
	int col;

	Smin = (dotQmin - (*b))/(*a);
	Smax = (dotQmax - (*b))/(*a);

	for (int i=0; i < a->rows(); i++) {
		//switch
		if (Smin(i) > Smax(i)){
			temp=Smin(i);
			Smin(i)=Smax(i);
			Smax(i)=temp;
		}
		//remove saturated
		if (((*W)(i,i)<0.2)|| ((*a)(i)==0)) { // if it is not 1 (safer)
			Smin(i)=-INF;
			Smax(i)=INF;
		}
	}

	smax=Smax.minCoeff(mostCriticalJoint,&col);
	smin=Smin.maxCoeff();

	if ( (smin>smax) || (smax<0.0) || (smin>1.0)|| (smax== inf )){
		(*scalingFactor)=-1.0; // the task is not feasible
	}else{
		//if (smax>maxScalingFactor) (*scalingFactor)=maxScalingFactor;
		//else
			(*scalingFactor)=smax;
	}

}

void IKL::getTaskScalingFactor_basic(Array<Scalar,Dynamic,1> *a, Array<Scalar,Dynamic,1> *b,Scalar *scalingFactor,int *mostCriticalJoint,Scalar maxScalingFactor){
	Array<Scalar,Dynamic,1> Smin,Smax;
	Scalar temp,smax,smin;
	Scalar inf=INF;
	int col;

	Smin=(dotQmin - (*b))/(*a);
	Smax=(dotQmax - (*b))/(*a);

	for (int i=0; i < a->rows(); i++){
		//switch
		if (Smin(i) > Smax(i)){
			temp=Smin(i);
			Smin(i)=Smax(i);
			Smax(i)=temp;
		}
		//remove saturated
		if (((*a)(i)==0)) { // if it is not 0 (safer)
			Smin(i)=-inf;
			Smax(i)=inf;

		}

	}

	smax=Smax.minCoeff(mostCriticalJoint,&col);
	smin=Smin.maxCoeff();

/*
	stringstream log;
	log<<"/na "<<a->transpose()<<endl;
	log<<"/nb "<<b->transpose()<<endl;
	log<<"/nSmax "<<Smax.transpose()<<endl;
	log<<"/nSmin "<<Smin.transpose()<<endl;
	log<<"/nsmax - smin  "<<smax<<" - "<<smin<<endl;
	ROS_INFO("%s",log.str().c_str());
*/


	if ( (smin>smax) || (smax<0.0) || (smin>1.0) || (smax== inf )) {
		(*scalingFactor)=-1.0; // the task is not feasible
	}else{
		//if (smax>maxScalingFactor) (*scalingFactor)=maxScalingFactor;
		//else
			(*scalingFactor)=smax;
	}

}


void IKL::getTaskScalingFactor_fast(Array<Scalar,Dynamic,1> *a, Array<Scalar,Dynamic,1> *b,VectorXi *S,Scalar *scalingFactor,int *mostCriticalJoint,Scalar maxScalingFactor){

	Array<Scalar,Dynamic,1> Smin,Smax;
	Scalar temp,smax,smin;
	Scalar inf=INF;
	int col;

	Smin=(dotQmin - (*b))/(*a);
	Smax=(dotQmax - (*b))/(*a);

	for (int i=0; i < a->rows(); i++){
		//switch
		if (Smin(i) > Smax(i)){
			temp=Smin(i);
			Smin(i)=Smax(i);
			Smax(i)=temp;
		}
		//remove saturated
		if ( ((*S)(i)) || ((*a)(i)==0)) { // if it is not 0 (safer)
			Smin(i)=-inf;
			Smax(i)=inf;

		}
	}


	smax=Smax.minCoeff(mostCriticalJoint,&col);
	smin=Smin.maxCoeff();
	if ( (smin>smax) || (smax<0.0) || (smin>1.0) || (smax==inf) ){
		(*scalingFactor)=-1.0; // the task is not feasible
	}else{
		//if (smax>maxScalingFactor) (*scalingFactor)=maxScalingFactor;
		//else
			(*scalingFactor)=smax;
	}

}

bool IKL::isOptimal(int priority, VectorD *dotQ, MatrixD *tildeP, MatrixD *W ,VectorD *dotQn ,double eps){

	VectorD barMu;
	bool isOptimal=true;

	barMu=tildeP->transpose()*(*dotQ);

	for (int i=0; i<n_dof; i++){
		if ((*W)(i,i) < 0.1 ) { //equal to 0.0 (safer)

			if ( abs((*dotQ)(i) - dotQmax(i)) < eps ) barMu(i)=-barMu(i);

			if ( barMu(i) < 0.0 ){
				(*W)(i,i)=1.0;
				(*dotQn)(i)=0.0;
				isOptimal=false;
			}
		}

	}
	return isOptimal;
}

Scalar IKL::SNSsingle(int priority,VectorD *higherPriorityJointVelocity,MatrixD *higherPriorityNull, MatrixD *jacobian, VectorD *task, VectorD *jointVelocity,MatrixD *nullSpaceProjector,double start_scale){

	//INITIALIZATION
	VectorD tildeDotQ;
	MatrixD projectorSaturated; 			//(((I-W_k)*P_{k-1})^#
	MatrixD JPinverse;					//(J_k P_{k-1})^#
	MatrixD temp;
	bool isW_identity;
	MatrixD barP=*higherPriorityNull;  // remove the pointer arguments advantage... but I need to modify it
	Array<Scalar,Dynamic,1> a,b;					   // used to compute the task scaling factor
	bool limit_excedeed;
	bool singularTask=false;
	bool reachedSingularity=false;
	Scalar scalingFactor=1.0;
	int mostCriticalJoint;
	//best solution
	Scalar bestScale=-1.0;
	VectorD bestTildeDotQ;
	MatrixD bestInvJP;
	VectorD bestDotQn;
    VectorD dotQn;		//saturate velocity in the null space

	int k;

//INIT
	W[priority]=I;
	isW_identity=true;
	dotQn=VectorD::Zero(n_dof);


	//SNS
	int count=0;
	do{
		count++;
		//ROS_INFO("%d",count);
		if (count>2*n_dof){
#ifndef _ONLY_WARNING_ON_ERROR
			ROS_ERROR("Infinite loop on SNS for task (%d)",priority);

			ROS_INFO("p:%d  scale:%f  mc:%d  sing:%d",priority,scalingFactor,mostCriticalJoint,(int)reachedSingularity);
			//##############################
			string s;
			stringstream buffer;
			streambuf * old = std::cout.rdbuf(buffer.rdbuf());
			cout << "W\n" << W[priority] << endl << "P(k-1)\n" <<(*higherPriorityNull)<< endl<<"Psat\n" << projectorSaturated <<endl << "dotQ\n" << dotQ <<endl << "JPinverse\n" << JPinverse <<endl;
			s = buffer.str();
			ROS_INFO("\n p %d \n %s",priority,s.c_str());
			//#############################

			exit(1);
#else
			ROS_WARN("Infinite loop on SNS for task (%d)",priority);
			ROS_INFO("p:%d  scale:%f  mc:%d  sing:%d",priority,scalingFactor,mostCriticalJoint,(int)reachedSingularity);
			// the task is not executed
			*jointVelocity= *higherPriorityJointVelocity;
			*nullSpaceProjector=*higherPriorityNull;
			limit_excedeed=false;
			continue;
#endif
		}
		limit_excedeed=false;

		// remember that in the SNS W==I always and only on the first loop
		if (isW_identity){
			tildeDotQ=(*higherPriorityJointVelocity);//VectorD::Zero(n_dof);
			//compute (J P)^#
			temp=(*jacobian)*(*higherPriorityNull);
			singularTask=!pinv_damped_P(&temp, &JPinverse, nullSpaceProjector);
		}else{
			//JPinverse is already computed
			tildeDotQ=(*higherPriorityJointVelocity) + projectorSaturated*dotQn;
		}
		dotQ= tildeDotQ + JPinverse *( (*task)  - (*jacobian)*tildeDotQ);


		a=(JPinverse * (*task)).array();
		b=dotQ.array() - a;
		getTaskScalingFactor(&a, &b,&W[priority], &scalingFactor, &mostCriticalJoint);

		if (scalingFactor >= 1.0){
			// task accomplished

		}else {

			limit_excedeed=true;
	/*		if (scalingFactor < 0.0){
				//the task is not feasible
				if (isW_identity) singularTask=true;
				else reachedSingularity=true;
			}
			*/
			if(singularTask){
				// the task is singular so return a scaled damped solution (no SNS possible)
				//ROS_WARN("task %d is singular, scaling factor: %f",priority,scalingFactor);
				if (scalingFactor >= 0.0){
					(*jointVelocity)= tildeDotQ + JPinverse *( scalingFactor*(*task)  - (*jacobian)*tildeDotQ);
				}else{
					// the task is not executed
					//ROS_INFO("task not executed: J sing");
					//W[priority]=I;
					//dotQn=VectorD::Zero(n_dof);
					*jointVelocity= *higherPriorityJointVelocity;
					*nullSpaceProjector=*higherPriorityNull;
				}


				return scalingFactor;
			}

			if ((scalingFactor>bestScale)){
				//save best solution so far
				bestScale=scalingFactor;
				bestTildeDotQ=tildeDotQ;
				bestInvJP=JPinverse;
				bestDotQn=dotQn;
			}


				//saturate the most critical join
			W[priority](mostCriticalJoint,mostCriticalJoint)=0.0;
			isW_identity=false;
			if (dotQ(mostCriticalJoint)>dotQmax(mostCriticalJoint)){
				dotQn(mostCriticalJoint)=dotQmax(mostCriticalJoint)-(*higherPriorityJointVelocity)(mostCriticalJoint);
			}else{
				dotQn(mostCriticalJoint)=dotQmin(mostCriticalJoint)-(*higherPriorityJointVelocity)(mostCriticalJoint);
			}


			if (priority==0){ //for the primary task higherPriorityNull==I
				barP=W[0];
				projectorSaturated=(I-W[0]);
			}else{
				temp=(I-W[priority]);//*(*higherPriorityNull);
				reachedSingularity|=!pinv_forBarP(&temp,higherPriorityNull,&projectorSaturated);
				//temp=(I-W[priority])*(*higherPriorityNull);
				//reachedSingularity|=!Spinv(&temp,&projectorSaturated);

				barP=(I-projectorSaturated)*(*higherPriorityNull);
			}

			temp=(*jacobian)*barP;

			reachedSingularity|=!pinv(&temp, &JPinverse);



			if (reachedSingularity){
				if (bestScale >= 0.0){

					//ROS_INFO("best solution %f",bestScale);
					dotQn=bestDotQn;
					dotQ=bestTildeDotQ + bestInvJP *( bestScale*(*task)  - (*jacobian)*bestTildeDotQ);
					//use the best solution found... no further saturation possible
					(*jointVelocity)= dotQ;
				}else{
					// the task is not executed
					//ROS_WARN("task not executed: reached sing");
					for(int i=0; i< 1e8; i++){};
					*jointVelocity= *higherPriorityJointVelocity;
					*nullSpaceProjector=*higherPriorityNull;
				}

				return bestScale;
			}

		}

	}while(limit_excedeed);

	(*jointVelocity)=dotQ;
	return 1.0;
}

Scalar IKL::OSNSsingle(int priority,VectorD *higherPriorityJointVelocity,MatrixD *higherPriorityNull, MatrixD *jacobian, VectorD *task, VectorD *jointVelocity,MatrixD *nullSpaceProjector,double start_scale){

	//INITIALIZATION
	//VectorD tildeDotQ;
	MatrixD projectorSaturated; 			//(((I-W_k)*P_{k-1})^#
	MatrixD JPinverse;					//(J_k P_{k-1})^#
	MatrixD temp;
	bool isW_identity;
	MatrixD barP=*higherPriorityNull;  // remove the pointer arguments advantage... but I need to modify it
	Array<Scalar,Dynamic,1> a,b;					   // used to compute the task scaling factor
	bool limit_excedeed;
	bool singularTask=false;
	bool reachedSingularity=false;
	Scalar scalingFactor=1.0;
	int mostCriticalJoint;
	bool exceed;
	bool push_out;
	//best solution
	Scalar bestScale=-0.1;
	MatrixD bestW; //(only in OSNS)
	//MatrixD bestPS;
	//VectorD bestTildeDotQ;
	MatrixD bestInvJP;
	VectorD bestDotQn;
	MatrixD bestTildeP;

    VectorD dotQn;		//saturate velocity in the null space
	VectorD dotQs;
	MatrixD tildeP;       // used in the  OSNS
	int k;

	//these two are needed to consider also W=I in case of a non feasible task
	//bool nonSaturatedTried=false;
	bool invJPcomputed=false;

	nSat[priority]=0;
	//bool optimal=true;

//	ROS_INFO("NEW");
//Compute the solution with W=I it is needed anyway to obtain nullSpaceProjector-----------------
	//tildeDotQ=(*higherPriorityJointVelocity);//VectorD::Zero(n_dof);
	//compute (J P)^#
	temp=(*jacobian)*(*higherPriorityNull);
	singularTask=!pinv_damped_P(&temp, &JPinverse, nullSpaceProjector);
	//nonSaturatedTried=true;
	tildeP=MatrixD::Zero(n_dof,n_dof);
	dotQs=  (*higherPriorityJointVelocity) + JPinverse *( (*task)  - (*jacobian)*(*higherPriorityJointVelocity));
	a=(JPinverse * (*task)).array();
	b=dotQs.array() - a;
	getTaskScalingFactor(&a, &b,&I, &scalingFactor, &mostCriticalJoint);

	//double scalingI=scalingFactor;
	if (scalingFactor >= 1.0){
		// then is clearly the optimum since all joints velocity are computed with the pseudoinverse
//ROS_INFO("solved I");
		(*jointVelocity)=dotQs;
		W[priority]=I;
		dotQopt[priority]=dotQs;
		return scalingFactor;
	}

	if (singularTask){
		// the task is singular so return a scaled damped solution (no SNS possible)
		//ROS_ERROR("sing");
		if (scalingFactor >= 0.0){
			W[priority]=I;
			(*jointVelocity)= (*higherPriorityJointVelocity) + scalingFactor*JPinverse *(*task)  + tildeP*(*higherPriorityJointVelocity);
			dotQopt[priority]=(*jointVelocity);
		}else{
			// the task is not executed
			W[priority]=I;
			//dotQn=VectorD::Zero(n_dof);
			*jointVelocity= *higherPriorityJointVelocity;
			dotQopt[priority]=(*jointVelocity);
			*nullSpaceProjector=*higherPriorityNull;
		}
		return scalingFactor;
	}

	if (scalingFactor>bestScale){
//		ROS_INFO("best scale I %f",scalingFactor);

		//save best solution so far
		bestScale=scalingFactor;
		//bestTildeDotQ=tildeDotQ;
		bestInvJP=JPinverse;
		bestW=I;
		bestTildeP=tildeP;
		//bestPS=projectorSaturated;
		bestDotQn=VectorD::Zero(n_dof);
	}

	W[priority]=I; //test: do not use the warm start
//----------------------------------------------------------------------- END W=I

	//INIT
	dotQn=VectorD::Zero(n_dof);
	if (isIdentity(&(W[priority]))) {
		isW_identity=true;
		dotQopt[priority]=dotQs; // use the one computed above
	}
	else {
		isW_identity=false;
		for (int i=0; i<n_dof; i++){
			if ((W[priority])(i,i) <0.1){ // equal to 0.0 (safer)
				nSat[priority]++;
				//I'm not considering a different dotQ for each task... this could be a little improvement
				if (dotQopt[priority](i)>=0.0){
					dotQn(i)=dotQmax(i)-(*higherPriorityJointVelocity)(i);
				}else{
					dotQn(i)=dotQmin(i)-(*higherPriorityJointVelocity)(i);
				}
			}else dotQn(i)=0.0;
		}
	}



	//SNS
	int count=0;
	do{
		reachedSingularity=false;
		count++;
		//ROS_INFO("%d",count);
		if (count>2*n_dof){
#ifndef _ONLY_WARNING_ON_ERROR
			ROS_ERROR("Infinite loop on SNS for task (%d)",priority);
/*
			ROS_INFO("p:%d  scale:%f  mc:%d  sing:%d",priority,scalingFactor,mostCriticalJoint,(int)reachedSingularity);
			//##############################
			string s;
			stringstream buffer;
			streambuf * old = std::cout.rdbuf(buffer.rdbuf());
			cout << "W\n" << W[priority] << endl << "P(k-1)\n" <<(*higherPriorityNull)<< endl<<"Psat\n" << projectorSaturated <<endl << "dotQ\n" << dotQopt[priority] <<endl << "JPinverse\n" << JPinverse <<endl;
			s = buffer.str();
			ROS_INFO("\n p %d \n %s",priority,s.c_str());
			//#############################
*/
			exit(1);
#else
			ROS_WARN("Infinite loop on SNS for task (%d)",priority);
/*
			ROS_INFO("p:%d  scale:%f  mc:%d  sing:%d",priority,scalingFactor,mostCriticalJoint,(int)reachedSingularity);
			//##############################
			string s;
			stringstream buffer;
			streambuf * old = std::cout.rdbuf(buffer.rdbuf());
			cout << "W\n" << W[priority]  <<endl;
			s = buffer.str();
			ROS_INFO("\n bs %f \n %s",scalingFactor,s.c_str());
			//#############################
*/
			// the task is not executed
			if (bestScale >= 0.0){
				W[priority]=bestW;
				dotQopt[priority]=(*higherPriorityJointVelocity) + bestInvJP *( bestScale*(*task)  - (*jacobian)*(*higherPriorityJointVelocity)) +bestTildeP*bestDotQn;
			}else{
				W[priority]=I;
				dotQopt[priority]=(*higherPriorityJointVelocity);
				//exit(1);
			}

			(*jointVelocity)=dotQopt[priority];
			return bestScale;
#endif
		}
		limit_excedeed=false;
		// If W=I everything is done --> go to the saturation phase

		if (!isW_identity && !invJPcomputed){
			if (priority==0){ //for the primary task higherPriorityNull==I
				barP=W[0];
				projectorSaturated=(I-W[0]);
			}else{
				temp=(I-W[priority]);//*(*higherPriorityNull);
				reachedSingularity|=!pinv_forBarP(&temp,higherPriorityNull,&projectorSaturated);
				//temp=(I-W[priority])*(*higherPriorityNull);
				//reachedSingularity|=pinv_damped(&temp,&projectorSaturated);
				barP=(I-projectorSaturated)*(*higherPriorityNull);
			}
			temp=(*jacobian)*barP;
			reachedSingularity|=!pinv(&temp, &JPinverse);

			invJPcomputed=false; //it needs to be computed for the next step
		}
		if (!isW_identity){

			tildeP=(I-JPinverse*(*jacobian)) * projectorSaturated;

			//compute the joint velocity
			dotQopt[priority]= (*higherPriorityJointVelocity) + JPinverse *( (*task)  - (*jacobian)*(*higherPriorityJointVelocity)) +tildeP*dotQn;

			//compute the scaling factor
			a=(JPinverse * (*task)).array();
			b=dotQopt[priority].array() - a;
			getTaskScalingFactor(&a, &b,&W[priority], &scalingFactor, &mostCriticalJoint);


			if ((scalingFactor >= 1.0)||(scalingFactor < 0)){

				//check optimality
				if(!isOptimal(priority, &dotQopt[priority], &tildeP,&W[priority],&dotQn)){
					//modified W and dotQn
					limit_excedeed=true;
					//ROS_INFO("non OPT");
					continue;
				}
			}
			if (scalingFactor >= 1.0){ //solution found

				//here limit_excedeed=false
				continue;
			}

			//if no solution found
			limit_excedeed=true;

			// is scaled an optimum
			if (scalingFactor >= 0){
				dotQs=(*higherPriorityJointVelocity) + JPinverse *( scalingFactor*(*task)  - (*jacobian)*(*higherPriorityJointVelocity)) +tildeP*dotQn;
				if(!isOptimal(priority, &dotQs, &tildeP,&W[priority],&dotQn)){
					//ROS_INFO("non OPT s");
					//modified W and dotQn
					limit_excedeed=true;
					continue;
				}
			}

			if (scalingFactor>bestScale){

				//save best solution so far
				bestScale=scalingFactor;
				//bestTildeDotQ=tildeDotQ;
				bestInvJP=JPinverse;
				bestW=W[priority];
				bestTildeP=tildeP;
				//bestPS=projectorSaturated;
				bestDotQn=dotQn;
			}

		}else limit_excedeed=true; //if it was W=I, to be here the limit is exceeded

		nSat[priority]++;
		W[priority](mostCriticalJoint,mostCriticalJoint)=0.0;
		isW_identity=false;
		if (dotQopt[priority](mostCriticalJoint)>dotQmax(mostCriticalJoint)){
			dotQn(mostCriticalJoint)=dotQmax(mostCriticalJoint)-(*higherPriorityJointVelocity)(mostCriticalJoint);
		}else{
			dotQn(mostCriticalJoint)=dotQmin(mostCriticalJoint)-(*higherPriorityJointVelocity)(mostCriticalJoint);
		}


		//compute JPinverse
		if (priority==0){ //for the primary task higherPriorityNull==I
			barP=W[0];
			projectorSaturated=(I-W[0]);
		}else{
			temp=(I-W[priority]);//*(*higherPriorityNull);
			reachedSingularity|=!pinv_forBarP(&temp,higherPriorityNull,&projectorSaturated);
			//temp=(I-W[priority])*(*higherPriorityNull);
			//reachedSingularity|=pinv_damped(&temp,&projectorSaturated);
			barP=(I-projectorSaturated)*(*higherPriorityNull);
		}
		temp=(*jacobian)*barP;
		reachedSingularity|=!pinv(&temp, &JPinverse);

		invJPcomputed=true;
		//if reachedSingularity then take the best solution

		if ((reachedSingularity)||(scalingFactor < 1e-12)){
			//ROS_INFO("sing best scale  %f",bestScale);

			if (bestScale >= 0.0){
				W[priority]=bestW;
				dotQopt[priority]=(*higherPriorityJointVelocity) + bestInvJP *( bestScale*(*task)  - (*jacobian)*(*higherPriorityJointVelocity)) +bestTildeP*bestDotQn;
			}else{
				//W[priority]=I;
				dotQopt[priority]=(*higherPriorityJointVelocity);
/*
				//second chance
				W[priority]=I;
				dotQn=VectorD::Zero(n_dof);
				invJPcomputed=false;
				limit_excedeed=true;
				continue;
*/
			}

			(*jointVelocity)=dotQopt[priority];
			return bestScale;
			//limit_excedeed=false;
			//continue;
		}


	}while(limit_excedeed);
	(*jointVelocity)=dotQopt[priority];
	return scalingFactor;
}

Scalar IKL::FSNSsingle(int priority,VectorD *higherPriorityJointVelocity,MatrixD *higherPriorityNull, MatrixD *jacobian, VectorD *task, VectorD *jointVelocity,MatrixD *nullSpaceProjector,double start_scale){
//THERE IS A PROBLEM if we use 3 tasks... to be checked



	//INITIALIZATION
	MatrixD JPinverse;					//(J_k P_{k-1})^#
	Array<Scalar,Dynamic,1> a,b;					   // used to compute the task scaling factor
	bool limit_excedeed;
	Scalar scalingFactor=1.0;
	int mostCriticalJoint;
	bool singularTask=false;

	Scalar best_Scale=-1.0;
	VectorD best_dq1;
	VectorD best_dq2;
	VectorD best_dqw;
	int best_nSat;


	MatrixD tildeZ;
	VectorD dq1,dq2,dqw;

	MatrixD bin,zin;
	Scalar dqw_in;
	Scalar error;


//initialization
	nSat[priority]=0;
	S[priority]=VectorXi::Zero(n_dof);



//compute the base solution

	singularTask=!pinv_QR_Z(jacobian, higherPriorityNull, &JPinverse, nullSpaceProjector);
	dq1=JPinverse * (*task);
	dq2=-JPinverse * (*jacobian)*(*higherPriorityJointVelocity);
	dqw=VectorD::Zero(n_dof);
	tildeZ=*nullSpaceProjector;

	dotQ=  (*higherPriorityJointVelocity) + dq1 +dq2;
	a=dq1.array();
	b=dotQ.array() - a;
	getTaskScalingFactor_fast(&a, &b,&S[priority], &scalingFactor, &mostCriticalJoint);

	if (scalingFactor >= 1.0){
		// then is clearly the optimum since all joints velocity are computed with the pseudoinverse
		(*jointVelocity)=dotQ;
		nSat[priority]=0;
		dotQopt[priority]=dotQ;
		return scalingFactor;
	}

	if (singularTask){
		// the task is singular so return a scaled damped solution (no SNS possible)
		//ROS_ERROR("sing");
		if (scalingFactor >= 0.0){
			nSat[priority]=0;
			(*jointVelocity)= (*higherPriorityJointVelocity) +  scalingFactor*dq1  + dq2;
			dotQopt[priority]=(*jointVelocity);
		}else{
			// the task is not executed
			nSat[priority]=0;
			*jointVelocity= *higherPriorityJointVelocity;
			dotQopt[priority]=(*jointVelocity);
			*nullSpaceProjector=*higherPriorityNull;
		}
		return scalingFactor;
	}

	if (scalingFactor>best_Scale){
//		ROS_INFO("best scale I %f",scalingFactor);
		//save best solution so far
		best_Scale=scalingFactor;
		best_dq1=dq1;
		best_dq2=dq2;
		best_dqw=dqw;
		best_nSat=0;

	}

//_______________________________________________________END of base computation


	//SNS
	int count=0;
	do{
		count++;
		//ROS_INFO("%d",count);
		if (count>2*n_dof){
#ifndef _ONLY_WARNING_ON_ERROR
			ROS_ERROR("Infinite loop on SNS for task (%d)",priority);

			exit(1);
#else
			ROS_WARN("Infinite loop on SNS for task (%d)",priority);

			// the task is not executed
			nSat[priority]=0;
			*jointVelocity= *higherPriorityJointVelocity;
			dotQopt[priority]=(*jointVelocity);
			*nullSpaceProjector=*higherPriorityNull;
			limit_excedeed=false;
			//continue;
			return -1.0;
#endif
		}
		limit_excedeed=true;

		//saturate the most critical joint
		zin=tildeZ.row(mostCriticalJoint);



		//solution non correct --> singularity
//		error=((*jacobian)*dotQ - (*task)).norm();
/*
		//##############################
		string s;
		stringstream buffer;
		streambuf * old = std::cout.rdbuf(buffer.rdbuf());
		cout << "err= " << error << " norm(zin)= " << zin.norm()<< "\n bin \n" << tildeZ*(zin/zin.squaredNorm()) << "\n Z \n" << tildeZ <<endl;
		s = buffer.str();
		ROS_INFO("\n bs %f idxW %d \n %s",scalingFactor,mostCriticalJoint,s.c_str());
		//#############################

*/
//		if ((error>1e-10) || (error!=error)){
		if ((zin.norm()<1e-8)||(scalingFactor<1e-12)){
			if (best_Scale>=0){
				//take the best solution
				*jointVelocity=(*higherPriorityJointVelocity) + best_Scale*best_dq1 + best_dq2 + best_dqw;
				dotQopt[priority]=(*jointVelocity);
				//nSat[priority]=best_nSat;
				*nullSpaceProjector=tildeZ; //if start net task from previous saturations
				return best_Scale;
			}else{
				//no solution
				//nSat[priority]=0;
				*jointVelocity= *higherPriorityJointVelocity;
				dotQopt[priority]=(*jointVelocity);
				*nullSpaceProjector=*higherPriorityNull;
				limit_excedeed=false;
				//continue;
				return -1.0;
			}
		}

		//if we do not use norm(zin) then this part can go first
		bin=tildeZ * (zin.transpose()/zin.squaredNorm());

		dq1-= bin*dq1(mostCriticalJoint);
		dq2-= bin*dq2(mostCriticalJoint);

		if (dotQ(mostCriticalJoint)>=0.0){
			dqw_in=dotQmax(mostCriticalJoint)-(*higherPriorityJointVelocity)(mostCriticalJoint);
		}else{
			dqw_in=dotQmin(mostCriticalJoint)-(*higherPriorityJointVelocity)(mostCriticalJoint);
		}
		dqw+= bin*(dqw_in - dqw(mostCriticalJoint));

		dotQ=(*higherPriorityJointVelocity)+dq1+dq2+dqw;

		nSat[priority]++;
		S[priority](mostCriticalJoint)=nSat[priority];
		tildeZ-=bin*zin;

		a=dq1.array();
		b=dotQ.array() - a;
		getTaskScalingFactor_fast(&a, &b,&S[priority], &scalingFactor, &mostCriticalJoint);
		if (scalingFactor >= 1.0){
			// task accomplished
			*jointVelocity= dotQ;
			dotQopt[priority]=(*jointVelocity);
			*nullSpaceProjector=tildeZ; //if start net task from previous saturations
			return scalingFactor;
		}else {
			if ((scalingFactor>best_Scale)){
				//save best solution so far
				best_Scale=scalingFactor;
				best_dq1=dq1;
				best_dq2=dq2;
				best_dqw=dqw;
				best_nSat=nSat[priority];
			}
		}


	}while(limit_excedeed); //actually in this implementation if we use while(1) it would be the same

	(*jointVelocity)=dotQ;
	return 1.0;
}

//#define LOG_ACTIVE

Scalar IKL::FOSNSsingle(int priority,VectorD *higherPriorityJointVelocity,MatrixD *higherPriorityNull, MatrixD *jacobian, VectorD *task, VectorD *jointVelocity,MatrixD *nullSpaceProjector,double start_scale){

	//INITIALIZATION
	MatrixD JPinverse;					//(J_k P_{k-1})^#
	Array<Scalar,Dynamic,1> a,b;					   // used to compute the task scaling factor
	bool limit_excedeed;
	Scalar scalingFactor=1.0;
	int mostCriticalJoint;
	bool singularTask=false;

	Scalar base_Scale;
	Scalar best_Scale=-1.0;
	VectorD best_dq1;
	VectorD best_dq2;
	VectorD best_dqw;
	int best_nSat;


	MatrixD tildeZ;
	VectorD dq1,dq2;
	VectorD dqw;
	VectorD dq1_base,dq2_base;
	VectorD dqn=VectorD::Zero(n_dof);
	VectorD dotQs;

	MatrixD zin;
	VectorD bin,bout;
	Scalar dqw_in;
	Scalar error;

	Scalar mu_in;
	Scalar mu_in1;
	Scalar mu_inp2w;
	Scalar mu_out;
	Scalar mu_out1;
	Scalar mu_outp2w;
	Scalar min_mu;
	int id_min_mu=n_dof+1; //just to be not possible
	VectorD scaledMU;

	bool computedScalingFactor;

#ifdef LOG_ACTIVE
	int n_in=0,n_out=0;
	stringstream log;
#endif

//compute the base solution

	singularTask=!pinv_QR_Z(jacobian, higherPriorityNull, &JPinverse, nullSpaceProjector);
	dq1=JPinverse * (*task);
	dq2=-JPinverse * (*jacobian)*(*higherPriorityJointVelocity);
	dqw=VectorD::Zero(n_dof);
	tildeZ=*nullSpaceProjector;
	dq1_base=dq1;
	dq2_base=dq2;
	dotQ=  (*higherPriorityJointVelocity) + dq1 +dq2;
	a=dq1.array();
	b=dotQ.array() - a;
	getTaskScalingFactor_basic(&a, &b, &scalingFactor, &mostCriticalJoint);

#ifdef LOG_ACTIVE
	log<<"task "<<priority<<endl;
	log<<"base Z norm "<<higherPriorityNull->norm()<<endl;
	log<<"base J*Z norm"<<((*jacobian)*(*higherPriorityNull)).norm()<<endl;
	log<<"scale factor at 0 "<<scalingFactor<<endl;
	log<<"base S "<<S[priority].transpose()<<endl;
#endif

	if (scalingFactor >= 1.0){
		// then is clearly the optimum since all joints velocity are computed with the pseudoinverse
		(*jointVelocity)=dotQ;
		//dotQopt[priority]=dotQ;
		nSat[priority]=0;
		satList[priority].clear();
		S[priority]=VectorXi::Zero(n_dof);
#ifdef LOG_ACTIVE
		log<<"task accomplished without saturations"<<endl;
		log<<"scale "<<scalingFactor<<endl;
		//log<<"last dotQ "<<higherPriorityJointVelocity->transpose()<<endl;
		//log<<"dotQ "<<dotQ.transpose()<<endl;
		//log<<"dotQmin "<<dotQmin.transpose()<<endl;
		//log<<"dotQmax "<<dotQmax.transpose()<<endl;
		if (singularTask) log<<"THE TASK IS SINGULAR"<<endl;
		ROS_WARN("\n%s\n\n",log.str().c_str());
#endif
		return scalingFactor;
	}else{
		base_Scale=scalingFactor;
		if ((scalingFactor>best_Scale)){
			//save best solution so far
			best_Scale=scalingFactor;
			best_dq1=dq1;
			best_dq2=dq2;
			best_dqw=dqw;
			best_nSat=nSat[priority];
		}
	}

	if ((singularTask)||(base_Scale<0)){
		// the task is singular so return a scaled damped solution (no SNS possible)


		if (scalingFactor > 0.0){
			(*jointVelocity)= (*higherPriorityJointVelocity) +  scalingFactor*dq1  + dq2;
			//dotQopt[priority]=(*jointVelocity);
		}else{
			// the task is not executed
			*jointVelocity= *higherPriorityJointVelocity;
			//dotQopt[priority]=(*jointVelocity);
			*nullSpaceProjector=*higherPriorityNull;
		}
#ifdef LOG_ACTIVE
		if (singularTask) log<<"the task is singular"<<endl;
		if (base_Scale<0) log<<"base scale < 0"<<endl;
		log<<"scale "<<scalingFactor<<endl;
	//	log<<"dotQ "<<dotQ.transpose();
		//log<<"dotQmin "<<dotQmin.transpose();
	//	log<<"dotQmax "<<dotQmax.transpose();
		ROS_WARN("\n%s\n\n",log.str().c_str());
#endif
		nSat[priority]=0;
		satList[priority].clear();
		S[priority]=VectorXi::Zero(n_dof);
		return scalingFactor;
	}


//_______________________________________________________END of base computation

#define WARM_START
#ifdef	WARM_START
//############### WARM START
//B=MatrixD::Zero(n_dof,n_dof);
	if (nSat[priority]){
//		log<<"started with "<<nSat[priority]<< " saturated\n";
		MatrixD Zws=MatrixD::Zero(nSat[priority],tildeZ.cols());
		MatrixD invZws=MatrixD::Zero(tildeZ.cols(),nSat[priority]);
		MatrixD forMu=MatrixD::Zero(nSat[priority],nSat[priority]);
		MatrixD Bws=MatrixD::Zero(n_dof,nSat[priority]);
		VectorD dq1_ws=VectorD::Zero(nSat[priority]);
		VectorD dq2_ws=VectorD::Zero(nSat[priority]);
		VectorD dqw_ws=VectorD::Zero(nSat[priority]);
		int idws=0;
		bool invertibelZws;

		for (it=satList[priority].begin(); it!=satList[priority].end(); ){
			int id=*it;

			if (tildeZ.row(id).norm()<1e-10){
				//it means that the joint has been already saturated by the previous task
				if (it==satList[priority].begin()){
					it=satList[priority].erase_after(satList[priority].before_begin());
				}else{
					it=satList[priority].erase_after(prev_it);
				}
				S[priority](id)=0;
				nSat[priority]--;
			}else {

#ifdef LOG_ACTIVE
				n_in++;
				log<<id<< " ";
#endif

				Zws.row(idws)=tildeZ.row(id);
				dq1_ws(idws)=dq1_base(id);
				dq2_ws(idws)=dq2_base(id);
				if (S[priority](id)>0.0){
					dqw_ws(idws)=dotQmax(id)-(*higherPriorityJointVelocity)(id);
				}else{
					dqw_ws(idws)=dotQmin(id)-(*higherPriorityJointVelocity)(id);
				}
				dqn(id)=dqw_ws(idws);

				idws++;
				prev_it=it;
				it++;
			}

		}


		if (nSat[priority]){
			Zws.conservativeResize(nSat[priority],tildeZ.cols());
			invertibelZws=pinv_QR(&Zws,&invZws);



			if (!invertibelZws){
			//	ROS_WARN("Zws is not invertible... what should I do?");
#ifdef LOG_ACTIVE_NO
					//log<<endl<<"Z \n"<<tildeZ;
					//log<<endl<<"Zws \n"<<Zws;
					log<<"\nS\n"<<S[priority].transpose()<<endl<<"norm Zws ";
					for (int i=0;i<nSat[priority];i++){
					//for (it=satList[priority].begin(); it!=satList[priority].end(); ++it){
					//	int id=*it;
					//	log<<" "<<tildeZ.row(id).norm();
						log<<" "<<Zws.row(i).norm();
					}
					log<<"\n Zws rows "<<Zws.rows();
					ROS_WARN("%s",log.str().c_str());
					exit(1);
#endif
				satList[priority].clear();
				nSat[priority]=0;
				S[priority]=VectorXi::Zero(n_dof);

			}else{

				Bws=tildeZ*invZws;
				forMu=invZws.transpose()*invZws;
	//			forMu=Bws.transpose()*Bws;
				idws=0;
				for (it=satList[priority].begin(); it!=satList[priority].end(); ++it){
					int id=*it;
					B.col(id)=Bws.col(idws);

					lagrangeMu1(id)=-forMu.row(idws)*dq1_ws;
					lagrangeMup2w(id)=forMu.row(idws)*(dqw_ws-dq2_ws);
					lagrangeMu(id)=lagrangeMu1(id)+lagrangeMup2w(id);

					idws++;
				}
				dq1=dq1_base-Bws*dq1_ws;
				dq2=dq2_base-Bws*dq2_ws;
				dqw=Bws*dqw_ws;
				dotQ=(*higherPriorityJointVelocity)+dq1+dq2+dqw;


				tildeZ=*nullSpaceProjector - Bws*Zws;
	#ifdef LOG_ACTIVE
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
				a=dq1.array();
				b=dotQ.array() - a;
				getTaskScalingFactor_fast(&a, &b,&S[priority], &scalingFactor, &mostCriticalJoint);
				scaledMU=scalingFactor*lagrangeMu1 + lagrangeMup2w;
				//find the minimum negative mu
				for (it=satList[priority].begin(); it!=satList[priority].end(); ++it){
					int id=*it;
					if (dqn(id)>=0) scaledMU(id)=-scaledMU(id);
				}
				log<<"/nstart scale "<< scalingFactor<<endl;
				//log<<"/nstart scaled Mu "<< scaledMU.transpose()<<endl;
				log<<"/nstart scaled dq "<< ((*higherPriorityJointVelocity)+scalingFactor*dq1+dq2+dqw).transpose()<<endl;
				log<<"/nstart S "<<S[priority].transpose();
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
	#endif
				//find the minimum negative mu
				min_mu=0;
				id_min_mu=n_dof+1;
				for (it=satList[priority].begin(); it!=satList[priority].end(); ++it){
					int id=*it;
					if (dqn(id)>=0) lagrangeMu(id)=-lagrangeMu(id);
					if ((lagrangeMu(id)<min_mu)&&(abs(dqn(id))> 1e-12)){
						min_mu=lagrangeMu(id);
						id_min_mu=id;
					}
				}
	#ifdef LOG_ACTIVE
				log<<"/nstart Mu "<< lagrangeMu.transpose()<<endl;
	#endif
			}
		}
	}
	computedScalingFactor=false;
#else
//######################### else
	satList[priority].clear();
    nSat[priority]=0;
	S[priority]=VectorXi::Zero(n_dof);

//	lagrangeMu1=VectorD::Zero(n_dof);
//	lagrangeMup2w=VectorD::Zero(n_dof);
//#########################
#endif

	//SNS
	int count=0;
	do{
		count++;
		//ROS_INFO("%d",count);
		if (count>2*n_dof){
#ifndef _ONLY_WARNING_ON_ERROR
			ROS_ERROR("Infinite loop on SNS for task (%d)",priority);

			exit(1);
#else
			ROS_WARN("Infinite loop on SNS for task (%d): nSat=%d ",priority,nSat[priority]);
/*			//##############################
			string s;
			stringstream buffer;
			streambuf * old = std::cout.rdbuf(buffer.rdbuf());
			cout << "nSat\n" << nSat[priority]<<endl;
			s = buffer.str();
			ROS_INFO("\n p %d \n %s",priority,s.c_str());
			//#############################
*/
			// the task is not executed
			//nSat[priority]=0;
			satList[priority].clear();
			nSat[priority]=0;
			S[priority]=VectorXi::Zero(n_dof);
			*jointVelocity= *higherPriorityJointVelocity;
			//dotQopt[priority]=(*jointVelocity);
			*nullSpaceProjector=*higherPriorityNull;
			limit_excedeed=false;
#ifdef LOG_ACTIVE
			//log<<endl<<"last dq "<<dotQ.transpose()<< "\ndq1"<<dq1.transpose()<< "\ndq2"<<dq2.transpose()  << "\ndqn"<<dqn.transpose()<< "\ndqw"<<dqw.transpose()<< "\n\nS"<<S[priority].transpose()<<endl;
			log<<"last scaling factor "<< scalingFactor;
			ROS_WARN("%s",log.str().c_str());
			exit(1);
#endif
			//continue;
			return -1.0;
#endif
		}
		limit_excedeed=true;

		if (!computedScalingFactor){
			a=dq1.array();
			b=dotQ.array() - a;
			getTaskScalingFactor_fast(&a, &b,&S[priority], &scalingFactor, &mostCriticalJoint);
		}
		computedScalingFactor=false;

#ifdef LOG_ACTIVE
		//if (n_in>2*n_dof-3){
//		log<<endl<<"last dq "<<dotQ.transpose()<< "\ndq1"<<dq1.transpose()<< "\ndq2"<<dq2.transpose()  << "\ndqn"<<dqn.transpose()<< "\ndqw"<<dqw.transpose()<< "\n\nS"<<S[priority].transpose()<<endl;
			log<<" last scaling factor "<< scalingFactor;
			log<<"\n n_in "<< n_in<<" n_out"<<n_out<<"  -> S "<<S[priority].transpose();
//			log<<"\n nSat"<< nSat[priority];

			log<<"\n best scale factor "<< best_Scale<<"\n\n";

			scaledMU=scalingFactor*lagrangeMu1 + lagrangeMup2w;
			//find the minimum negative mu
			for (it=satList[priority].begin(); it!=satList[priority].end(); ++it){
				int id=*it;
				if (dqn(id)>=0) scaledMU(id)=-scaledMU(id);
			}
			//log<<"/nstop scaled Mu "<< scaledMU.transpose()<<endl;
			//log<<"stop scaled dq "<< (*higherPriorityJointVelocity)+scalingFactor*dq1+dq2+dqw<<endl;
			log<<"error"<<((*jacobian)*dotQ - (*task)).norm();

			ROS_WARN("\n%s\n\n",log.str().c_str());
			log.str("");
	//		exit(1);
	//	}
#endif


		if ((scalingFactor >= 1.0)||(scalingFactor <0.0)){
			//check the optimality of the solution
			if (id_min_mu<n_dof){
#ifdef LOG_ACTIVE
				n_out++;
				log<<" O"<<id_min_mu;
#endif
				//remove the id_min_mu joint from saturation and update B
				bout=B.col(id_min_mu);
				mu_out1=lagrangeMu1(id_min_mu);
				mu_outp2w=lagrangeMup2w(id_min_mu);
				for (it=satList[priority].begin(); it!=satList[priority].end(); ){
					int id=*it;
					if (id==id_min_mu){
						//remove id-min_mu from list of saturated joints
						if (it==satList[priority].begin()){
							it=satList[priority].erase_after(satList[priority].before_begin());
						}else{
							it=satList[priority].erase_after(prev_it);
						}
					}else {
						//update B
						Scalar baux=(Scalar) bout.dot(B.col(id))/bout.squaredNorm();
						B.col(id)-=bout*baux;
						prev_it=it;
						it++;
					}
				}
				//update the solution
				Scalar dq1X=dq1_base(id_min_mu);
				Scalar dq2X=dq2_base(id_min_mu);
				Scalar dqwX=dqn(id_min_mu);
				MatrixD ZX=nullSpaceProjector->row(id_min_mu);
				for (it=satList[priority].begin(); it!=satList[priority].end(); ++it){
					int id=*it;
					lagrangeMu1(id)+=bout(id)*mu_out1;
					lagrangeMup2w(id)+=bout(id)*mu_outp2w;
					lagrangeMu(id) =lagrangeMu1(id) + lagrangeMup2w(id);
					dq1X-=B(id_min_mu,id)*dq1_base(id);
					dq2X-=B(id_min_mu,id)*dq2_base(id);
					dqwX-=B(id_min_mu,id)*dqn(id);
					ZX-=B(id_min_mu,id)*nullSpaceProjector->row(id);

				}
				dq1+=bout*dq1X;
				dq2+=bout*dq2X;
				dqw-=bout*dqwX;
				tildeZ+=bout*ZX;
				nSat[priority]--;
				S[priority](id_min_mu)=0;
				dotQ=(*higherPriorityJointVelocity)+dq1+dq2+dqw;

				//find the minimum negative mu
				min_mu=0;
				id_min_mu=n_dof+1;
				for (it=satList[priority].begin(); it!=satList[priority].end(); ++it){
					int id=*it;
					if (dqn(id)>=0) lagrangeMu(id)=-lagrangeMu(id);
					if ((lagrangeMu(id)<min_mu)&&(abs(dqn(id))> 1e-12)){
						min_mu=lagrangeMu(id);
						id_min_mu=id;
					}
				}
				continue;
			}
		}
		if (scalingFactor >= 1.0){
			// task accomplished
			*jointVelocity= dotQ;
			//dotQopt[priority]=(*jointVelocity);
			*nullSpaceProjector=tildeZ; //if start net task from previous saturations
//    		ROS_INFO("n_in %d n_out %d",n_in,n_out);

			return scalingFactor;

		}else if ((best_Scale > 0.0)||(scalingFactor > 0.0)){
			if ((scalingFactor>best_Scale)){
				//save best solution so far
				best_Scale=scalingFactor;
				best_dq1=dq1;
				best_dq2=dq2;
				best_dqw=dqw;
				best_nSat=nSat[priority];
			}
			//check if the scaled solution is the optimum
			scaledMU=scalingFactor*lagrangeMu1 + lagrangeMup2w;
			//find the minimum negative mu
			min_mu=0;
			id_min_mu=n_dof+1;
			for (it=satList[priority].begin(); it!=satList[priority].end(); ++it){
				int id=*it;
				if (dqn(id)>=0) scaledMU(id)=-scaledMU(id);
				//if ((scaledMU(id)<0) && (scaledMU(id)>min_mu)){
				if ((scaledMU(id)<min_mu)&&(abs(dqn(id))> 1e-12)){
					min_mu=scaledMU(id);
					id_min_mu=id;
				}
			}

			if (id_min_mu<n_dof){
#ifdef LOG_ACTIVE
				n_out++;
				log<<" Os"<<id_min_mu;
#endif
				//remove the id_min_mu joint from saturation and update B
				bout=B.col(id_min_mu);
				mu_out1=lagrangeMu1(id_min_mu);
				mu_outp2w=lagrangeMup2w(id_min_mu);
				for (it=satList[priority].begin(); it!=satList[priority].end(); ){
					int id=*it;
					if (id==id_min_mu){
						//remove id-min_mu from list of saturated joints
						if (it==satList[priority].begin()){
							it=satList[priority].erase_after(satList[priority].before_begin());
						}else{
							it=satList[priority].erase_after(prev_it);
						}
					}else {
						//update B
						Scalar baux=(Scalar) bout.dot(B.col(id))/bout.squaredNorm();
						B.col(id)-=bout*baux;
						prev_it=it;
						it++;
					}
				}
				//update the solution
				Scalar dq1X=dq1_base(id_min_mu);
				Scalar dq2X=dq2_base(id_min_mu);
				Scalar dqwX=dqn(id_min_mu);
				MatrixD ZX=nullSpaceProjector->row(id_min_mu);
				for (it=satList[priority].begin(); it!=satList[priority].end(); ++it){
					int id=*it;
					lagrangeMu1(id)+=bout(id)*mu_out1;
					lagrangeMup2w(id)+=bout(id)*mu_outp2w;
					lagrangeMu(id) =lagrangeMu1(id) + lagrangeMup2w(id);
					dq1X-=B(id_min_mu,id)*dq1_base(id);
					dq2X-=B(id_min_mu,id)*dq2_base(id);
					dqwX-=B(id_min_mu,id)*dqn(id);
					ZX-=B(id_min_mu,id)*nullSpaceProjector->row(id);

				}
				dq1+=bout*dq1X;
				dq2+=bout*dq2X;
				dqw-=bout*dqwX;
				tildeZ+=bout*ZX;
				nSat[priority]--;
				S[priority](id_min_mu)=0;
				dotQ=(*higherPriorityJointVelocity)+dq1+dq2+dqw;

				//find the minimum negative mu
				min_mu=0;
				id_min_mu=n_dof+1;
				for (it=satList[priority].begin(); it!=satList[priority].end(); ++it){
					int id=*it;
					if (dqn(id)>=0) lagrangeMu(id)=-lagrangeMu(id);
					if ((lagrangeMu(id)<min_mu)&&(abs(dqn(id))> 1e-12)){
						min_mu=lagrangeMu(id);
						id_min_mu=id;
					}
				}
				continue;
			}

		}


		int idxW=mostCriticalJoint;

		//saturate the most critical joint
		zin=tildeZ.row(idxW);
		//if we do not use norm(zin) then this part can go first
		bin=tildeZ * (zin.transpose()/zin.squaredNorm());

		dq1-= bin*dq1(idxW);
		dq2-= bin*dq2(idxW);

		if (dotQ(idxW)>=0.0){
			dqw_in=dotQmax(idxW)-(*higherPriorityJointVelocity)(idxW);
			S[priority](idxW)=+1;
		}else{
			dqw_in=dotQmin(idxW)-(*higherPriorityJointVelocity)(idxW);
			S[priority](idxW)=-1;
		}
		dqw+= bin*(dqw_in - dqw(idxW));
		dqn(idxW)=dqw_in;

		dotQ=(*higherPriorityJointVelocity)+dq1+dq2+dqw;


		a=dq1.array();
		b=dotQ.array() - a;
		getTaskScalingFactor_fast(&a, &b,&S[priority], &scalingFactor, &mostCriticalJoint);
		computedScalingFactor=true;

#ifdef LOG_ACTIVE
		log<<" I "<<idxW<< "norm zin "<<zin.norm()<<" at "<<dqw_in<<" with scale "<<scalingFactor<<endl;
		if (scalingFactor<0) log<<"IT WAS < 0 "<<endl;
		if (scalingFactor<1e-12) log<<"IT WAS < eps "<<endl;
		n_in++;
	//	ROS_WARN("\n%s\n\n",log.str().c_str());
	//	log.str("");

#endif
		if ((zin.norm()<1e-8) ||(scalingFactor<1e-12)){
			S[priority](idxW)=0; //put it back to the right value
			if (best_Scale>=0){
				//take the best solution
				*jointVelocity=(*higherPriorityJointVelocity) + best_Scale*best_dq1 + best_dq2 + best_dqw;
				//dotQopt[priority]=(*jointVelocity);
				//nSat[priority]=best_nSat;
				*nullSpaceProjector=tildeZ; //if start net task from previous saturations
//				ROS_INFO("n_in %d n_out %d",n_in,n_out);
				if (best_Scale==base_Scale){
					//no saturation was needed to obtain the best scale
					satList[priority].clear();
					nSat[priority]=0;
					S[priority]=VectorXi::Zero(n_dof);

				}
				//if (priority==1) *jointVelocity=(*higherPriorityJointVelocity);
				return best_Scale;
			}else{
				//no solution
				//nSat[priority]=0;
				*jointVelocity= *higherPriorityJointVelocity;
				//dotQopt[priority]=(*jointVelocity);
				*nullSpaceProjector=*higherPriorityNull;
				limit_excedeed=false;
				//continue;
#ifdef LOG_ACTIVE_UNFEASIBLE
				log<<endl<<"last dq "<<dotQ.transpose()<< "\ndq1"<<dq1.transpose()<< "\ndq2"<<dq2.transpose()  << "\ndqn"<<dqn.transpose()<< "\ndqw"<<dqw.transpose()<< "\n\nS"<<S[priority].transpose()<<endl;
				log<<"last scaling factor "<< scalingFactor;
				ROS_WARN("%s",log.str().c_str());
				exit(1);
#endif
				return -1.0;
			}
		}



		nSat[priority]++;
		tildeZ-=bin*zin;

		//update mu and B
		min_mu=0;
		id_min_mu=n_dof+1;
		mu_in1=bin.dot(dq1);
		mu_inp2w=bin.dot(dotQ-dq1);
		for (it=satList[priority].begin(); it!=satList[priority].end(); ++it){
			int id=*it;
			lagrangeMu1(id)-=B(idxW,id)*mu_in1;
			lagrangeMup2w(id)-=B(idxW,id)*mu_inp2w;
			lagrangeMu(id)=lagrangeMu1(id)+lagrangeMup2w(id);
			B.col(id)-=bin*B(idxW,id);
			if (dqn(id)>=0) lagrangeMu(id)=-lagrangeMu(id);
			if ((lagrangeMu(id)<min_mu)&&(abs(dqn(id))> 1e-12)){
				min_mu=lagrangeMu(id);
				id_min_mu=id;
			}
		}
		satList[priority].push_front(idxW);
		lagrangeMu1(idxW)=mu_in1;
		lagrangeMup2w(idxW)=mu_inp2w;
		//lagrangeMu(idxW)=lagrangeMu1(idxW)+lagrangeMup2w(idxW);

	}while(limit_excedeed); //actually in this implementation if we use while(1) it would be the same

	(*jointVelocity)=dotQ;
	//dotQopt[priority]=(*jointVelocity);
	return scalingFactor;
}



Scalar QPsingle(int priority,VectorD *higherPriorityJointVelocity,MatrixD *higherPriorityAugmentedJacobian, MatrixD *jacobian, VectorD *task, VectorD *jointVelocity,MatrixD *augmentedJacobian,double start_scale){

	//QuadProgPP::Matrix<double> G,CE,CI;
	//QuadProgPP::Vector<double> g0, ce0, ci0, x;

return 1.0;

}


Scalar IKL::getJointVelocity_SNS(VectorD *jointVelocity, StackOfTasks *sot, VectorD *jointConfiguration){


	  int n_task=sot->size();
	  int robotDOF=(*sot)[0].jacobian.cols();

	  //P_0=I
	  //dq_0=0
	  MatrixD P = MatrixD::Identity(robotDOF,robotDOF);
	  *jointVelocity = VectorD::Zero(robotDOF,1);
	  VectorD higherPriorityJointVelocity;
	  MatrixD higherPriorityNull;

	  shapeJointVelocityBound(jointConfiguration);

	  for(int i_task=0; i_task<n_task; i_task++){ //consider all tasks

		  higherPriorityJointVelocity=*jointVelocity;
		  higherPriorityNull=P;
		  scaleFactors[i_task]=SNSsingle(i_task,&higherPriorityJointVelocity,&higherPriorityNull, &((*sot)[i_task].jacobian), &((*sot)[i_task].desired),jointVelocity,&P);

	  }

/*
	  VectorD obtained= ((*sot)[0].jacobian)* (*jointVelocity);
	  VectorD desired= scaleFactors[0]*((*sot)[0].desired);
	  if (obtained(0)-desired(0) > 1e-8){
			  ROS_INFO("SNS error high\n scale factors: %f \n errors : %f ",scaleFactors[0],( ((*sot)[0].jacobian)* (*jointVelocity) - scaleFactors[0]*((*sot)[0].desired)  ).norm()  );
				string s;
				stringstream buffer;
				streambuf * old = std::cout.rdbuf(buffer.rdbuf());
				cout << "W\n" << W[0]  <<endl;
				s = buffer.str();
				ROS_INFO("\n %s",s.c_str());
		  }
*/
//ROS_WARN("%f - %f",scaleFactors[0],scaleFactors[1]);
	  return 1.0;

}


Scalar IKL::getJointVelocity_OSNS(VectorD *jointVelocity, StackOfTasks *sot, VectorD *jointConfiguration){


	  int n_task=sot->size();
	  int robotDOF=(*sot)[0].jacobian.cols();

	  //P_0=I
	  //dq_0=0
	  MatrixD P = MatrixD::Identity(robotDOF,robotDOF);
	  *jointVelocity = VectorD::Zero(robotDOF,1);
	  VectorD higherPriorityJointVelocity;
	  MatrixD higherPriorityNull;

	  shapeJointVelocityBound(jointConfiguration);

	  for(int i_task=0; i_task<n_task; i_task++){ //consider all tasks
		  higherPriorityJointVelocity=*jointVelocity;
		  higherPriorityNull=P;
		  scaleFactors[i_task]=OSNSsingle(i_task,&higherPriorityJointVelocity,&higherPriorityNull, &((*sot)[i_task].jacobian), &((*sot)[i_task].desired),jointVelocity,&P);
/*
		  if (scaleFactors[i_task]<0){
			  //second chance
			  W[i_task]=I;
			  P=higherPriorityNull;
			  scaleFactors[i_task]=OSNSsingle(i_task,&higherPriorityJointVelocity,&higherPriorityNull, &((*sot)[i_task].jacobian), &((*sot)[i_task].desired),jointVelocity,&P);

		  }
*/
		  if ( scaleFactors[i_task]>1)  scaleFactors[i_task]=1;
	  }

/*
	  if (( ((*sot)[0].jacobian)* (*jointVelocity) - scaleFactors[0]*((*sot)[0].desired)  ).norm() > 1e-12){
		  ROS_INFO("SNS error high\n scale factors: %f , %f\n errors : %f , %f",scaleFactors[0],scaleFactors[1],( ((*sot)[0].jacobian)* (*jointVelocity) - scaleFactors[0]*((*sot)[0].desired)  ).norm(),( ((*sot)[1].jacobian)* (*jointVelocity) - scaleFactors[1]*((*sot)[1].desired)  ).norm()  );
			string s;
			stringstream buffer;
			streambuf * old = std::cout.rdbuf(buffer.rdbuf());
			cout << "W\n" << W[0]  <<endl;
			s = buffer.str();
			ROS_INFO("\n %s",s.c_str());
	  }
*/
/*
	  if (( ((*sot)[0].jacobian)* (*jointVelocity) - scaleFactors[0]*((*sot)[0].desired)  ).norm() > 1e-8){
			  ROS_INFO("SNS error high\n scale factors: %f \n errors : %f ",scaleFactors[0],( ((*sot)[0].jacobian)* (*jointVelocity) - scaleFactors[0]*((*sot)[0].desired)  ).norm()  );
				string s;
				stringstream buffer;
				streambuf * old = std::cout.rdbuf(buffer.rdbuf());
				cout << "W\n" << W[0]  <<endl;
				s = buffer.str();
				ROS_INFO("\n %s",s.c_str());
		  }
*/

//ROS_WARN("%f - %f",scaleFactors[0],scaleFactors[1]);
	  return (double)1.0*nSat[0]+nSat[1]+nSat[2]+nSat[3]+nSat[4];

}


Scalar IKL::getJointVelocity_OSNSsm(VectorD *jointVelocity, StackOfTasks *sot, VectorD *jointConfiguration){

	  int n_task=sot->size();
	  int robotDOF=(*sot)[0].jacobian.cols();

	  //P_0=I
	  //dq_0=0
	  MatrixD P = MatrixD::Identity(robotDOF,robotDOF);
	  *jointVelocity = VectorD::Zero(robotDOF,1);
	  VectorD higherPriorityJointVelocity;
	  MatrixD higherPriorityNull;

	  MatrixD PS = MatrixD::Identity(robotDOF,robotDOF);

	  shapeJointVelocityBound(jointConfiguration);

	  for(int i_task=0; i_task<n_task; i_task++){ //consider all tasks
		  higherPriorityJointVelocity=*jointVelocity;
		  higherPriorityNull=P;


		  scaleFactors[i_task]=OSNSsingle(i_task,&higherPriorityJointVelocity,&higherPriorityNull, &((*sot)[i_task].jacobian), &((*sot)[i_task].desired),jointVelocity,&PS);
		  if (scaleFactors[i_task]<0){
			  //second chance
			  W[i_task]=I;
			  PS=higherPriorityNull;
			  scaleFactors[i_task]=OSNSsingle(i_task,&higherPriorityJointVelocity,&higherPriorityNull, &((*sot)[i_task].jacobian), &((*sot)[i_task].desired),jointVelocity,&PS);

		  }
		  if (scaleFactors[i_task] > 0.0){
			  if (scaleFactors[i_task]*scaleMargin < (1.0 )){
				  double taskScale=scaleFactors[i_task]*scaleMargin;
				  VectorD scaledTask=((*sot)[i_task].desired)*taskScale;
				  scaleFactors[i_task]=OSNSsingle(i_task,&higherPriorityJointVelocity,&higherPriorityNull, &((*sot)[i_task].jacobian), &scaledTask,jointVelocity,&P);
				  scaleFactors[i_task]=taskScale;

			  }else{
				  scaleFactors[i_task]=1.0;
				  P=PS;
		 	}
		  }
/*
		  if (scaleFactors[i_task] > scaleMargin){
			  if (scaleFactors[i_task] < (1.0 + scaleMargin)){
				  double taskScale=scaleFactors[i_task]-scaleMargin;
				  VectorD scaledTask=((*sot)[i_task].desired)*taskScale;
				  scaleFactors[i_task]=OSNSsingle(i_task,&higherPriorityJointVelocity,&higherPriorityNull, &((*sot)[i_task].jacobian), &scaledTask,jointVelocity,&P);
				  scaleFactors[i_task]=taskScale;

			  }else{
				  scaleFactors[i_task]=1.0;
		 	}
		  }
*/
	  }

/*
	  if (( ((*sot)[1].jacobian)* (*jointVelocity) - scaleFactors[1]*((*sot)[1].desired)  ).norm() > 1e-12){
		  ROS_INFO("SNS error high\n scale factors: %f , %f\n errors : %f , %f",scaleFactors[0],scaleFactors[1],( ((*sot)[0].jacobian)* (*jointVelocity) - scaleFactors[0]*((*sot)[0].desired)  ).norm(),( ((*sot)[1].jacobian)* (*jointVelocity) - scaleFactors[1]*((*sot)[1].desired)  ).norm()  );
			string s;
			stringstream buffer;
			streambuf * old = std::cout.rdbuf(buffer.rdbuf());
			cout << "W\n" << W[0]  <<endl;
			s = buffer.str();
			ROS_INFO("\n %s",s.c_str());
	  }
*/
	  /*
	  if (( ((*sot)[0].jacobian)* (*jointVelocity) - scaleFactors[0]*((*sot)[0].desired)  ).norm() > 1e-8){
			  ROS_INFO("SNS error high\n scale factors: %f \n errors : %f ",scaleFactors[0],( ((*sot)[0].jacobian)* (*jointVelocity) - scaleFactors[0]*((*sot)[0].desired)  ).norm()  );
				string s;
				stringstream buffer;
				streambuf * old = std::cout.rdbuf(buffer.rdbuf());
				cout << "W\n" << W[0]  <<endl;
				s = buffer.str();
				ROS_INFO("\n %s",s.c_str());
		  }
*/

//ROS_WARN("%f - %f",scaleFactors[0],scaleFactors[1]);
	  return 1.0;

}

Scalar IKL::getJointVelocity_FSNS(VectorD *jointVelocity, StackOfTasks *sot, VectorD *jointConfiguration){


	  int n_task=sot->size();
	  int robotDOF=(*sot)[0].jacobian.cols();

	  //P_0=I
	  //dq_0=0
	  MatrixD P = MatrixD::Identity(robotDOF,robotDOF);
	  *jointVelocity = VectorD::Zero(robotDOF,1);
	  VectorD higherPriorityJointVelocity;
	  MatrixD higherPriorityNull;

	  shapeJointVelocityBound(jointConfiguration);


	  for(int i_task=0; i_task<n_task; i_task++){ //consider all tasks
		  higherPriorityJointVelocity=*jointVelocity;
		  higherPriorityNull=P;
		  scaleFactors[i_task]=FSNSsingle(i_task,&higherPriorityJointVelocity,&higherPriorityNull, &((*sot)[i_task].jacobian), &((*sot)[i_task].desired),jointVelocity,&P);
		  if ( scaleFactors[i_task]>1)  scaleFactors[i_task]=1;
	  }


//	    ROS_INFO("s0 %f, s1 %f",scaleFactors[0],scaleFactors[1]);
/*

	  if (( ((*sot)[1].jacobian)* (*jointVelocity) - scaleFactors[1]*((*sot)[1].desired)  ).norm() > 1e-12){
		  ROS_INFO("SNS error high\n scale factors: %f , %f\n errors : %f , %f",scaleFactors[0],scaleFactors[1],( ((*sot)[0].jacobian)* (*jointVelocity) - scaleFactors[0]*((*sot)[0].desired)  ).norm(),( ((*sot)[1].jacobian)* (*jointVelocity) - scaleFactors[1]*((*sot)[1].desired)  ).norm()  );
	  }
*/
/*
	  if (( ((*sot)[0].jacobian)* (*jointVelocity) - scaleFactors[0]*((*sot)[0].desired)  ).norm() > 1e-10){
			  ROS_INFO("SNS error high\n scale factors: %f \n errors : %f  \n nSat %d",scaleFactors[0],( ((*sot)[0].jacobian)* (*jointVelocity) - scaleFactors[0]*((*sot)[0].desired)  ).norm(),nSat[0]  );
				string s;
				stringstream buffer;
				streambuf * old = std::cout.rdbuf(buffer.rdbuf());
				cout << "B\n" << B  <<endl;
				s = buffer.str();
				ROS_INFO("\n %s",s.c_str());
		  }
*/

//ROS_WARN("%f - %f",scaleFactors[0],scaleFactors[1]);
	  return (double)1.0*nSat[0]+nSat[1];

}

Scalar IKL::getJointVelocity_FOSNS(VectorD *jointVelocity, StackOfTasks *sot, VectorD *jointConfiguration){


	  int n_task=sot->size();
	  int robotDOF=(*sot)[0].jacobian.cols();

	  //P_0=I
	  //dq_0=0
	  MatrixD P = MatrixD::Identity(robotDOF,robotDOF);
	  MatrixD PS = MatrixD::Identity(robotDOF,robotDOF);
	  *jointVelocity = VectorD::Zero(robotDOF,1);
	  VectorD higherPriorityJointVelocity;
	  MatrixD higherPriorityNull;

	  shapeJointVelocityBound(jointConfiguration);

	  // this is not the best solution... the scale margin should be computed inside FOSNSsingle

	  for(int i_task=0; i_task<n_task; i_task++){ //consider all tasks
		  higherPriorityJointVelocity=*jointVelocity;
		  higherPriorityNull=P;
		  scaleFactors[i_task]=FOSNSsingle(i_task,&higherPriorityJointVelocity,&higherPriorityNull, &((*sot)[i_task].jacobian), &((*sot)[i_task].desired),jointVelocity,&PS);

		  //if ( scaleFactors[i_task]>1)  scaleFactors[i_task]=1;

		  if (scaleFactors[i_task] > 0.0){
			  if (scaleFactors[i_task]*scaleMargin < (1.0 )){
				  double taskScale=scaleFactors[i_task]*scaleMargin;
				  VectorD scaledTask=((*sot)[i_task].desired)*taskScale;
				  scaleFactors[i_task]=FOSNSsingle(i_task,&higherPriorityJointVelocity,&higherPriorityNull, &((*sot)[i_task].jacobian), &scaledTask,jointVelocity,&P);
				  scaleFactors[i_task]=taskScale;

			  }else{
				  scaleFactors[i_task]=1.0;
				  P=PS;
		 	}
		  }
	  }




/*
	  Scalar error0=( ((*sot)[0].jacobian)* (dotQopt[0]) - scaleFactors[0]*((*sot)[0].desired)  ).norm();
	  Scalar error1=( ((*sot)[0].jacobian)* (dotQopt[1]) - scaleFactors[0]*((*sot)[0].desired)  ).norm();
	  Scalar error2=( ((*sot)[0].jacobian)* (dotQopt[2]) - scaleFactors[0]*((*sot)[0].desired)  ).norm();
	  Scalar error3=( ((*sot)[0].jacobian)* (dotQopt[3]) - scaleFactors[0]*((*sot)[0].desired)  ).norm();

	  stringstream log;
	  log<<"errors "<<error0<<" "<<error1<<" "<<error2<<" "<<error3<<" \n";
	  log<<"dq0 "<<	dotQopt[0].transpose()<<endl;
	  log<<"dq1 "<<	dotQopt[1].transpose()<<endl;
	  log<<"dq2 "<<	dotQopt[2].transpose()<<endl;
	  log<<"dq3 "<<	dotQopt[3].transpose()<<endl;
	  ROS_INFO("%s",log.str().c_str());
*/
//	    ROS_INFO("s0 %f, s1 %f",scaleFactors[0],scaleFactors[1]);
/*

	  if (( ((*sot)[1].jacobian)* (*jointVelocity) - scaleFactors[1]*((*sot)[1].desired)  ).norm() > 1e-12){
		  ROS_INFO("SNS error high\n scale factors: %f , %f\n errors : %f , %f",scaleFactors[0],scaleFactors[1],( ((*sot)[0].jacobian)* (*jointVelocity) - scaleFactors[0]*((*sot)[0].desired)  ).norm(),( ((*sot)[1].jacobian)* (*jointVelocity) - scaleFactors[1]*((*sot)[1].desired)  ).norm()  );
	  }
*/

	//  if (( ((*sot)[0].jacobian)* (*jointVelocity) - scaleFactors[0]*((*sot)[0].desired)  ).norm() > 1e-8){
	//		  ROS_INFO("SNS error high\n scale factors: %f  %f\n errors : %f  \n nSat %d",scaleFactors[0],scaleFactors[0],( ((*sot)[0].jacobian)* (*jointVelocity) - scaleFactors[0]*((*sot)[0].desired)  ).norm(),nSat[0]  );
			  //Scalar error0=( ((*sot)[0].jacobian)* (dotQopt[0]) - scaleFactors[0]*((*sot)[0].desired)  ).norm();
			 // Scalar error1=( ((*sot)[0].jacobian)* (dotQopt[1]) - scaleFactors[0]*((*sot)[0].desired)  ).norm();
			  //Scalar error2=( ((*sot)[0].jacobian)* (dotQopt[2]) - scaleFactors[0]*((*sot)[0].desired)  ).norm();
			  //Scalar error3=( ((*sot)[0].jacobian)* (dotQopt[3]) - scaleFactors[0]*((*sot)[0].desired)  ).norm();

			  //stringstream log;
			 // log<<"errors "<<error0<<" "<<error1<<" "<<error2<<" "<<error3<<" \n";
			  //ROS_INFO("%s",log.str().c_str());

			 // exit(1);
		  //}


//ROS_WARN("%f - %f",scaleFactors[0],scaleFactors[1]);
	  double nSatTot=0.0;
	  for (int i=0;i<n_task;i++) nSatTot+=nSat[i];
	  return nSatTot;

}

Scalar IKL::getJointVelocity_QP(VectorD *jointVelocity, StackOfTasks *sot, VectorD *jointConfiguration){



	  return (double)1.0*nSat[0]+nSat[1]+nSat[2]+nSat[3]+nSat[4];

}


