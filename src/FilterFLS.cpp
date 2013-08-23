/*!
* @file 	FilterFLS.cpp
* @author 	Michael Blösch
* @date		10.10.2012
 */

#include "FilterFLS.hpp"
#include "Manager.hpp"
#include "tinyxml.h"
#include <map>
#include <Eigen/Cholesky>
#include <iostream>

using namespace std;

namespace LSE {

FilterFLS::FilterFLS(Manager* pManager,const char* pFilename): FilterBase(){
	pManager_ = pManager;

	// Init all parameters
	xInit_.t_ = 0;
	xInit_.x_.r_.setZero();
	xInit_.x_.v_.setZero();
	xInit_.x_.q_ = Rotations::quatIdentity();
	xInit_.x_.bf_.setZero();
	xInit_.x_.bw_.setZero();
	xInit_.x_.nr_.setZero();
	xInit_.x_.w_.setZero();
	xInit_.x_.f_ = -pManager_->g_;
	xInit_.x_.nbf_.setZero();
	xInit_.x_.nbw_.setZero();
	for(int i=0;i<LSE_N_LEG;i++){
		xInit_.CFC_[i] = 0;
		xInit_.LegArray_[i] = 0;
	}
	for(int i=0;i<1+4*(LSE_VUKF_N);i++){
		xInit_.X_[i] = xInit_.x_;
	}
	xInit_.mbSigmaSampled_ = false;
	xInit_.P_.setIdentity();
	xInit_.y_.setZero();
	Wr_ = 0*Eigen::Matrix3d::Identity();
	Wbf_ = 0.0001*Eigen::Matrix3d::Identity();
	Wbw_ = 0.000618*Eigen::Matrix3d::Identity();
	kinOutTh_ = 7.82;
	restorationFactor_ = 1;

	// Flags
	mbEstimateRotBias_ = true;
	mbEstimateAccBias_ = true;
	mbUseKin_ = true;

	// Time stepping
	mbFixedTimeStepping_ = false;
	timeStep_ = 0.0;

	// UKF Parameters
	UKFAlpha_ = 1e-3;
	UKFKappa_ = 0;
	UKFBeta_ = 2;

	loadParam(pFilename);

	// Compute derived UKF parameters (weights and lambda)
	int L = 2*LSE_VUKF_N+3*LSE_N_LEG;
	UKFLambda_ = UKFAlpha_*UKFAlpha_*(L+UKFKappa_)-L;
	UKFGamma_ = sqrt(UKFLambda_ + L);
	UKFWs_ = UKFLambda_/(L+UKFLambda_);
	UKFWc_ = UKFLambda_/(L+UKFLambda_)+(1-UKFAlpha_*UKFAlpha_+UKFBeta_);
	UKFWi_ = 1/(2*(L+UKFLambda_));

	resetEstimate(0);
}

FilterFLS::~FilterFLS(){
}

void FilterFLS::update(const double& t){
	// Find safe time and filter safe state
	double tsNew = t;
	if(!pManager_->imuMeasList_.empty() && !pManager_->encMeasList_.empty()){
		tsNew = std::min(tsNew,pManager_->imuMeasList_.rbegin()->first+pManager_->tImu_);
		tsNew = std::min(tsNew,pManager_->encMeasList_.rbegin()->first+pManager_->tEnc_);
		if(xs_.t_<tsNew){
			filterState(xs_,tsNew);
			if(pManager_->isLogging_){
				logState();
			}
		}
	}

	// Predict state
	xp_ = xs_;
	filterState(xp_,t);
}

void FilterFLS::update(){
	double tmax = 0;
	bool gotMeas = false;
	if(!pManager_->imuMeasList_.empty()){
		if(gotMeas){
			tmax = std::max(tmax,pManager_->imuMeasList_.rbegin()->first+pManager_->tImu_);
		} else {
			tmax = pManager_->imuMeasList_.rbegin()->first+pManager_->tImu_;
		}
		gotMeas = true;
	}
	if(!pManager_->encMeasList_.empty()){
		if(gotMeas){
			tmax = std::max(tmax,pManager_->encMeasList_.rbegin()->first+pManager_->tImu_);
		} else {
			tmax = pManager_->encMeasList_.rbegin()->first+pManager_->tImu_;
		}
		gotMeas = true;
	}

	if(gotMeas){
		update(tmax);
	}
}

State FilterFLS::getEst(){
	Eigen::Matrix3d R_WI,R_IB;
	R_WI = Rotations::quatToRotMat(xp_.x_.q_);
	R_IB = Rotations::quatToRotMat(pManager_->q_IB_);
	State x = State();
	x.t_ = xp_.t_;
	x.r_ = R_WI*(-xp_.x_.r_-R_IB*pManager_->B_r_BI_);
	x.v_ = R_WI*(-xp_.x_.v_-Rotations::vecToSqew(xp_.x_.w_-xp_.x_.bw_)*R_IB*pManager_->B_r_BI_);
	x.q_ = Rotations::quatInverse(Rotations::quatL(xp_.x_.q_)*pManager_->q_IB_);
	x.w_ = R_IB.transpose()*(xp_.x_.w_-xp_.x_.bw_);
	x.P_.setZero();
	x.P_.block(0,0,9,9) = xp_.P_.block(0,0,9,9);
	x.P_.block(9,9,3,3) = xp_.P_.block(12,12,3,3)+pManager_->Rw_;
	return x;
}

SlippageDetection FilterFLS::getSlippage(){
	SlippageDetection x = SlippageDetection();
	return x;
}

void FilterFLS::resetEstimate(const double& t){
	xs_ = xInit_;
	xs_.t_ = t;
	xp_ = xs_;
}

void FilterFLS::filterState(InternState& x,const double& tEnd){
	std::map<double,ImuMeas>::iterator itImu;
	std::map<double,EncMeas>::iterator itEnc;
	ImuMeas imuMeas;
	double tNext;

	// Get current measurements
	itImu = pManager_->imuMeasList_.upper_bound(x.t_-pManager_->tImu_);
	itEnc = pManager_->encMeasList_.upper_bound(x.t_-pManager_->tEnc_);

	// Get current IMU measurement
	if(itImu != pManager_->imuMeasList_.end()){
		imuMeas = itImu->second;
	} else if(!pManager_->imuMeasList_.empty()){
		itImu--;
		imuMeas = itImu->second;
		itImu++;
	} else {
		imuMeas.f_ = -pManager_->g_;
		imuMeas.w_.setZero();
	}

	while(x.t_<tEnd){
		// Get next measurement time
		tNext = tEnd;
		if(itImu != pManager_->imuMeasList_.end()) tNext = std::min(tEnd,itImu->first+pManager_->tImu_);
		if(itEnc != pManager_->encMeasList_.end()) tNext = std::min(tNext,itEnc->first+pManager_->tEnc_);

		// Predict state
		predictState(x,tNext,imuMeas);

		// Correct state if necessary
		if(itEnc != pManager_->encMeasList_.end() && tNext >= itEnc->first+pManager_->tEnc_){
			if(mbUseKin_){
				encUpdateState(x,itEnc->second);
			}
			itEnc++;
		}

		// Get next IMU measuerement if necessary
		if(itImu != pManager_->imuMeasList_.end() && tNext >= itImu->first+pManager_->tImu_){
			itImu++;
			if(itImu != pManager_->imuMeasList_.end()){
				imuMeas = itImu->second;
			}
		}
	}
}


void FilterFLS::predictState(InternState& x, const double& tPre, const ImuMeas& m){
	double dt = tPre-x.t_;
	if(mbFixedTimeStepping_){
		dt = timeStep_;
	}
	if(!mbEstimateAccBias_) x.x_.bf_.setZero();
	if(!mbEstimateRotBias_) x.x_.bw_.setZero();
	x.x_.w_ = m.w_;
	x.x_.f_ = m.f_;

	Matrix30d PA;
	PA.setZero();
	PA.block(0,0,15,15) = x.P_;
	PA.block(15,15,3,3) = Wr_/dt;
	PA.block(18,18,3,3) = pManager_->Rf_/dt;
	PA.block(21,21,3,3) = pManager_->Rw_/dt;
	PA.block(24,24,3,3) = Wbf_/dt;
	PA.block(27,27,3,3) = Wbw_/dt;

	// Cholesky decomposition
	Eigen::LLT<Matrix30d> lltOfPA(PA);
	Matrix30d SPA = lltOfPA.matrixL();


//	Eigen::Matrix<double,LSE_VUKF_N,LSE_VUKF_N> F;
//	F.setZero();
//	F.block(0,0,3,3) = -Rotations::vecToSqew(x.x_.w_-x.x_.bw_);
//	F.block(0,3,3,3) = Eigen::Matrix3d::Identity();
//	F.block(0,12,3,3) = -Rotations::vecToSqew(x.x_.r_);
//	F.block(3,3,3,3) = -Rotations::vecToSqew(x.x_.w_-x.x_.bw_);
//	F.block(3,6,3,3) = -Rotations::quatToRotMat(x.x_.q_)*Rotations::vecToSqew(pManager_->g_);
//	F.block(3,9,3,3) = Eigen::Matrix3d::Identity();
//	F.block(3,12,3,3) = -Rotations::vecToSqew(x.x_.v_);
//	F.block(6,12,3,3) = -Rotations::quatToRotMat(x.x_.q_);
//
//	F = Eigen::Matrix<double,LSE_VUKF_N,LSE_VUKF_N>::Identity() + dt*F;

	// Sample sigma points for prediction
	x.X_[0] = x.x_;
	for(int i=1;i<=2*LSE_VUKF_N;i++){
		x.X_[i].r_ = x.x_.r_ + UKFGamma_*SPA.block(0,i-1,3,1);
		x.X_[i].v_ = x.x_.v_ + UKFGamma_*SPA.block(3,i-1,3,1);
		x.X_[i].q_ = Rotations::quatL(Rotations::rotVecToQuat(UKFGamma_*SPA.block(6,i-1,3,1)))*x.x_.q_;
		if(!mbEstimateAccBias_){
			x.X_[i].bf_ = x.x_.bf_;
		} else {
			x.X_[i].bf_ = x.x_.bf_ + UKFGamma_*SPA.block(9,i-1,3,1);
		}
		if(!mbEstimateRotBias_){
			x.X_[i].bw_ = x.x_.bw_;
		} else {
			x.X_[i].bw_ = x.x_.bw_ + UKFGamma_*SPA.block(12,i-1,3,1);
		}
		x.X_[i].nr_ = UKFGamma_*SPA.block(15,i-1,3,1);
		x.X_[i].f_ = x.x_.f_ + UKFGamma_*SPA.block(18,i-1,3,1);
		x.X_[i].w_ = x.x_.w_ + UKFGamma_*SPA.block(21,i-1,3,1);
		x.X_[i].nbf_ = UKFGamma_*SPA.block(24,i-1,3,1);
		x.X_[i].nbw_ = UKFGamma_*SPA.block(27,i-1,3,1);

		x.X_[i+2*LSE_VUKF_N].r_ = x.x_.r_ - UKFGamma_*SPA.block(0,i-1,3,1);
		x.X_[i+2*LSE_VUKF_N].v_ = x.x_.v_ - UKFGamma_*SPA.block(3,i-1,3,1);
		x.X_[i+2*LSE_VUKF_N].q_ = Rotations::quatL(Rotations::rotVecToQuat(-UKFGamma_*SPA.block(6,i-1,3,1)))*x.x_.q_;
		if(!mbEstimateAccBias_){
			x.X_[i+2*LSE_VUKF_N].bf_ = x.x_.bf_;
		} else {
			x.X_[i+2*LSE_VUKF_N].bf_ = x.x_.bf_ - UKFGamma_*SPA.block(9,i-1,3,1);
		}
		if(!mbEstimateRotBias_){
			x.X_[i+2*LSE_VUKF_N].bw_ = x.x_.bw_;
		} else {
			x.X_[i+2*LSE_VUKF_N].bw_ = x.x_.bw_ - UKFGamma_*SPA.block(12,i-1,3,1);
		}
		x.X_[i+2*LSE_VUKF_N].nr_ = -UKFGamma_*SPA.block(15,i-1,3,1);
		x.X_[i+2*LSE_VUKF_N].f_ = x.x_.f_ - UKFGamma_*SPA.block(18,i-1,3,1);
		x.X_[i+2*LSE_VUKF_N].w_ = x.x_.w_ - UKFGamma_*SPA.block(21,i-1,3,1);
		x.X_[i+2*LSE_VUKF_N].nbf_ = -UKFGamma_*SPA.block(24,i-1,3,1);
		x.X_[i+2*LSE_VUKF_N].nbw_ = -UKFGamma_*SPA.block(27,i-1,3,1);
	}

	// Propagate sigma points
	AugmentedState as = x.x_;
	Eigen::Matrix3d G0T = pManager_->gamma(0,-dt*(as.w_-as.bw_));
	Eigen::Matrix3d G1T = pManager_->gamma(1,-dt*(as.w_-as.bw_));
	Eigen::Matrix3d G2T = pManager_->gamma(2,-dt*(as.w_-as.bw_));
	Eigen::Matrix3d R_IW, R_WI;
	bool upGyroStuff = true;
	for(int i=0;i<=2*(2*LSE_VUKF_N);i++){
		if (as.w_-as.bw_ == x.X_[i].w_ - x.X_[i].bw_) upGyroStuff = false;
		as = x.X_[i];
		if(upGyroStuff){
			G0T = pManager_->gamma(0,-dt*(as.w_-as.bw_));
			G1T = pManager_->gamma(1,-dt*(as.w_-as.bw_));
			G2T = pManager_->gamma(2,-dt*(as.w_-as.bw_));
		}
		R_WI = Rotations::quatToRotMat(as.q_);
		R_IW = R_WI.transpose();

		as.r_ = G0T*(as.r_+dt*as.v_-dt*dt/2*(2*G2T*(as.f_-as.bf_)+R_IW*pManager_->g_))+as.nr_*dt;
		as.v_ = G0T*(as.v_-dt*(G1T*(as.f_-as.bf_)+R_IW*pManager_->g_));
		as.q_ = Rotations::quatL(as.q_)*Rotations::rotVecToQuat(dt*(as.w_-as.bw_));
		as.bf_ = as.bf_ + as.nbf_*dt;
		as.bw_ = as.bw_ + as.nbw_*dt;
		x.X_[i] = as;
		bool upGyroStuff = true;
	}

	// Compute predicted state and covariance
	Eigen::Vector3d vec;
	vec.setZero();
	x.x_.r_ = (UKFWs_+2*3*LSE_N_LEG*UKFWi_)*x.X_[0].r_;
	x.x_.v_ = (UKFWs_+2*3*LSE_N_LEG*UKFWi_)*x.X_[0].v_;
	x.x_.bf_ = (UKFWs_+2*3*LSE_N_LEG*UKFWi_)*x.X_[0].bf_;
	x.x_.bw_ = (UKFWs_+2*3*LSE_N_LEG*UKFWi_)*x.X_[0].bw_;
	for(int i=1;i<=2*(2*LSE_VUKF_N);i++){
		x.x_.r_ += UKFWi_*x.X_[i].r_;
		x.x_.v_ += UKFWi_*x.X_[i].v_;
		vec += UKFWi_*Rotations::quatToRotVec(Rotations::quatL(x.X_[i].q_)*Rotations::quatInverse(x.X_[0].q_));
		x.x_.bf_ += UKFWi_*x.X_[i].bf_;
		x.x_.bw_ += UKFWi_*x.X_[i].bw_;
	}
	x.x_.q_ = Rotations::quatL(Rotations::rotVecToQuat(vec))*x.X_[0].q_;
	Eigen::Matrix<double,LSE_VUKF_N,1> vec15;
	vec15.block(0,0,3,1) = x.X_[0].r_-x.x_.r_ ;
	vec15.block(3,0,3,1) = x.X_[0].v_-x.x_.v_ ;
	vec15.block(6,0,3,1) = Rotations::quatToRotVec(Rotations::quatL(x.X_[0].q_)*Rotations::quatInverse(x.x_.q_));
	vec15.block(9,0,3,1) = x.X_[0].bf_-x.x_.bf_ ;
	vec15.block(12,0,3,1) = x.X_[0].bw_-x.x_.bw_ ;
	x.P_ = (UKFWc_+2*3*LSE_N_LEG*UKFWi_)*vec15*vec15.transpose();
	for(int i=1;i<=2*(2*LSE_VUKF_N);i++){
		vec15.block(0,0,3,1) = x.X_[i].r_-x.x_.r_ ;
		vec15.block(3,0,3,1) = x.X_[i].v_-x.x_.v_ ;
		vec15.block(6,0,3,1) = Rotations::quatToRotVec(Rotations::quatL(x.X_[i].q_)*Rotations::quatInverse(x.x_.q_));
		vec15.block(9,0,3,1) = x.X_[i].bf_-x.x_.bf_ ;
		vec15.block(12,0,3,1) = x.X_[i].bw_-x.x_.bw_ ;
		x.P_ += UKFWi_*vec15*vec15.transpose();
	}

//	//EKF
//	x.x_ = x.X_[0];
//	x.P_ = F*x.P_*F.transpose() + PA.block(15,15,15,15)*dt*dt;


	// Avoid singular P
	if(!mbEstimateAccBias_) x.P_.block(9,9,3,3).setIdentity();
	if(!mbEstimateRotBias_) x.P_.block(12,12,3,3).setIdentity();

	x.t_ = tPre;
	x.mbSigmaSampled_ = true;
}

void FilterFLS::encUpdateState(InternState& x, const EncMeas& m){
	if(x.mbSigmaSampled_){
		// Update Contact count
		for(int i=0;i<LSE_N_LEG;i++){
			if(m.CF_[i]){
				x.CFC_[i]++;
			} else {
				x.CFC_[i] = 0;
			}
		}

		// Determine contact foots, reject outliers
		for(int i=0;i<LSE_N_LEG;i++){
			if(x.CFC_[i]>0){
				x.LegArray_[i] = 1;
			} else {
				x.LegArray_[i] = 0;
			}
		}

		// Compute forward kinematics and Jacobians
		Eigen::Matrix<double,3*LSE_N_LEG,1> s;
		Eigen::Matrix<double,3*LSE_N_LEG,LSE_DOF_LEG> J;
		for(int i=0;i<LSE_N_LEG;i++){
			// I_r_IF = C(q_IB)*(-B_r_BI + B_r_BK + C'(q_KB)*K_r_KF
			s.block(3*i,0,3,1) = Rotations::quatToRotMat(pManager_->q_IB_)*(-pManager_->B_r_BI_+pManager_->B_r_BK_+Rotations::quatToRotMat(pManager_->q_KB_).transpose()*(*pManager_->legKin)(m.e_.col(i),i));
			J.block(3*i,0,3,LSE_DOF_LEG) = Rotations::quatToRotMat(pManager_->q_IB_)*Rotations::quatToRotMat(pManager_->q_KB_).transpose()*(*pManager_->legKinJac)(m.e_.col(i),i);
		}

		// Compute measurement noise sigma points using Cholesky decomposition
		Eigen::LLT<Eigen::Matrix3d> lltOfRda(pManager_->Rda_);
		Eigen::Matrix3d SRda = lltOfRda.matrixL();
		SRda = UKFGamma_*SRda;

		// Project through measurement function
		AugmentedState as = x.x_;
		Eigen::Matrix<double,3*LSE_N_LEG,1> nf;
		Eigen::Matrix<double,3*LSE_N_LEG,2*(2*LSE_VUKF_N+3*LSE_N_LEG)+1> Y;
		int noiColInd = 0;
		for(int i=0;i<=2*(2*LSE_VUKF_N+3*LSE_N_LEG);i++){
			if(i<=2*(2*LSE_VUKF_N)){
				as = x.X_[i];
				nf.setZero();
			} else if (i<=2*(2*LSE_VUKF_N)+3*LSE_N_LEG) {
				as = x.X_[0];
				nf.setZero();
				noiColInd = std::floor((i-2*(2*LSE_VUKF_N)-1)/3);
				nf.block(3*noiColInd,0,3,1) = SRda.col((i-2*(2*LSE_VUKF_N)-1)%3);
			} else {
				as = x.X_[0];
				nf.setZero();
				noiColInd = std::floor((i-2*(2*LSE_VUKF_N)-3*LSE_N_LEG-1)/3);
				nf.block(3*noiColInd,0,3,1) = -SRda.col((i-2*(2*LSE_VUKF_N)-3*LSE_N_LEG-1)%3);
			}
			for(int j=0;j<LSE_N_LEG;j++){
				Y.block(j*3,i,3,1) = -as.v_ + Rotations::vecToSqew(as.w_-as.bw_)*s.block(3*j,0,3,1) + J.block(3*j,0,3,LSE_DOF_LEG)*m.v_.col(j) - nf.block(3*j,0,3,1);
			}
		}

		// Update: compute innovation, corresponding covariance, cross-covariance, Kalman gain
		x.y_ = UKFWs_*Y.col(0);
		for(int i=1;i<=2*(2*LSE_VUKF_N+3*LSE_N_LEG);i++){
			x.y_ += UKFWi_*Y.col(i);
		}

		Eigen::Matrix<double,3*LSE_N_LEG,3*LSE_N_LEG> Py = UKFWc_*(Y.col(0)-x.y_)*(Y.col(0)-x.y_).transpose();
		for(int i=1;i<=2*(2*LSE_VUKF_N+3*LSE_N_LEG);i++){
			Py += UKFWi_*(Y.col(i)-x.y_)*(Y.col(i)-x.y_).transpose();
		}

		Eigen::Matrix<double,LSE_VUKF_N,1> vec15;
		vec15.block(0,0,3,1) = x.X_[0].r_-x.x_.r_ ;
		vec15.block(3,0,3,1) = x.X_[0].v_-x.x_.v_ ;
		vec15.block(6,0,3,1) = Rotations::quatToRotVec(Rotations::quatL(x.X_[0].q_)*Rotations::quatInverse(x.x_.q_));
		vec15.block(9,0,3,1) = x.X_[0].bf_-x.x_.bf_ ;
		vec15.block(12,0,3,1) = x.X_[0].bw_-x.x_.bw_ ;
		Eigen::Matrix<double,LSE_VUKF_N,3*LSE_N_LEG> Pxy = UKFWc_*vec15*(Y.col(0)-x.y_).transpose();
		for(int i=1;i<=2*(2*LSE_VUKF_N+3*LSE_N_LEG);i++){
			if(i<=2*(2*LSE_VUKF_N)){
				as = x.X_[i];
			} else {
				as = x.X_[0];
			}
			vec15.block(0,0,3,1) = as.r_-x.x_.r_ ;
			vec15.block(3,0,3,1) = as.v_-x.x_.v_ ;
			vec15.block(6,0,3,1) = Rotations::quatToRotVec(Rotations::quatL(as.q_)*Rotations::quatInverse(x.x_.q_));
			vec15.block(9,0,3,1) = as.bf_-x.x_.bf_ ;
			vec15.block(12,0,3,1) = as.bw_-x.x_.bw_ ;
			Pxy += UKFWi_*vec15*(Y.col(i)-x.y_).transpose();
		}

		// Compute inverse of innovation covariance and reject outliers (the probability to find y out of the 3-sigma bound is about 0.25%
		Eigen::Matrix<double,3*LSE_N_LEG,3*LSE_N_LEG> Pyinv = Py.inverse();
		outlierDetection(x,Pyinv);
		for(int i=0;i<LSE_N_LEG;i++){
			if(x.LegArray_[i]==0){
				Pyinv.block(3*i,0,3,3*LSE_N_LEG).setZero();
				Pyinv.block(0,3*i,3*LSE_N_LEG,3).setZero();
			}
		}

		Eigen::Matrix<double,LSE_VUKF_N,3*LSE_N_LEG> K = Pxy*Pyinv;


//		// EKF
//		y = Y.col(0);
//		Eigen::Matrix<double,3*LSE_N_LEG,LSE_VUKF_N> G;
//		G.setZero();
//		for(int i=0;i<LSE_N_LEG;i++){
//			G.block(3*i,3,3,3) = -Eigen::Matrix3d::Identity();
//			G.block(3*i,12,3,3) = Rotations::vecToSqew(s.block(3*i,0,3,1));
//		}
//		Py = G*x.P_*G.transpose();
//		for(int i=0;i<LSE_N_LEG;i++){
//			Py.block(3*i,3*i,3,3) += pManager_->Rda_;
//		}
//		Pyinv = Py.inverse();
//		for(int i=0;i<LSE_N_LEG;i++){
//			if(x.LegArray_[i]==0){
//				Pyinv.block(3*i,0,3,3*LSE_N_LEG).setZero();
//				Pyinv.block(0,3*i,3*LSE_N_LEG,3).setZero();
//			}
//		}
//		K = x.P_*G.transpose()*Pyinv;

		// Update state and covariance matrix
		vec15 = -K*x.y_;
		x.x_.r_ = x.x_.r_ + vec15.block(0,0,3,1);
		x.x_.v_ = x.x_.v_ + vec15.block(3,0,3,1);
		x.x_.q_ = Rotations::quatL(Rotations::rotVecToQuat(vec15.block(6,0,3,1)))*x.x_.q_;
		x.x_.bf_ = x.x_.bf_ + vec15.block(9,0,3,1);
		x.x_.bw_ = x.x_.bw_ + vec15.block(12,0,3,1);
		x.P_ = x.P_ - K*Py*K.transpose();

		// TODO check unobservability
	}
}

void FilterFLS::outlierDetection(InternState& x,const Eigen::Matrix<double,12,12>& Pyinv){
	bool outliers[LSE_N_LEG];
	double ratio[LSE_N_LEG];

	// Compute ratios (weighted errors) and check if any contact available (for a potential outlier restoration)
	bool restoreOutliers = false;
	for(int i=0;i<LSE_N_LEG;i++){
		ratio[i] = (x.y_.block(i*3,0,3,1).transpose()*Pyinv.block(3*i,3*i,3,3)*x.y_.block(i*3,0,3,1))(0,0);
		if(x.LegArray_[i]!=0) restoreOutliers = true;
	}

	// Initial outlier detection based on predicted measurement covariance
	for(int i=0;i<LSE_N_LEG;i++){
		outliers[i] = 0;
		if(x.LegArray_[i]!=0){
			if(ratio[i]>(kinOutTh_)){
				outliers[i] = 1;
			} else {
				restoreOutliers = false;
			}
		}
	}

	// Rather keep an outlier in than having no measurement at all
	if(restoreOutliers){
		double minRatio = 0;
		for(int i=0;i<LSE_N_LEG;i++){
			if(x.LegArray_[i]!=0){
				if(ratio[i] < minRatio || minRatio == 0){
					minRatio = ratio[i];
				}
			}
		}

		// Recheck with higher outlier threshold such that at least one measurement is kept
		for(int i=0;i<LSE_N_LEG;i++){
			outliers[i] = 0;
			if(x.LegArray_[i]!=0){
				if(ratio[i]>minRatio*restorationFactor_){
					outliers[i] = 1;
				}
			}
		}
	}



	// Store in state x
	for(int i=0;i<LSE_N_LEG;i++){
		if(outliers[i]){
			x.LegArray_[i]=0;
//			std::cout << "Detected Outlier" << std::endl;
		}
	}
}

void FilterFLS::loadParam(const char* pFilename){
	// Open parameter file
	TiXmlDocument doc(pFilename);
	if (!doc.LoadFile()) return;

	// Define handles and elements
	TiXmlHandle hDoc(&doc);
	TiXmlElement* pElem;
	TiXmlHandle hRoot(0);

	int mInt;

	// Get root
	pElem=hDoc.FirstChildElement("LeggedStateEstimator").Element();
	if (pElem){
		hRoot=TiXmlHandle(pElem);

		// MeasurementsSettings
		pElem=hRoot.FirstChild("VUKFSettings").FirstChild("Position").FirstChild("Init").Element();
		if (pElem){
			pElem->QueryDoubleAttribute("x", &xInit_.x_.r_(0));
			pElem->QueryDoubleAttribute("y", &xInit_.x_.r_(1));
			pElem->QueryDoubleAttribute("z", &xInit_.x_.r_(2));
		}
		pElem=hRoot.FirstChild("VUKFSettings").FirstChild("Position").FirstChild("InitStd").Element();
		if (pElem){
			pElem->QueryDoubleAttribute("x", &xInit_.P_(0,0));
			pElem->QueryDoubleAttribute("y", &xInit_.P_(1,1));
			pElem->QueryDoubleAttribute("z", &xInit_.P_(2,2));
		}

		pElem=hRoot.FirstChild("VUKFSettings").FirstChild("Velocity").FirstChild("Init").Element();
		if (pElem){
			pElem->QueryDoubleAttribute("x", &xInit_.x_.v_(0));
			pElem->QueryDoubleAttribute("y", &xInit_.x_.v_(1));
			pElem->QueryDoubleAttribute("z", &xInit_.x_.v_(2));
		}
		pElem=hRoot.FirstChild("VUKFSettings").FirstChild("Velocity").FirstChild("InitStd").Element();
		if (pElem){
			pElem->QueryDoubleAttribute("x", &xInit_.P_(3,3));
			pElem->QueryDoubleAttribute("y", &xInit_.P_(4,4));
			pElem->QueryDoubleAttribute("z", &xInit_.P_(5,5));
		}

		pElem=hRoot.FirstChild("VUKFSettings").FirstChild("Attitude").FirstChild("Init").Element();
		if (pElem){
			pElem->QueryDoubleAttribute("x", &xInit_.x_.q_(0));
			pElem->QueryDoubleAttribute("y", &xInit_.x_.q_(1));
			pElem->QueryDoubleAttribute("z", &xInit_.x_.q_(2));
			pElem->QueryDoubleAttribute("w", &xInit_.x_.q_(3));
		}
		pElem=hRoot.FirstChild("VUKFSettings").FirstChild("Attitude").FirstChild("InitStd").Element();
		if (pElem){
			pElem->QueryDoubleAttribute("x", &xInit_.P_(6,6));
			pElem->QueryDoubleAttribute("y", &xInit_.P_(7,7));
			pElem->QueryDoubleAttribute("z", &xInit_.P_(8,8));
		}

		pElem=hRoot.FirstChild("VUKFSettings").FirstChild("AccelerometerBias").FirstChild("Init").Element();
		if (pElem){
			pElem->QueryDoubleAttribute("x", &xInit_.x_.bf_(0));
			pElem->QueryDoubleAttribute("y", &xInit_.x_.bf_(1));
			pElem->QueryDoubleAttribute("z", &xInit_.x_.bf_(2));
		}
		pElem=hRoot.FirstChild("VUKFSettings").FirstChild("AccelerometerBias").FirstChild("InitStd").Element();
		if (pElem){
			pElem->QueryDoubleAttribute("x", &xInit_.P_(9,9));
			pElem->QueryDoubleAttribute("y", &xInit_.P_(10,10));
			pElem->QueryDoubleAttribute("z", &xInit_.P_(11,11));
		}

		pElem=hRoot.FirstChild("VUKFSettings").FirstChild("GyroscopeBias").FirstChild("Init").Element();
		if (pElem){
			pElem->QueryDoubleAttribute("x", &xInit_.x_.bw_(0));
			pElem->QueryDoubleAttribute("y", &xInit_.x_.bw_(1));
			pElem->QueryDoubleAttribute("z", &xInit_.x_.bw_(2));
		}
		pElem=hRoot.FirstChild("VUKFSettings").FirstChild("GyroscopeBias").FirstChild("InitStd").Element();
		if (pElem){
			pElem->QueryDoubleAttribute("x", &xInit_.P_(12,12));
			pElem->QueryDoubleAttribute("y", &xInit_.P_(13,13));
			pElem->QueryDoubleAttribute("z", &xInit_.P_(14,14));
		}
		pElem=hRoot.FirstChild("VUKFSettings").FirstChild("Position").FirstChild("PreNoiStd").Element();
		if (pElem){
			pElem->QueryDoubleAttribute("x", &Wr_(0,0));
			pElem->QueryDoubleAttribute("y", &Wr_(1,1));
			pElem->QueryDoubleAttribute("z", &Wr_(2,2));
		}
		pElem=hRoot.FirstChild("VUKFSettings").FirstChild("AccelerometerBias").FirstChild("PreNoiStd").Element();
		if (pElem){
			pElem->QueryDoubleAttribute("x", &Wbf_(0,0));
			pElem->QueryDoubleAttribute("y", &Wbf_(1,1));
			pElem->QueryDoubleAttribute("z", &Wbf_(2,2));
		}
		pElem=hRoot.FirstChild("VUKFSettings").FirstChild("GyroscopeBias").FirstChild("PreNoiStd").Element();
		if (pElem){
			pElem->QueryDoubleAttribute("x", &Wbw_(0,0));
			pElem->QueryDoubleAttribute("y", &Wbw_(1,1));
			pElem->QueryDoubleAttribute("z", &Wbw_(2,2));
		}
		pElem=hRoot.FirstChild("VUKFSettings").Element();
		if (pElem){
			pElem->QueryDoubleAttribute("alpha", &UKFAlpha_);
			pElem->QueryDoubleAttribute("beta", &UKFBeta_);
			pElem->QueryDoubleAttribute("kappa", &UKFKappa_);
			pElem->QueryDoubleAttribute("timeStepping", &timeStep_);
			pElem->QueryIntAttribute("useKin", &mInt);
			mbUseKin_ = mInt;
		}
		pElem=hRoot.FirstChild("VUKFSettings").FirstChild("AccelerometerBias").Element();
		if (pElem){
			pElem->QueryIntAttribute("estimate", &mInt);
			mbEstimateAccBias_ = mInt;
		}
		pElem=hRoot.FirstChild("VUKFSettings").FirstChild("GyroscopeBias").Element();
		if (pElem){
			pElem->QueryIntAttribute("estimate", &mInt);
			mbEstimateRotBias_ = mInt;
		}
		pElem=hRoot.FirstChild("VUKFSettings").FirstChild("Foothold").Element();
		if (pElem){
			pElem->QueryDoubleAttribute("outlierThreshold", &kinOutTh_);
			pElem->QueryDoubleAttribute("restorationFactor", &restorationFactor_);
		}
	}

	xInit_.P_ = xInit_.P_*xInit_.P_;
	Wr_ = Wr_*Wr_;
	Wbf_ = Wbf_*Wbf_;
	Wbw_ = Wbw_*Wbw_;

	if(timeStep_==0){
		mbFixedTimeStepping_ = false;
	} else {
		mbFixedTimeStepping_ = true;
	}

	if(restorationFactor_<1){
		restorationFactor_ = 1;
	}
}

void FilterFLS::logState(){
	  pManager_->ofsLog_ << xs_.t_ << "\t";
	  pManager_->ofsLog_ << xs_.x_.r_(0) << "\t" << xs_.x_.r_(1) << "\t" << xs_.x_.r_(2) << "\t";
	  pManager_->ofsLog_ << xs_.x_.v_(0) << "\t" << xs_.x_.v_(1) << "\t" << xs_.x_.v_(2) << "\t";
	  pManager_->ofsLog_ << xs_.x_.q_(0) << "\t" << xs_.x_.q_(1) << "\t" << xs_.x_.q_(2) << "\t" << xs_.x_.q_(3) << "\t";
	  pManager_->ofsLog_ << xs_.x_.bf_(0) << "\t" << xs_.x_.bf_(1) << "\t" << xs_.x_.bf_(2) << "\t";
	  pManager_->ofsLog_ << xs_.x_.bw_(0) << "\t" << xs_.x_.bw_(1) << "\t" << xs_.x_.bw_(2) << "\t";
	  pManager_->ofsLog_ << xs_.LegArray_[0] << "\t" << xs_.LegArray_[1] << "\t" << xs_.LegArray_[2] << "\t" << xs_.LegArray_[3] << "\t";
	  for(int i=0;i<LSE_VUKF_N;i++){
		  pManager_->ofsLog_ << xs_.P_(i,i) << "\t";
	  }
	  for(int i=0;i<LSE_DOF_LEG*LSE_N_LEG;i++){
		  pManager_->ofsLog_ << xs_.y_(i) << "\t";
	  }
	  pManager_->ofsLog_ << endl;
}

std::string FilterFLS::getKeyString(){
	std::ostringstream oss (std::ostringstream::out);
	oss << pManager_->Rda_(0,0) << "_" << kinOutTh_ << "_" << restorationFactor_;
	return oss.str();
}

}




