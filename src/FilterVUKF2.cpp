/*!
* @file 	FilterVUKF2.cpp
* @author 	Michael Blösch
* @date		10.10.2012
 */

#include "FilterVUKF2.hpp"
#include "Manager.hpp"
#include "tinyxml.h"
#include <map>
#include <Eigen/Cholesky>
#include <iostream>

//TODO: make timelimit on safe state
//TODO: increase efficiency of PN and PU computation
//TODO: reduce to euler forward

using namespace std;

namespace LSE {

FilterVUKF2::FilterVUKF2(Manager* pManager,const char* pFilename): FilterBase(){
	pManager_ = pManager;

	// Init all parameters
	xInit_.t_ = 0;
	xInit_.x_.r_.setZero();
	xInit_.x_.v_.setZero();
	xInit_.x_.q_ = Rotations::quatIdentity();
	xInit_.x_.bf_.setZero();
	xInit_.x_.bw_.setZero();
	xInit_.w_.setZero();
	xInit_.f_ = -pManager_->g_;
	for(int i=0;i<LSE_N_LEG;i++){
		xInit_.CFC_[i] = 0;
		xInit_.slippageDetection_.flag_[i] = 0;
	}
	for(int i=0;i<1+2*(VUKFF_state_dim+VUKFF_preNoise_dim+VUKFF_upNoise_dim);i++){
		xInit_.X_[i] = xInit_.x_;
	}
	xInit_.mbSigmaSampled_ = false;
	xInit_.P_.setIdentity();
	xInit_.PN_.setZero();
	xInit_.UN_.setZero();

	// This is the square root of the final value!!!
	Wr_ = 0.0001*Eigen::Matrix3d::Identity();
	Wv_ = 0.003*Eigen::Matrix3d::Identity();
	Wq_ = 0.00001*Eigen::Matrix3d::Identity();
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
	int L = VUKFF_state_dim+VUKFF_preNoise_dim+VUKFF_upNoise_dim;
	UKFLambda_ = UKFAlpha_*UKFAlpha_*(L+UKFKappa_)-L;
	UKFGamma_ = sqrt(UKFLambda_ + L);
	UKFWs_ = UKFLambda_/(L+UKFLambda_);
	UKFWc_ = UKFLambda_/(L+UKFLambda_)+(1-UKFAlpha_*UKFAlpha_+UKFBeta_);
	UKFWi_ = 1/(2*(L+UKFLambda_));

	resetEstimate(0);
}

FilterVUKF2::~FilterVUKF2(){
}

void FilterVUKF2::update(const double& t){
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

void FilterVUKF2::update(){
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

State FilterVUKF2::getEst(){
	Eigen::Matrix3d R_WI,R_IB;
	R_WI = Rotations::quatToRotMat(xp_.x_.q_);
	R_IB = Rotations::quatToRotMat(pManager_->q_IB_);
	State x = State();
	x.t_ = xp_.t_;
	x.r_ = R_WI*(-xp_.x_.r_-R_IB*pManager_->B_r_BI_);
	x.v_ = R_WI*(-xp_.x_.v_-Rotations::vecToSqew(xp_.w_-xp_.x_.bw_)*R_IB*pManager_->B_r_BI_);
	x.q_ = Rotations::quatInverse(Rotations::quatL(xp_.x_.q_)*pManager_->q_IB_);
	x.w_ = R_IB.transpose()*(xp_.w_-xp_.x_.bw_);
	x.P_.setZero();
	x.P_.block(0,0,9,9) = xp_.P_.block(0,0,9,9);
	x.P_.block(9,9,3,3) = xp_.P_.block(12,12,3,3)+pManager_->Rw_;
	return x;
}

SlippageDetection FilterVUKF2::getSlippage(){
	SlippageDetection x = SlippageDetection();
	return x;
}

void FilterVUKF2::resetEstimate(const double& t){
	xs_ = xInit_;
	xs_.t_ = t;
	xp_ = xs_;
}

void FilterVUKF2::filterState(InternState& x,const double& tEnd){
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


void FilterVUKF2::predictState(InternState& x, const double& tPre, const ImuMeas& m){
	double dt = tPre-x.t_;
	if(mbFixedTimeStepping_){
		dt = timeStep_;
	}
	if(dt<1e-6){
		dt = 1e-6;
	}
	if(!mbEstimateAccBias_) x.x_.bf_.setZero();
	if(!mbEstimateRotBias_) x.x_.bw_.setZero();
	x.w_ = m.w_;
	x.f_ = m.f_;

	// Prediction noise covariance matrix
	samplePredictionNoise(x,dt);

	// Sample Sigma Points
	// Cholesky decomposition of covariance matrix
	Eigen::LLT<MatrixP> lltOfP(x.P_);
	if(lltOfP.info()==Eigen::NumericalIssue) std::cout << "Numerical issues while computing Cholesky of P" << std::endl;
	SP_ = lltOfP.matrixL();
	SP_ = SP_*UKFGamma_;
	x.X_[0] = x.x_;
	for(int i=1;i<=VUKFF_state_dim;i++){
		x.X_[i] = x.x_+SP_.col(i-1);
		x.X_[i+VUKFF_state_dim] = x.x_+(-1*SP_.col(i-1));
	}

	// Propagate Sigma Points
	// State Part
	for(int i=0;i<=2*VUKFF_state_dim;i++){
		predict(x.X_[i],dt,m);
	}
	Eigen::Matrix<double,VUKFF_preNoise_dim,1> n;
	// Prediction noise Part
	for(int i=1;i<=2*VUKFF_preNoise_dim;i++){
		n = x.PN_.col(i-1);
		// Handle case where Bias estimation disabled
		if(!mbEstimateAccBias_) n.block<3,1>(9,0).setZero();
		if(!mbEstimateRotBias_) n.block<3,1>(12,0).setZero();

		x.X_[2*VUKFF_state_dim+i] = x.x_;
		predict(x.X_[2*VUKFF_state_dim+i],dt,m,n);
	}

	// Compute predicted state and covariance
	Eigen::Matrix<double,VUKFF_state_dim,1> vec;
	vec.setZero();
	for(int i=1;i<=2*(VUKFF_state_dim+VUKFF_preNoise_dim);i++){
		vec = vec + (UKFWi_*(x.X_[i]-x.X_[0]));
	}
	x.x_ = x.X_[0]+vec;
	vec = x.X_[0]-x.x_;
	x.P_ = (UKFWc_+2*VUKFF_upNoise_dim*UKFWi_)*vec*vec.transpose();
	for(int i=1;i<=2*(VUKFF_state_dim+VUKFF_preNoise_dim);i++){
		vec = x.X_[i]-x.x_;
		x.P_ += UKFWi_*vec*vec.transpose();
	}

	// Avoid singular P
	if(!mbEstimateAccBias_) x.P_.block<3,3>(9,9) = xInit_.P_.block<3,3>(9,9);
	if(!mbEstimateRotBias_) x.P_.block<3,3>(12,12) = xInit_.P_.block<3,3>(12,12);

	x.t_ = tPre;
	x.mbSigmaSampled_ = true;
}

void FilterVUKF2::encUpdateState(InternState& x, const EncMeas& m){
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
				x.slippageDetection_.flag_[i] = 1;
			} else {
				x.slippageDetection_.flag_[i] = 0;
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

		// Update noise covariance matrix
		sampleUpdateNoise(x);

		// Project through measurement function
		VUKFFilterState filterState;
		Eigen::Matrix<double,VUKFF_upNoise_dim,1> upNoi;
		Eigen::Vector3d wNoise;
		Eigen::Matrix<double,VUKFF_upNoise_dim,1+2*(VUKFF_state_dim+VUKFF_preNoise_dim+VUKFF_upNoise_dim)> Y;
		for(int i=0;i<=2*(VUKFF_state_dim+VUKFF_preNoise_dim+VUKFF_upNoise_dim);i++){
			if(i<=2*(VUKFF_state_dim)){
				filterState = x.X_[i];
				wNoise.setZero();
				upNoi.setZero();
			} else if (i<=2*(VUKFF_state_dim+VUKFF_preNoise_dim)) {
				filterState = x.X_[0];
				wNoise = x.PN_.block<3,1>(18,(i-2*(VUKFF_state_dim)-1));
				upNoi.setZero();
			}  else {
				filterState = x.X_[0];
				wNoise.setZero();
				upNoi = x.UN_.col((i-2*(VUKFF_state_dim+VUKFF_preNoise_dim)-1));
			}
			for(int j=0;j<LSE_N_LEG;j++){
				Y.block(j*3,i,3,1) = -filterState.v_ + Rotations::vecToSqew(x.w_-filterState.bw_-wNoise)*s.block(3*j,0,3,1) + J.block(3*j,0,3,LSE_DOF_LEG)*m.v_.col(j)-upNoi.block<3,1>(j*3,0);
			}
		}

		// Compute innovation and corresponding covariance
		x.y_ = UKFWs_*Y.col(0);
		for(int i=1;i<=2*(VUKFF_state_dim+VUKFF_preNoise_dim+VUKFF_upNoise_dim);i++){
			x.y_ += UKFWi_*Y.col(i);
		}
		Eigen::Matrix<double,VUKFF_upNoise_dim,VUKFF_upNoise_dim> Py;
		Py = UKFWc_*(Y.col(0)-x.y_)*(Y.col(0)-x.y_).transpose();
		for(int i=1;i<=2*(VUKFF_state_dim+VUKFF_preNoise_dim+VUKFF_upNoise_dim);i++){
			Py += UKFWi_*(Y.col(i)-x.y_)*(Y.col(i)-x.y_).transpose();
		}

		// Compute cross-correlation
		Eigen::Matrix<double,VUKFF_state_dim,1> vec;
		vec = x.X_[0]-x.x_;
		Eigen::Matrix<double,VUKFF_state_dim,VUKFF_upNoise_dim> Pxy;
		Pxy = UKFWc_*vec*(Y.col(0)-x.y_).transpose();
		for(int i=1;i<=2*(VUKFF_state_dim+VUKFF_preNoise_dim+VUKFF_upNoise_dim);i++){
			if(i<=2*(VUKFF_state_dim+VUKFF_preNoise_dim)){
				filterState = x.X_[i];
			} else {
				filterState = x.X_[0];
			}
			vec = filterState-x.x_;
			Pxy += UKFWi_*vec*(Y.col(i)-x.y_).transpose();
		}

		// Compute inverse of innovation covariance
		Eigen::Matrix<double,VUKFF_upNoise_dim,VUKFF_upNoise_dim> Pyinv = Py.inverse();
		outlierDetection(x,Pyinv);
		for(int i=0;i<LSE_N_LEG;i++){
			if(x.slippageDetection_.flag_[i]==0){
				Pyinv.block<3,VUKFF_upNoise_dim>(3*i,0).setZero();
				Pyinv.block<VUKFF_upNoise_dim,3>(0,3*i).setZero();
			}
		}

		// Compute Kalman Gain
		Eigen::Matrix<double,VUKFF_state_dim,VUKFF_upNoise_dim> K = Pxy*Pyinv;

		// Update state and covariance matrix
		vec = -K*x.y_;
		x.x_ = x.x_+vec;
		x.P_ = x.P_ - K*Py*K.transpose();
	}
}

void FilterVUKF2::outlierDetection(InternState& x,const Eigen::Matrix<double,12,12>& Pyinv){
	bool outliers[LSE_N_LEG];
	double ratio[LSE_N_LEG];

	// Compute ratios (weighted errors) and check if any contact available (for a potential outlier restoration)
	bool restoreOutliers = false;
	for(int i=0;i<LSE_N_LEG;i++){
		ratio[i] = (x.y_.block(i*3,0,3,1).transpose()*Pyinv.block(3*i,3*i,3,3)*x.y_.block(i*3,0,3,1))(0,0);
		if(x.slippageDetection_.flag_[i]!=0) restoreOutliers = true;
	}

	// Initial outlier detection based on predicted measurement covariance
	for(int i=0;i<LSE_N_LEG;i++){
		outliers[i] = 0;
		if(x.slippageDetection_.flag_[i]!=0){
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
			if(x.slippageDetection_.flag_[i]!=0){
				if(ratio[i] < minRatio || minRatio == 0){
					minRatio = ratio[i];
				}
			}
		}

		// Recheck with higher outlier threshold such that at least one measurement is kept
		for(int i=0;i<LSE_N_LEG;i++){
			outliers[i] = 0;
			if(x.slippageDetection_.flag_[i]!=0){
				if(ratio[i]>minRatio*restorationFactor_){
					outliers[i] = 1;
				}
			}
		}
	}



	// Store in state x
	for(int i=0;i<LSE_N_LEG;i++){
		if(outliers[i]){
			x.slippageDetection_.flag_[i]=0;
			x.slippageDetection_.slip_.col(i)=x.y_.block<3,1>(3*i,0);
		}
	}
}

void FilterVUKF2::loadParam(const char* pFilename){
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
		pElem=hRoot.FirstChild("VUKFSettings").FirstChild("Velocity").FirstChild("PreNoiStd").Element();
		if (pElem){
			pElem->QueryDoubleAttribute("x", &Wv_(0,0));
			pElem->QueryDoubleAttribute("y", &Wv_(1,1));
			pElem->QueryDoubleAttribute("z", &Wv_(2,2));
		}
		pElem=hRoot.FirstChild("VUKFSettings").FirstChild("Attitude").FirstChild("PreNoiStd").Element();
		if (pElem){
			pElem->QueryDoubleAttribute("x", &Wq_(0,0));
			pElem->QueryDoubleAttribute("y", &Wq_(1,1));
			pElem->QueryDoubleAttribute("z", &Wq_(2,2));
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
	Wv_ = Wv_*Wv_;
	Wq_ = Wq_*Wq_;
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

void FilterVUKF2::logState(){
	  pManager_->ofsLog_ << xs_.t_ << "\t";
	  pManager_->ofsLog_ << xs_.x_.r_(0) << "\t" << xs_.x_.r_(1) << "\t" << xs_.x_.r_(2) << "\t";
	  pManager_->ofsLog_ << xs_.x_.v_(0) << "\t" << xs_.x_.v_(1) << "\t" << xs_.x_.v_(2) << "\t";
	  pManager_->ofsLog_ << xs_.x_.q_(0) << "\t" << xs_.x_.q_(1) << "\t" << xs_.x_.q_(2) << "\t" << xs_.x_.q_(3) << "\t";
	  pManager_->ofsLog_ << xs_.x_.bf_(0) << "\t" << xs_.x_.bf_(1) << "\t" << xs_.x_.bf_(2) << "\t";
	  pManager_->ofsLog_ << xs_.x_.bw_(0) << "\t" << xs_.x_.bw_(1) << "\t" << xs_.x_.bw_(2) << "\t";
	  pManager_->ofsLog_ << xs_.slippageDetection_.flag_[0] << "\t" << xs_.slippageDetection_.flag_[1] << "\t" << xs_.slippageDetection_.flag_[2] << "\t" << xs_.slippageDetection_.flag_[3] << "\t";
	  for(int i=0;i<LSE_VUKF_N;i++){
		  pManager_->ofsLog_ << xs_.P_(i,i) << "\t";
	  }
	  for(int i=0;i<LSE_DOF_LEG*LSE_N_LEG;i++){
		  pManager_->ofsLog_ << xs_.y_(i) << "\t";
	  }
	  pManager_->ofsLog_ << endl;
}

std::string FilterVUKF2::getKeyString(){
	std::ostringstream oss (std::ostringstream::out);
	oss << pManager_->Rda_(0,0) << "_" << kinOutTh_ << "_" << restorationFactor_;
	return oss.str();
}

VUKFFilterState VUKFFilterState::operator +(const Eigen::Matrix<double,VUKFF_state_dim,1> &y) const{
	VUKFFilterState x;
	x.r_ = r_+y.block<3,1>(0,0);
	x.v_ = v_+y.block<3,1>(3,0);
	x.q_ = Rotations::quatL(Rotations::rotVecToQuat(y.block<3,1>(6,0)))*q_;
	x.bf_ = bf_+y.block<3,1>(9,0);
	x.bw_ = bw_+y.block<3,1>(12,0);
	return x;
}

Eigen::Matrix<double,VUKFF_state_dim,1> VUKFFilterState::operator -(const VUKFFilterState &y) const{
	Eigen::Matrix<double,VUKFF_state_dim,1> x;
	x.block<3,1>(0,0) = r_-y.r_;
	x.block<3,1>(3,0) = v_-y.v_;
	x.block<3,1>(6,0) = Rotations::quatToRotVec(Rotations::quatL(q_)*Rotations::quatInverse(y.q_));
	x.block<3,1>(9,0) = bf_-y.bf_;
	x.block<3,1>(12,0) = bw_-y.bw_;
	return x;
}

void FilterVUKF2::predict(VUKFFilterState& x,double dt,ImuMeas imuMeas){
	Eigen::Matrix3d G0T = pManager_->gamma(0,-dt*(imuMeas.w_-x.bw_));
	Eigen::Matrix3d G1T = pManager_->gamma(1,-dt*(imuMeas.w_-x.bw_));
	Eigen::Matrix3d G2T = pManager_->gamma(2,-dt*(imuMeas.w_-x.bw_));
	Eigen::Matrix3d R_IW, R_WI;

	R_WI = Rotations::quatToRotMat(x.q_);
	R_IW = R_WI.transpose();

	x.r_ = G0T*(x.r_+dt*x.v_-dt*dt/2*(2*G2T*(imuMeas.f_-x.bf_)+R_IW*pManager_->g_));
	x.v_ = G0T*(x.v_-dt*(G1T*(imuMeas.f_-x.bf_)+R_IW*pManager_->g_));
	x.q_ = Rotations::quatL(x.q_)*Rotations::rotVecToQuat(dt*(imuMeas.w_-x.bw_));
}


void FilterVUKF2::predict(VUKFFilterState& x,double dt,ImuMeas imuMeas,Eigen::Matrix<double,VUKFF_preNoise_dim,1> n){
	Eigen::Matrix3d G0T = pManager_->gamma(0,-dt*(imuMeas.w_-x.bw_-n.block<3,1>(18,0)));
	Eigen::Matrix3d G1T = pManager_->gamma(1,-dt*(imuMeas.w_-x.bw_-n.block<3,1>(18,0)));
	Eigen::Matrix3d G2T = pManager_->gamma(2,-dt*(imuMeas.w_-x.bw_-n.block<3,1>(18,0)));
	Eigen::Matrix3d R_IW, R_WI;

	R_WI = Rotations::quatToRotMat(x.q_);
	R_IW = R_WI.transpose();

	x.r_ = G0T*(x.r_+dt*x.v_-dt*dt/2*(2*G2T*(imuMeas.f_-x.bf_-n.block<3,1>(15,0))+R_IW*pManager_->g_))+dt*n.block<3,1>(0,0);
	x.v_ = G0T*(x.v_-dt*(G1T*(imuMeas.f_-x.bf_-n.block<3,1>(15,0))+R_IW*pManager_->g_))+dt*n.block<3,1>(3,0);
	x.q_ = Rotations::quatL(x.q_)*Rotations::quatL(Rotations::rotVecToQuat(dt*(imuMeas.w_-x.bw_-n.block<3,1>(18,0))))*Rotations::rotVecToQuat(dt*n.block<3,1>(6,0));
	x.bf_ = x.bf_+dt*n.block<3,1>(9,0);
	x.bw_ = x.bw_+dt*n.block<3,1>(12,0);
}

void FilterVUKF2::samplePredictionNoise(InternState& x,double dt){
	// Prediction noise covariance matrix
	Npre_.setZero();
	Npre_.block<3,3>(0,0) = Wr_/dt;
	Npre_.block<3,3>(3,3) = Wv_/dt;
	Npre_.block<3,3>(6,6) = Wq_/dt;
	Npre_.block<3,3>(9,9) = Wbf_/dt;
	Npre_.block<3,3>(12,12) = Wbw_/dt;
	Npre_.block<3,3>(15,15) = pManager_->Rf_/dt;
	Npre_.block<3,3>(18,18) = pManager_->Rw_/dt;
	Eigen::LLT<MatrixPreCov> lltOfNpre(Npre_);
	SNpre_ = lltOfNpre.matrixL();
	if(lltOfNpre.info()==Eigen::NumericalIssue) std::cout << "Numerical issues while computing Cholesky of Npre_" << std::endl;
	SNpre_ = SNpre_*UKFGamma_;
	x.PN_.block<VUKFF_preNoise_dim,VUKFF_preNoise_dim>(0,0) = SNpre_;
	x.PN_.block<VUKFF_preNoise_dim,VUKFF_preNoise_dim>(0,VUKFF_preNoise_dim) = -SNpre_;
}

void FilterVUKF2::sampleUpdateNoise(InternState& x){
	// Update noise covariance matrix
	Nup_.setZero();
	for(int i=0;i<LSE_N_LEG;i++){
		Nup_.block<3,3>(3*i,3*i) = pManager_->Rda_;
	}
	Eigen::LLT<MatrixUpCov> lltOfNup(Nup_);
	SNup_ = lltOfNup.matrixL();
	if(lltOfNup.info()==Eigen::NumericalIssue) std::cout << "Numerical issues while computing Cholesky of Nup_" << std::endl;
	SNup_ = SNup_*UKFGamma_;
	x.UN_.block<VUKFF_upNoise_dim,VUKFF_upNoise_dim>(0,0) = SNup_;
	x.UN_.block<VUKFF_upNoise_dim,VUKFF_upNoise_dim>(0,VUKFF_upNoise_dim) = -SNup_;
}

}




