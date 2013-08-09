/*!
* @file 	FilterSync.cpp
* @author 	Michael Blösch
* @date		10.10.2012
 */

#include "FilterSync.hpp"
#include "Manager.hpp"
#include "tinyxml.h"
#include <map>
#include <Eigen/Cholesky>

namespace LSE {

FilterSync::FilterSync(Manager* pManager,const char* pFilename){
	pManager_ = pManager;

	// Init all parameters
	xInit_.t_ = 0;
	xInit_.x_.r_.setZero();
	xInit_.x_.v_.setZero();
	xInit_.x_.q_ = Rotations::quatIdentity();
	xInit_.w_.setZero();
	xInit_.x_.p_.setZero();
	xInit_.x_.bf_.setZero();
	xInit_.x_.bw_.setZero();
	for(int i=0;i<LSE_N_LEG;i++){
		xInit_.CFC_[i] = 0;
	}
	xInit_.P_.setIdentity();
	xInit_.f_ = -pManager_->g_;
	Wr_ = 0*Eigen::Matrix3d::Identity();
	Wv_ = 0.003*Eigen::Matrix3d::Identity();
	Wq_ = 0*Eigen::Matrix3d::Identity();
	Wp_ = 0.01*Eigen::Matrix3d::Identity();
	Wbf_ = 0.0001*Eigen::Matrix3d::Identity();
	Wbw_ = 0.000618*Eigen::Matrix3d::Identity();
	Ts_ = 0.01;

	// Flags
	mbEstimateRotBias_ = true;
	mbEstimateAccBias_ = true;

	// UKF Parameters
	UKFAlpha_ = 1e-3;
	UKFKappa_ = 0;
	UKFBeta_ = 2;

	// Compute derived UKF parameters (weights and lambda)
	int L = SF_state_dim+SF_preNoise_dim+SF_upNoise_dim;
	UKFLambda_ = UKFAlpha_*UKFAlpha_*(L+UKFKappa_)-L;
	UKFGamma_ = sqrt(UKFLambda_ + L);
	UKFWs_ = UKFLambda_/(L+UKFLambda_);
	UKFWc_ = UKFLambda_/(L+UKFLambda_)+(1-UKFAlpha_*UKFAlpha_+UKFBeta_);
	UKFWi_ = 1/(2*(L+UKFLambda_));

	loadParam(pFilename);

	compDiscretizedNoiseMat();

}

FilterSync::~FilterSync(){
}

void FilterSync::update(const double& t){
	// Nothing to handle
	update();
}

void FilterSync::update(){
	// Get actual measurements
	ImuMeas imuMeas;
	imuMeas.f_ = pManager_->g_;
	imuMeas.w_.setZero();
	EncMeas encMeas;
	encMeas.CF_[0] = 0;
	encMeas.CF_[1] = 0;
	encMeas.CF_[2] = 0;
	encMeas.CF_[3] = 0;
	encMeas.e_.setZero();
	if(!pManager_->imuMeasList_.empty()){
		std::map<double,ImuMeas>::iterator itImu;
		itImu = pManager_->imuMeasList_.end();
		itImu--;
		imuMeas = itImu->second;
	}
	if(!pManager_->encMeasList_.empty()){
		std::map<double,EncMeas>::iterator itEnc;
		itEnc = pManager_->encMeasList_.end();
		itEnc--;
		encMeas = itEnc->second;
	}

	// Handle case where Bias estimation disabled
	if(!mbEstimateAccBias_) x_.x_.bf_.setZero();
	if(!mbEstimateRotBias_) x_.x_.bw_.setZero();

	// Compute forward kinematics
	Eigen::Matrix<double,SF_upNoise_dim,1> s;
	for(int i=0;i<LSE_N_LEG;i++){
		// I_r_IF = C(q_IB)*(-B_r_BI + B_r_BK + C'(q_KB)*K_r_KF
		s.block(3*i,0,3,1) = Rotations::quatToRotMat(pManager_->q_IB_)*(-pManager_->B_r_BI_+pManager_->B_r_BK_+Rotations::quatToRotMat(pManager_->q_KB_).transpose()*(*pManager_->legKin)(encMeas.e_.col(i),i));
	}

	// Update Contact count
	for(int i=0;i<LSE_N_LEG;i++){
		if(encMeas.CF_[i]){
			x_.CFC_[i]++;
		} else {
			x_.CFC_[i] = 0;
		}
	}

	// Handle Initialization of new contacts (first part)
	for(int i=0;i<LSE_N_LEG;i++){
		if(x_.CFC_[i] <= 1){
			x_.x_.p_.col(i) = x_.x_.r_+Rotations::quatToRotMat(x_.x_.q_).transpose()*s.block<3,1>(3*i,0);
			x_.P_.block(0,15+3*i,SF_state_dim,3).setZero();
			x_.P_.block(15+3*i,0,3,SF_state_dim).setZero();
			x_.P_.block(15+3*i,15+3*i,3,3) = Eigen::Matrix3d::Identity()*1e6;
		}
	}

	// Sample Sigma Points
	// Cholesky decomposition of covariance matrix
	Eigen::LLT<MatrixP> lltOfP(x_.P_);
	if(lltOfP.info()==Eigen::NumericalIssue) std::cout << "Numerical issues while computing Cholesky of P" << std::endl;
	SP_ = lltOfP.matrixL();
	SP_ = SP_*UKFGamma_;
	X_[0] = x_.x_;
	for(int i=1;i<=SF_state_dim;i++){
		X_[i] = x_.x_+SP_.col(i-1);
		X_[i+SF_state_dim] = x_.x_+(-1*SP_.col(i-1));
	}

	// Propagate Sigma Points
	// State Part
	for(int i=0;i<=2*SF_state_dim;i++){
		predict(X_[i],Ts_,imuMeas);
	}
	// Prediction noise Part
	for(int i=1;i<=SF_preNoise_dim;i++){
		Eigen::Matrix<double,SF_preNoise_dim,1> n = SNpre_.col(i-1);

		// Handle case where Bias estimation disabled
		if(!mbEstimateAccBias_) n.block<3,1>(9,0).setZero();
		if(!mbEstimateRotBias_) n.block<3,1>(12,0).setZero();

		X_[2*SF_state_dim+i] = x_.x_;
		predict(X_[2*SF_state_dim+i],Ts_,imuMeas,n);
		n = -n;
		X_[2*SF_state_dim+SF_preNoise_dim+i] = x_.x_;
		predict(X_[2*SF_state_dim+SF_preNoise_dim+i],Ts_,imuMeas,n);
	}

	// Compute predicted state and covariance
	Eigen::Matrix<double,SF_state_dim,1> vec;
	vec.setZero();
	for(int i=1;i<=2*(SF_state_dim+SF_preNoise_dim);i++){
		vec = vec + (UKFWi_*(X_[i]-X_[0]));
	}
	x_.x_ = X_[0]+vec;
	vec = X_[0]-x_.x_;
	x_.P_ = (UKFWc_+2*SF_upNoise_dim*UKFWi_)*vec*vec.transpose();
	for(int i=1;i<=2*(SF_state_dim+SF_preNoise_dim);i++){
		vec = X_[i]-x_.x_;
		x_.P_ += UKFWi_*vec*vec.transpose();
	}
	x_.t_ += Ts_;

	// Project through measurement function
	SyncFilterState filterState;
	Eigen::Matrix<double,SF_upNoise_dim,1> upNoi;
	Eigen::Matrix<double,SF_upNoise_dim,1+2*(SF_state_dim+SF_preNoise_dim+SF_upNoise_dim)> Y;
	for(int i=0;i<=2*(SF_state_dim+SF_preNoise_dim+SF_upNoise_dim);i++){
		if(i<=2*(SF_state_dim+SF_preNoise_dim)){
			filterState = X_[i];
			upNoi.setZero();
		} else if (i<=2*(SF_state_dim+SF_preNoise_dim)+SF_upNoise_dim) {
			filterState = X_[0];
			upNoi = SNup_.col((i-2*(SF_state_dim+SF_preNoise_dim)-1));
		} else {
			filterState = X_[0];
			upNoi = -SNup_.col((i-2*(SF_state_dim+SF_preNoise_dim)-1-SF_upNoise_dim));
		}
		for(int j=0;j<LSE_N_LEG;j++){
			Y.block<3,1>(j*3,i) = s.block<3,1>(j*3,0)-Rotations::quatToRotMat(filterState.q_)*(filterState.p_.col(j)-filterState.r_)+upNoi.block<3,1>(j*3,0);
		}
	}

	// Compute innovation and corresponding covariance
	Eigen::Matrix<double,SF_upNoise_dim,1> y;
	y = UKFWs_*Y.col(0);
	for(int i=1;i<=2*(SF_state_dim+SF_preNoise_dim+SF_upNoise_dim);i++){
		y += UKFWi_*Y.col(i);
	}
	Eigen::Matrix<double,SF_upNoise_dim,SF_upNoise_dim> Py;
	Py = UKFWc_*(Y.col(0)-y)*(Y.col(0)-y).transpose();
	for(int i=1;i<=2*(SF_state_dim+SF_preNoise_dim+SF_upNoise_dim);i++){
		Py += UKFWi_*(Y.col(i)-y)*(Y.col(i)-y).transpose();
	}

	// Compute cross-correlation
	vec = X_[0]-x_.x_;
	Eigen::Matrix<double,SF_state_dim,SF_upNoise_dim> Pxy;
	Pxy = UKFWc_*vec*(Y.col(0)-y).transpose();
	for(int i=1;i<=2*(SF_state_dim+SF_preNoise_dim+SF_upNoise_dim);i++){
		if(i<=2*(SF_state_dim+SF_preNoise_dim)){
			filterState = X_[i];
		} else {
			filterState = X_[0];
		}
		vec = filterState-x_.x_;
		Pxy += UKFWi_*vec*(Y.col(i)-y).transpose();
	}

	// Compute inverse of innovation covariance
	Eigen::Matrix<double,SF_upNoise_dim,SF_upNoise_dim> Pyinv = Py.inverse();
	for(int i=0;i<LSE_N_LEG;i++){
		if(x_.CFC_[i]==0){
			Pyinv.block<3,SF_upNoise_dim>(3*i,0).setZero();
			Pyinv.block<SF_upNoise_dim,3>(0,3*i).setZero();
		}
	}

	// Compute Kalman Gain
	Eigen::Matrix<double,SF_state_dim,SF_upNoise_dim> K = Pxy*Pyinv;

	// Update state and covariance matrix
	vec = -K*y;
	x_.x_ = x_.x_+vec;
	x_.P_ = x_.P_ - K*Py*K.transpose();

	// Update accelerometer and gyroscope
	x_.f_ = imuMeas.f_-x_.x_.bf_;
	x_.w_ = imuMeas.w_-x_.x_.bw_;

	// Avoid singular P
	if(!mbEstimateAccBias_) x_.P_.block<3,3>(9,9) = xInit_.P_.block<3,3>(9,9);
	if(!mbEstimateRotBias_) x_.P_.block<3,3>(12,12) = xInit_.P_.block<3,3>(12,12);
}

State FilterSync::getEst(){
	Eigen::Matrix3d R_WI,R_IB;
	R_WI = Rotations::quatToRotMat(x_.x_.q_).transpose();
	R_IB = Rotations::quatToRotMat(pManager_->q_IB_);
	State x = State();
	x.t_ = x_.t_;
	x.r_ = x_.x_.r_-R_WI*R_IB*pManager_->B_r_BI_;
	x.v_ = x_.x_.v_-R_WI*(Rotations::vecToSqew(x_.w_)*R_IB*pManager_->B_r_BI_);
	x.q_ = Rotations::quatL(pManager_->q_IB_).transpose()*x_.x_.q_;
	x.w_ = R_IB.transpose()*x_.w_;
	x.P_.setZero();
	x.P_.block(0,0,9,9) = x_.P_.block(0,0,9,9);
	x.P_.block(9,9,3,3) = x_.P_.block(12,12,3,3)+pManager_->Rw_;
	return x;
}

void FilterSync::resetEstimate(const double& t){
	x_ = xInit_;
	x_.t_ = t;
}

void FilterSync::loadParam(const char* pFilename){
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
		pElem=hRoot.FirstChild("SyncSettings").FirstChild("Position").FirstChild("Init").Element();
		if (pElem){
			pElem->QueryDoubleAttribute("x", &xInit_.x_.r_(0));
			pElem->QueryDoubleAttribute("y", &xInit_.x_.r_(1));
			pElem->QueryDoubleAttribute("z", &xInit_.x_.r_(2));
		}
		pElem=hRoot.FirstChild("SyncSettings").FirstChild("Position").FirstChild("InitStd").Element();
		if (pElem){
			pElem->QueryDoubleAttribute("x", &xInit_.P_(0,0));
			pElem->QueryDoubleAttribute("y", &xInit_.P_(1,1));
			pElem->QueryDoubleAttribute("z", &xInit_.P_(2,2));
		}

		pElem=hRoot.FirstChild("SyncSettings").FirstChild("Velocity").FirstChild("Init").Element();
		if (pElem){
			pElem->QueryDoubleAttribute("x", &xInit_.x_.v_(0));
			pElem->QueryDoubleAttribute("y", &xInit_.x_.v_(1));
			pElem->QueryDoubleAttribute("z", &xInit_.x_.v_(2));
		}
		pElem=hRoot.FirstChild("SyncSettings").FirstChild("Velocity").FirstChild("InitStd").Element();
		if (pElem){
			pElem->QueryDoubleAttribute("x", &xInit_.P_(3,3));
			pElem->QueryDoubleAttribute("y", &xInit_.P_(4,4));
			pElem->QueryDoubleAttribute("z", &xInit_.P_(5,5));
		}

		pElem=hRoot.FirstChild("SyncSettings").FirstChild("Attitude").FirstChild("Init").Element();
		if (pElem){
			pElem->QueryDoubleAttribute("x", &xInit_.x_.q_(0));
			pElem->QueryDoubleAttribute("y", &xInit_.x_.q_(1));
			pElem->QueryDoubleAttribute("z", &xInit_.x_.q_(2));
			pElem->QueryDoubleAttribute("w", &xInit_.x_.q_(3));
		}
		pElem=hRoot.FirstChild("SyncSettings").FirstChild("Attitude").FirstChild("InitStd").Element();
		if (pElem){
			pElem->QueryDoubleAttribute("x", &xInit_.P_(6,6));
			pElem->QueryDoubleAttribute("y", &xInit_.P_(7,7));
			pElem->QueryDoubleAttribute("z", &xInit_.P_(8,8));
		}

		pElem=hRoot.FirstChild("SyncSettings").FirstChild("AccelerometerBias").FirstChild("Init").Element();
		if (pElem){
			pElem->QueryDoubleAttribute("x", &xInit_.x_.bf_(0));
			pElem->QueryDoubleAttribute("y", &xInit_.x_.bf_(1));
			pElem->QueryDoubleAttribute("z", &xInit_.x_.bf_(2));
		}
		pElem=hRoot.FirstChild("SyncSettings").FirstChild("AccelerometerBias").FirstChild("InitStd").Element();
		if (pElem){
			pElem->QueryDoubleAttribute("x", &xInit_.P_(9,9));
			pElem->QueryDoubleAttribute("y", &xInit_.P_(10,10));
			pElem->QueryDoubleAttribute("z", &xInit_.P_(11,11));
		}

		pElem=hRoot.FirstChild("SyncSettings").FirstChild("GyroscopeBias").FirstChild("Init").Element();
		if (pElem){
			pElem->QueryDoubleAttribute("x", &xInit_.x_.bw_(0));
			pElem->QueryDoubleAttribute("y", &xInit_.x_.bw_(1));
			pElem->QueryDoubleAttribute("z", &xInit_.x_.bw_(2));
		}
		pElem=hRoot.FirstChild("SyncSettings").FirstChild("GyroscopeBias").FirstChild("InitStd").Element();
		if (pElem){
			pElem->QueryDoubleAttribute("x", &xInit_.P_(12,12));
			pElem->QueryDoubleAttribute("y", &xInit_.P_(13,13));
			pElem->QueryDoubleAttribute("z", &xInit_.P_(14,14));
		}
		pElem=hRoot.FirstChild("SyncSettings").FirstChild("Position").FirstChild("PreNoiStd").Element();
		if (pElem){
			pElem->QueryDoubleAttribute("x", &Wr_(0,0));
			pElem->QueryDoubleAttribute("y", &Wr_(1,1));
			pElem->QueryDoubleAttribute("z", &Wr_(2,2));
		}
		pElem=hRoot.FirstChild("SyncSettings").FirstChild("Velocity").FirstChild("PreNoiStd").Element();
		if (pElem){
			pElem->QueryDoubleAttribute("x", &Wv_(0,0));
			pElem->QueryDoubleAttribute("y", &Wv_(1,1));
			pElem->QueryDoubleAttribute("z", &Wv_(2,2));
		}
		pElem=hRoot.FirstChild("SyncSettings").FirstChild("Attitude").FirstChild("PreNoiStd").Element();
		if (pElem){
			pElem->QueryDoubleAttribute("x", &Wq_(0,0));
			pElem->QueryDoubleAttribute("y", &Wq_(1,1));
			pElem->QueryDoubleAttribute("z", &Wq_(2,2));
		}
		pElem=hRoot.FirstChild("SyncSettings").FirstChild("Foothold").FirstChild("PreNoiStd").Element();
		if (pElem){
			pElem->QueryDoubleAttribute("x", &Wp_(0,0));
			pElem->QueryDoubleAttribute("y", &Wp_(1,1));
			pElem->QueryDoubleAttribute("z", &Wp_(2,2));
		}
		pElem=hRoot.FirstChild("SyncSettings").FirstChild("AccelerometerBias").FirstChild("PreNoiStd").Element();
		if (pElem){
			pElem->QueryDoubleAttribute("x", &Wbf_(0,0));
			pElem->QueryDoubleAttribute("y", &Wbf_(1,1));
			pElem->QueryDoubleAttribute("z", &Wbf_(2,2));
		}
		pElem=hRoot.FirstChild("SyncSettings").FirstChild("GyroscopeBias").FirstChild("PreNoiStd").Element();
		if (pElem){
			pElem->QueryDoubleAttribute("x", &Wbw_(0,0));
			pElem->QueryDoubleAttribute("y", &Wbw_(1,1));
			pElem->QueryDoubleAttribute("z", &Wbw_(2,2));
		}
		pElem=hRoot.FirstChild("SyncSettings").Element();
		if (pElem){
			pElem->QueryDoubleAttribute("timeStepping", &Ts_);
		}
		pElem=hRoot.FirstChild("SyncSettings").FirstChild("AccelerometerBias").Element();
		if (pElem){
			pElem->QueryIntAttribute("estimate", &mInt);
			mbEstimateAccBias_ = mInt;
		}
		pElem=hRoot.FirstChild("SyncSettings").FirstChild("GyroscopeBias").Element();
		if (pElem){
			pElem->QueryIntAttribute("estimate", &mInt);
			mbEstimateRotBias_ = mInt;
		}
	}

	xInit_.P_ = xInit_.P_*xInit_.P_;
	Wr_ = Wr_*Wr_;
	Wv_ = Wv_*Wv_;
	Wq_ = Wq_*Wq_;
	Wp_ = Wp_*Wp_;
	Wbf_ = Wbf_*Wbf_;
	Wbw_ = Wbw_*Wbw_;
}

SyncFilterState SyncFilterState::operator +(const Eigen::Matrix<double,SF_state_dim,1> &y) const{
	SyncFilterState x;
	x.r_ = r_+y.block<3,1>(0,0);
	x.v_ = v_+y.block<3,1>(3,0);
	x.q_ = Rotations::quatL(Rotations::rotVecToQuat(y.block<3,1>(6,0)))*q_;
	x.bf_ = bf_+y.block<3,1>(9,0);
	x.bw_ = bw_+y.block<3,1>(12,0);
	for(int i=0;i<LSE_N_LEG;i++){
		x.p_.col(i) = p_.col(i)+y.block<3,1>(15+3*i,0);
	}
	return x;
}

Eigen::Matrix<double,SF_state_dim,1> SyncFilterState::operator -(const SyncFilterState &y) const{
	Eigen::Matrix<double,SF_state_dim,1> x;
	x.block<3,1>(0,0) = r_-y.r_;
	x.block<3,1>(3,0) = v_-y.v_;
	x.block<3,1>(6,0) = Rotations::quatToRotVec(Rotations::quatL(q_)*Rotations::quatInverse(y.q_));
	x.block<3,1>(9,0) = bf_-y.bf_;
	x.block<3,1>(12,0) = bw_-y.bw_;
	for(int i=0;i<LSE_N_LEG;i++){
		x.block<3,1>(15+3*i,0) = p_.col(i)-y.p_.col(i);
	}
	return x;
}

void FilterSync::predict(SyncFilterState& x,double Ts,ImuMeas imuMeas){
	x.r_ = x.r_+Ts*x.v_;
	x.v_ = x.v_+Ts*(Rotations::quatToRotMat(x.q_).transpose()*(imuMeas.f_-x.bf_)+pManager_->g_);
	x.q_ = Rotations::quatL(Rotations::rotVecToQuat(-Ts*(imuMeas.w_-x.bw_)))*x.q_;
}


void FilterSync::predict(SyncFilterState& x,double Ts,ImuMeas imuMeas,Eigen::Matrix<double,SF_preNoise_dim,1> n){
	x.r_ = x.r_+Ts*x.v_+Ts*n.block<3,1>(0,0);
	x.v_ = x.v_+Ts*(Rotations::quatToRotMat(x.q_).transpose()*(imuMeas.f_-x.bf_-n.block<3,1>(15,0))+pManager_->g_)+Ts*n.block<3,1>(3,0);
	x.q_ = Rotations::quatL(Rotations::rotVecToQuat(Ts*n.block<3,1>(6,0)))*Rotations::quatL(Rotations::rotVecToQuat(-Ts*(imuMeas.w_-x.bw_-n.block<3,1>(18,0))))*x.q_;
	x.bf_ = x.bf_+Ts*n.block<3,1>(9,0);
	x.bw_ = x.bw_+Ts*n.block<3,1>(12,0);
	for(int j=0;j<LSE_N_LEG;j++){
		x.p_.col(j) = x.p_.col(j)+Ts*n.block<3,1>(21+3*j,0);
	}
}

void FilterSync::setSamplingTime(double Ts){
	Ts_ = Ts;
	compDiscretizedNoiseMat();
}

void FilterSync::compDiscretizedNoiseMat(){
	// Prediction noise covariance matrix
	Npre_.setZero();
	Npre_.block<3,3>(0,0) = Wr_/Ts_;
	Npre_.block<3,3>(3,3) = Wv_/Ts_;
	Npre_.block<3,3>(6,6) = Wq_/Ts_;
	Npre_.block<3,3>(9,9) = Wbf_/Ts_;
	Npre_.block<3,3>(12,12) = Wbw_/Ts_;
	Npre_.block<3,3>(15,15) = pManager_->Rf_/Ts_;
	Npre_.block<3,3>(18,18) = pManager_->Rw_/Ts_;
	for(int i=0;i<LSE_N_LEG;i++){
		Npre_.block<3,3>(21+3*i,21+3*i) = Wp_/Ts_;
	}
	Eigen::LLT<MatrixPreCov> lltOfNpre(Npre_);
	SNpre_ = lltOfNpre.matrixL();
	if(lltOfNpre.info()==Eigen::NumericalIssue) std::cout << "Numerical issues while computing Cholesky of Npre_" << std::endl;
	SNpre_ = SNpre_*UKFGamma_;

	// Update noise covariance matrix
	Nup_.setZero();
	for(int i=0;i<LSE_N_LEG;i++){
		Nup_.block<3,3>(3*i,3*i) = pManager_->Rs_;
	}

	Eigen::LLT<MatrixUpCov> lltOfNup(Nup_);
	SNup_ = lltOfNup.matrixL();
	if(lltOfNup.info()==Eigen::NumericalIssue) std::cout << "Numerical issues while computing Cholesky of Nup_" << std::endl;
	SNup_ = SNup_*UKFGamma_;
}

}




