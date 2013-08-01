/*!
* @file 	DelayCalibration.cpp
* @author 	Michael Blösch
* @date		10.10.2012
 */

#include "DelayCalibration.hpp"
#include "Manager.hpp"
#include "tinyxml.h"

namespace LSE {

DelayCalibration::DelayCalibration(Manager* pManager,const char* pFilename){
	pManager_ = pManager;

	// Timesteps and stuff
	t1_ = 0;
	t2_ = 0;
	N_ = 0;
	TsImu_ = 0.01;

	// Norm vector and iterators
	imuRateNorm_.clear();
	encRateNorm_.clear();
	posRateNorm_.clear();

	// Parameters
	dt_ = 0.001;
	maxDelay_ = 0.5;
	difWindowEnc_ = 2;
	difWindowPos_ = 2;
	mbUseImu_ = 0;
	mbUseEnc_ = 0;
	mbUsePos_ = 0;

	loadParam(pFilename);
}

DelayCalibration::~DelayCalibration(){
}

int DelayCalibration::calibrateDelay(const double& t,const double& T){
	// Check conditions
	if(T < maxDelay_ || dt_ > T) return 0;
	if((int)mbUseImu_+(int)mbUseEnc_+(int)mbUsePos_<2) return 0;

	// Initialize and recheck
	if(initialize(t,T)==0) return 0;
	if(t2_-t1_ < maxDelay_ || dt_ > t2_-t1_) return 0;

	if(mbUseImu_) getImuNorms();
	if(mbUseEnc_) getEncNorms();
	if(mbUsePos_) getPosNorms();

	int M = std::max(floor(N_/2),N_-2*ceil(maxDelay_/dt_));
	double cc;
	double ccMax;
	int ccMaxIndex;

	// Align Encoders and IMU by evaluating crosscorrelation
	if(mbUseImu_ && mbUseEnc_){
		ccMax = 0;
		ccMaxIndex = 0;
		for(int i=0;i<N_-M+1;i++){
			cc = 0;
			for(int j=0;j<M;j++){
				cc += encRateNorm_[floor((N_-M)/2)+j]*imuRateNorm_[i+j];
			}
			if(cc>ccMax){
				ccMax = cc;
				ccMaxIndex = i-floor((N_-M)/2);
			}
		}
		pManager_->tImu_ = 0;
		pManager_->tEnc_ = ccMaxIndex*dt_;
	}

	// Align 6DOF pose and IMU by evaluating crosscorrelation
	if(mbUseImu_ && mbUsePos_){
		ccMax = 0;
		ccMaxIndex = 0;
		for(int i=0;i<N_-M+1;i++){
			cc = 0;
			for(int j=0;j<M;j++){
				cc += posRateNorm_[floor((N_-M)/2)+j]*imuRateNorm_[i+j];
			}
			if(cc>ccMax){
				ccMax = cc;
				ccMaxIndex = i-floor((N_-M)/2);
			}
		}
		pManager_->tImu_ = 0;
		pManager_->tPos_ = ccMaxIndex*dt_;
	} else if(mbUseEnc_ && mbUsePos_){
		ccMax = 0;
		ccMaxIndex = 0;
		for(int i=0;i<N_-M+1;i++){
			cc = 0;
			for(int j=0;j<M;j++){
				cc += posRateNorm_[floor((N_-M)/2)+j]*encRateNorm_[i+j];
			}
			if(cc>ccMax){
				ccMax = cc;
				ccMaxIndex = i-floor((N_-M)/2);
			}
		}
		pManager_->tEnc_ = 0;
		pManager_->tPos_ = ccMaxIndex*dt_;
	}

	return 1;
}

int DelayCalibration::initialize(const double& t,const double& T){
	t1_ = t-T;
	t2_ = t;
	N_ = 0;

	// Abort if no measurements
	if(mbUseImu_ && pManager_->imuMeasList_.empty())
		return 0;
	if(mbUseEnc_ && pManager_->encMeasList_.empty())
		return 0;
	if(mbUsePos_ && pManager_->posMeasList_.empty())
		return 0;


	// Estimate frequency of IMU measurements
	if(mbUseImu_){
		int countImu = -1;
		for(itImu_ = pManager_->imuMeasList_.lower_bound(t1_);itImu_ != pManager_->imuMeasList_.lower_bound(t2_);itImu_++){
			countImu++;
		}
		if(countImu>0){
			itImu_ = pManager_->imuMeasList_.lower_bound(t2_);
			itImu_--;
			TsImu_ = (itImu_->first-pManager_->imuMeasList_.lower_bound(t1_)->first)/countImu;
		}
	}

	// Maximal range -> should avoid any subsequent problem
	if(mbUseImu_){
		t1_ = std::max(t1_,(pManager_->imuMeasList_.begin())->first-0.5*TsImu_);
		t2_ = std::min(t2_,(pManager_->imuMeasList_.rbegin())->first-0.5*TsImu_);
	}
	if(mbUseEnc_){
		// min
		itEnc_ = pManager_->encMeasList_.begin();
		itEnc2_ = pManager_->encMeasList_.begin();
		for(int i=0;i<difWindowEnc_;i++){
			if(itEnc2_ != pManager_->encMeasList_.end()) itEnc2_++;
			while((*itEnc2_).second.CF_.sum() < 3 && itEnc2_ != pManager_->encMeasList_.end()) itEnc2_++;
		}
		if(itEnc2_ == pManager_->encMeasList_.end()) return 0;
		t1_ = std::max(t1_,0.5*(itEnc_->first+itEnc2_->first));

		//max
		ritEnc_ = pManager_->encMeasList_.rbegin();
		ritEnc2_ = pManager_->encMeasList_.rbegin();
		for(int i=0;i<difWindowEnc_;i++){
			if(ritEnc_ != pManager_->encMeasList_.rend()) ritEnc_++;
			while((*ritEnc_).second.CF_.sum() < 3 && ritEnc_ != pManager_->encMeasList_.rend()) ritEnc_++;
		}
		if(ritEnc_ == pManager_->encMeasList_.rend()) return 0;
		t2_ = std::min(t2_,0.5*(ritEnc_->first+ritEnc2_->first));

		// TODO: should also check for available full stance
	}
	if(mbUsePos_){
		// min
		itPos_ = pManager_->posMeasList_.begin();
		itPos2_ = pManager_->posMeasList_.begin();
		for(int i=0;i<difWindowPos_;i++){
			if(itPos2_ != pManager_->posMeasList_.end()) itPos2_++;
		}
		if(itPos2_ == pManager_->posMeasList_.end()) return 0;
		t1_ = std::max(t1_,0.5*(itPos_->first+itPos2_->first));

		//max
		ritPos_ = pManager_->posMeasList_.rbegin();
		ritPos2_ = pManager_->posMeasList_.rbegin();
		for(int i=0;i<difWindowPos_;i++){
			if(ritPos_ != pManager_->posMeasList_.rend()) ritPos_++;
		}
		if(ritPos_ == pManager_->posMeasList_.rend()) return 0;
		t2_ = std::min(t2_,0.5*(ritPos_->first+ritPos2_->first));
	}

	// Check range
	N_ = floor((t2_-t1_)/dt_)+1;
	t2_ = t1_+(N_-1)*dt_;
	if(t1_>=t2_) return 0;

	// Allocate space for rotational rate norm vectors
	if(mbUseImu_){
		imuRateNorm_.clear();
		imuRateNorm_.reserve(N_);
	}
	if(mbUseEnc_){
		encRateNorm_.clear();
		encRateNorm_.reserve(N_);
	}
	if(mbUsePos_){
		posRateNorm_.clear();
		posRateNorm_.reserve(N_);
	}
	return 1;
}

void DelayCalibration::getImuNorms(){
	// Interpolate IMU norms
	itImu_ = pManager_->imuMeasList_.upper_bound(t1_+0.5*TsImu_);
	itImu_--;
	lastTime_ = (*itImu_).first-0.5*TsImu_;
	lastNorm_ = (*itImu_).second.w_.norm();
	itImu_++;
	interpolTime_ = t1_;
	while(lastTime_ < t2_ && itImu_ != pManager_->imuMeasList_.end()){
		newTime_ = (*itImu_).first-0.5*TsImu_;
		newNorm_ = (*itImu_).second.w_.norm();
		while(newTime_ >= interpolTime_){
			imuRateNorm_.push_back(lastNorm_+(newNorm_-lastNorm_)/(newTime_-lastTime_)*(interpolTime_-lastTime_));
			interpolTime_ += dt_;
		}
		lastTime_ = newTime_;
		lastNorm_ = newNorm_;
		itImu_++;
	}
}

void DelayCalibration::getEncNorms(){
	// Get footholds based on first full stance
	Eigen::Matrix<double,3,LSE_N_LEG> p;
	for(int i=0;i<LSE_N_LEG;i++){
		itEnc_ = pManager_->encMeasList_.lower_bound(t1_);
		while((*itEnc_).second.CF_.sum() < LSE_N_LEG && itEnc_ != pManager_->encMeasList_.end()) itEnc_++;
		p.col(i) = (*pManager_->legKin)((*itEnc_).second.e_.col(i),i);
	}

	// Calculate and interpolate enc norms
	// Initialize iterators
	itEnc_ = pManager_->encMeasList_.lower_bound(t1_);
	if(itEnc_ != pManager_->encMeasList_.begin()) itEnc_--;
	while((*itEnc_).second.CF_.sum() < 3 && itEnc_ != pManager_->encMeasList_.begin()) itEnc_--;
	itEnc2_ = itEnc_;
	for(int i=0;i<difWindowEnc_;i++){
		if(itEnc2_ != pManager_->encMeasList_.end()) itEnc2_++;
		while((*itEnc2_).second.CF_.sum() < 3 && itEnc2_ != pManager_->encMeasList_.end()) itEnc2_++;
	}
	while(0.5*(itEnc_->first+itEnc2_->first)>t1_ && itEnc_ != pManager_->encMeasList_.begin() && itEnc2_ != pManager_->encMeasList_.begin()){
		itEnc_--;
		itEnc2_--;
		while((*itEnc_).second.CF_.sum() < 3 && itEnc_ != pManager_->encMeasList_.begin()) itEnc_--;
		while((*itEnc2_).second.CF_.sum() < 3 && itEnc2_ != pManager_->encMeasList_.begin()) itEnc2_--;
	}
	// Compute rotational rate and corresponding time
	lastTime_ = ((*itEnc_).first+(*itEnc2_).first)/2;
	q1_ = quatFromFootholds((*itEnc_).second,p);
	q2_ = quatFromFootholds((*itEnc2_).second,p);
	lastNorm_ = Rotations::rangePi(Rotations::quatToRotVec(Rotations::quatR(q1_).transpose()*q2_)).norm()/((*itEnc2_).first-(*itEnc_).first);
	// Start computation loop
	if(itEnc_ != pManager_->encMeasList_.end()) itEnc_++;
	if(itEnc2_ != pManager_->encMeasList_.end()) itEnc2_++;
	interpolTime_ = t1_;
	while(lastTime_ < t2_ && itEnc2_ != pManager_->encMeasList_.end()){
		// Jump over bad stances
		while((*itEnc_).second.CF_.sum() < 3 && itEnc_ != pManager_->encMeasList_.end()) itEnc_++;
		while((*itEnc2_).second.CF_.sum() < 3 && itEnc2_ != pManager_->encMeasList_.end()) itEnc2_++;
		// Compute rotational rate and corresponding time
		newTime_ = ((*itEnc_).first+(*itEnc2_).first)/2;
		q1_ = quatFromFootholds((*itEnc_).second,p);
		q2_ = quatFromFootholds((*itEnc2_).second,p);
		newNorm_ = (Rotations::rangePi(Rotations::quatToRotVec(Rotations::quatR(q1_).transpose()*q2_)).norm()/((*itEnc2_).first-(*itEnc_).first));
		while(newTime_ >= interpolTime_){
			encRateNorm_.push_back(lastNorm_+(newNorm_-lastNorm_)/(newTime_-lastTime_)*(interpolTime_-lastTime_));
			interpolTime_ += dt_;
		}
		lastTime_ = newTime_;
		lastNorm_ = newNorm_;
		if(itEnc_ != pManager_->encMeasList_.end()) itEnc_++;
		if(itEnc2_ != pManager_->encMeasList_.end()) itEnc2_++;
	}
}

void DelayCalibration::getPosNorms(){
	// Calculate and interpolate pos norms
	// Initialize iterators
	itPos_ = pManager_->posMeasList_.lower_bound(t1_);
	if(itPos_ != pManager_->posMeasList_.begin()) itPos_--;
	itPos2_ = itPos_;
	for(int i=0;i<difWindowPos_;i++){
		if(itPos2_ != pManager_->posMeasList_.end()) itPos2_++;
	}
	while(0.5*(itPos_->first+itPos2_->first)>t1_ && itPos_ != pManager_->posMeasList_.begin() && itPos2_ != pManager_->posMeasList_.begin()){
		itPos_--;
		itPos2_--;
	}
	// Compute rotational rate and corresponding time
	lastTime_ = ((*itPos_).first+(*itPos2_).first)/2;
	q1_ = (*itPos_).second.q_;
	q2_ = (*itPos2_).second.q_;
	lastNorm_ = Rotations::rangePi(Rotations::quatToRotVec(Rotations::quatR(q1_).transpose()*q2_)).norm()/((*itPos2_).first-(*itPos_).first);
	// Start computation loop
	if(itPos_ != pManager_->posMeasList_.end()) itPos_++;
	if(itPos2_ != pManager_->posMeasList_.end()) itPos2_++;
	interpolTime_ = t1_;
	while(lastTime_ < t2_ && itPos2_ != pManager_->posMeasList_.end()){
		// Compute rotational rate and corresponding time
		newTime_ = ((*itPos_).first+(*itPos2_).first)/2;
		q1_ = (*itPos_).second.q_;
		q2_ = (*itPos2_).second.q_;
		newNorm_ = (Rotations::rangePi(Rotations::quatToRotVec(Rotations::quatR(q1_).transpose()*q2_)).norm()/((*itPos2_).first-(*itPos_).first));
		while(newTime_ >= interpolTime_){
			posRateNorm_.push_back(lastNorm_+(newNorm_-lastNorm_)/(newTime_-lastTime_)*(interpolTime_-lastTime_));
			interpolTime_ += dt_;
		}
		lastTime_ = newTime_;
		lastNorm_ = newNorm_;
		if(itPos_ != pManager_->posMeasList_.end()) itPos_++;
		if(itPos2_ != pManager_->posMeasList_.end()) itPos2_++;
	}
}

Rotations::Quat DelayCalibration::quatFromFootholds(const EncMeas& m,const Eigen::Matrix<double,3,LSE_N_LEG>& p){
	Rotations::Quat q;
	Rotations::Quat a;
	Rotations::Quat b;
	Eigen::Matrix4d A;
	Eigen::Matrix4d B;
	A.setZero();
	for(int i=0;i<LSE_N_LEG;i++){
		if(m.CF_[i]){
			for(int j=i+1;j<LSE_N_LEG;j++){
				if(m.CF_[j]){
					a.block(0,0,3,1) = p.col(i)-p.col(j);
					a(3) = 0;
					b.block(0,0,3,1) = (*pManager_->legKin)(m.e_.col(i),i)-(*pManager_->legKin)(m.e_.col(j),j);
					b(3) = 0;
					B = Rotations::quatR(a)-Rotations::quatL(b);
					A = A + B.transpose()*B;
				}
			}
		}
	}
	Eigen::SelfAdjointEigenSolver<Eigen::Matrix4d> es;
	es.compute(A);
	Eigen::Vector4d D = es.eigenvalues();
	Eigen::Matrix4d V = es.eigenvectors();

	// Find minimal eigenvalue and corresponding eigenvector
	int minInd;
	D.minCoeff(&minInd);
	q = V.col(minInd);
	return q;
}

void DelayCalibration::loadParam(const char* pFilename){
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

		pElem=hRoot.FirstChild("DelayCalibrationSettings").Element();
		if (pElem){
			pElem->QueryDoubleAttribute("dt", &dt_);
			pElem->QueryDoubleAttribute("maxdelay", &maxDelay_);
		}

		pElem=hRoot.FirstChild("DelayCalibrationSettings").FirstChild("IMU").Element();
		if (pElem){
			pElem->QueryIntAttribute("use", &mInt);
			mbUseImu_ = mInt;
		}

		pElem=hRoot.FirstChild("DelayCalibrationSettings").FirstChild("Kinematic").Element();
		if (pElem){
			pElem->QueryIntAttribute("use", &mInt);
			mbUseEnc_ = mInt;
			pElem->QueryIntAttribute("difWindow", &difWindowEnc_);
		}

		pElem=hRoot.FirstChild("DelayCalibrationSettings").FirstChild("PoseSensor").Element();
		if (pElem){
			pElem->QueryIntAttribute("use", &mInt);
			mbUsePos_ = mInt;
			pElem->QueryIntAttribute("difWindow", &difWindowPos_);
		}
	}
}

}





