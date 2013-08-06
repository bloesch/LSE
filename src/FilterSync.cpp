/*!
* @file 	FilterSync.cpp
* @author 	Michael Blösch
* @date		10.10.2012
 */

#include "FilterSync.hpp"
#include "Manager.hpp"
#include "tinyxml.h"
#include <map>

namespace LSE {

FilterSync::FilterSync(Manager* pManager,const char* pFilename){
	pManager_ = pManager;

	// Init all parameters
	xInit_.t_ = 0;
	xInit_.r_.setZero();
	xInit_.v_.setZero();
	xInit_.q_ = Rotations::quatIdentity();
	xInit_.w_.setZero();
	xInit_.p_.setZero();
	xInit_.bf_.setZero();
	xInit_.bw_.setZero();
	for(int i=0;i<LSE_N_LEG;i++){
		xInit_.CFC_[i] = 0;
	}
	xInit_.P_.setIdentity();
	xInit_.rLast_.setZero();
	xInit_.vLast_.setZero();
	xInit_.qLast_ = Rotations::quatIdentity();
	xInit_.pLin_.setZero();
	xInit_.wLin_.setZero();
	xInit_.f1Lin_ = -pManager_->g_;
	xInit_.f2Lin_ = -pManager_->g_;
	xInit_.tLast_ = 0;
	xInit_.f_ = -pManager_->g_;
	Wr_ = 0*Eigen::Matrix3d::Identity();
	Wv_ = 0.003*Eigen::Matrix3d::Identity();
	Wq_ = 0*Eigen::Matrix3d::Identity();
	Wp_ = 0.01*Eigen::Matrix3d::Identity();
	Wbf_ = 0.0001*Eigen::Matrix3d::Identity();
	Wbw_ = 0.000618*Eigen::Matrix3d::Identity();

	// Flags
	mbEstimateRotBias_ = true;
	mbEstimateAccBias_ = true;
	mbUseImu_ = true;
	mbUseKin_ = true;
	mbAssumeFlatFloor_ = false;

	loadParam(pFilename);
}

FilterSync::~FilterSync(){
}

void FilterSync::update(const double& t){
	// Find safe time and filter safe state
	double tsNew = t;
	if(!pManager_->imuMeasList_.empty() && !pManager_->encMeasList_.empty()){
		tsNew = std::min(tsNew,pManager_->imuMeasList_.rbegin()->first+pManager_->tImu_);
		tsNew = std::min(tsNew,pManager_->encMeasList_.rbegin()->first+pManager_->tEnc_);
		if(xs_.t_<tsNew){
			filterState(xs_,tsNew);
		}
	}

	// Predict state
	xp_ = xs_;
	filterState(xp_,t);
}

void FilterSync::update(){
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

State FilterSync::getEst(){
	Eigen::Matrix3d R_WI,R_IB;
	R_WI = Rotations::quatToRotMat(xp_.q_).transpose();
	R_IB = Rotations::quatToRotMat(pManager_->q_IB_);
	State x = State();
	x.t_ = xp_.t_;
	x.r_ = xp_.r_-R_WI*R_IB*pManager_->B_r_BI_;
	x.v_ = xp_.v_-R_WI*(Rotations::vecToSqew(xp_.w_)*R_IB*pManager_->B_r_BI_);
	x.q_ = Rotations::quatL(pManager_->q_IB_).transpose()*xp_.q_;
	x.w_ = R_IB.transpose()*xp_.w_;
	x.P_.setZero();
	x.P_.block(0,0,9,9) = xp_.P_.block(0,0,9,9);
	x.P_.block(9,9,3,3) = xp_.P_.block(12,12,3,3)+pManager_->Rw_;
	return x;
}

void FilterSync::resetEstimate(const double& t){
	xs_ = xInit_;
	xs_.t_ = t;
	xs_.tLast_ = t;
	xp_ = xs_;
}

void FilterSync::filterState(InternState& x,const double& tEnd){
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
		if(itEnc != pManager_->encMeasList_.end()) tNext = std::min(tEnd,itEnc->first+pManager_->tEnc_);

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

void FilterSync::predictState(InternState& x, const double& tPre, const ImuMeas& m){
	double dt = tPre-x.t_;
	if(!mbEstimateAccBias_) x.bf_.setZero();
	if(!mbEstimateRotBias_) x.bw_.setZero();
	Eigen::Vector3d w = m.w_-x.bw_;
	if(!mbUseImu_) w.setZero();
	Eigen::Vector3d a = pManager_->g_;
	a = a + 0.5*Rotations::quatToRotMat(x.q_).transpose()*(m.f_-x.bf_);
	x.q_ = Rotations::quatL(Rotations::rotVecToQuat(-dt*w))*x.q_;
	a = a + 0.5*Rotations::quatToRotMat(x.q_).transpose()*(m.f_-x.bf_);
	if(!mbUseImu_){
		a.setZero();
	}
	x.r_ = x.r_ + dt*x.v_ + std::pow(dt,2)/2*a;
	x.v_ = x.v_ + dt*a;
	x.w_ = w;
	x.t_ = tPre;
	x.f_ = m.f_-x.bf_;
}

void FilterSync::encUpdateState(InternState& x, const EncMeas& m){
	// Compute forward kinematic
	Eigen::Matrix<double,3,LSE_N_LEG> s;
	for(int i=0;i<LSE_N_LEG;i++){
		// I_r_IF = C(q_IB)*(-B_r_BI + B_r_BK + C'(q_KB)*K_r_KF
		s.col(i) = Rotations::quatToRotMat(pManager_->q_IB_)*(-pManager_->B_r_BI_+pManager_->B_r_BK_+Rotations::quatToRotMat(pManager_->q_KB_).transpose()*(*pManager_->legKin)(m.e_.col(i),i));
	}

	// Update Contact count
	for(int i=0;i<LSE_N_LEG;i++){
		if(m.CF_[i]){
			x.CFC_[i]++;
		} else {
			x.CFC_[i] = 0;
		}
	}

	// Required rotation matrices
	Eigen::Matrix3d R_IW,R_WI,R_IW_last,R_WI_last;
	R_IW = Rotations::quatToRotMat(x.q_);
	R_WI = R_IW.transpose();
	R_IW_last = Rotations::quatToRotMat(x.qLast_);
	R_WI_last = R_IW_last.transpose();

	// Handle Initialization of new contacts
	for(int i=0;i<LSE_N_LEG;i++){
		if(x.CFC_[i] == 1){
			x.p_.col(i) = x.r_+R_WI*s.col(i);
			if(mbAssumeFlatFloor_){
				x.p_(2,i) = 0;
			}
			x.P_.block(0,15+3*i,15+3*LSE_N_LEG,3).setZero();
			x.P_.block(15+3*i,0,3,15+3*LSE_N_LEG).setZero();
			x.P_.block(15+3*i,15+3*i,3,3) = Eigen::Matrix3d::Identity()*1e10;
		}
	}

	// Evaluate linearization point
	for(int i=0;i<LSE_N_LEG;i++){
		if(x.CFC_[i] == 1){
			x.pLin_.col(i) = x.p_.col(i);
		}
	}
	double dt = x.t_-x.tLast_;
	if(dt>1e-5){
		x.wLin_ = -1/dt*Rotations::quatToRotVec(Rotations::quatL(x.q_)*Rotations::quatInverse(x.qLast_));
		x.f1Lin_ = 2/pow(dt,2)*R_IW_last*(x.r_-x.rLast_-dt*x.vLast_-0.5*pow(dt,2)*pManager_->g_);
		x.f2Lin_ = 1/dt*R_IW_last*(x.v_-x.vLast_-dt*pManager_->g_);
	} else {
		x.wLin_ = x.w_;
		x.f1Lin_ = x.f_;
		x.f2Lin_ = x.f_;
	}
	x.tLast_ = x.t_;
	x.rLast_ = x.r_;
	x.vLast_ = x.v_;
	x.qLast_ = x.q_;

	// Get Gamma Matrices
	Eigen::Matrix3d G0,G1,G2,G3;
	G0 = pManager_->gamma(0,x.wLin_,dt);
	G1 = pManager_->gamma(1,x.wLin_,dt);
	G2 = pManager_->gamma(2,x.wLin_,dt);
	G3 = pManager_->gamma(3,x.wLin_,dt);

	// Predict covariance
	Eigen::Matrix<double,15+3*LSE_N_LEG,15+3*LSE_N_LEG> F;
	F.setIdentity();
	F.block(0,3,3,3) = dt*Eigen::Matrix3d::Identity();
	F.block(0,6,3,3) = -R_WI_last*Rotations::vecToSqew(pow(dt,2)/2.0*x.f1Lin_);
	F.block(0,9,3,3) = -R_WI_last*pow(dt,2)/2.0;
	F.block(3,6,3,3) = -R_WI_last*Rotations::vecToSqew(dt*x.f2Lin_);
	F.block(3,9,3,3) = -R_WI_last*dt;
	F.block(6,6,3,3) = G0.transpose();
	F.block(6,12,3,3) = -G1.transpose();

	if(!mbEstimateAccBias_){
		F.block(0,9,3,3).setZero();
		F.block(3,9,3,3).setZero();
	}
	if(!mbEstimateRotBias_){
		F.block(6,12,3,3).setZero();
	}

	Eigen::Matrix<double,15+3*LSE_N_LEG,15+3*LSE_N_LEG> W;
	W.setZero();
	if(mbEstimateAccBias_){
		W.block(0,0,3,3) = pow(dt,3)/3.0*pManager_->Rf_+pow(dt,5)/20.0*Wbf_+dt*Wr_;
		W.block(0,3,3,3) = pow(dt,2)/2.0*pManager_->Rf_+pow(dt,4)/8.0*Wbf_;
		W.block(0,9,3,3) = -pow(dt,3)/6.0*R_WI_last*Wbf_;
		W.block(3,0,3,3) = W.block(0,3,3,3).transpose();
		W.block(3,3,3,3) = pow(dt,1)/1.0*pManager_->Rf_+pow(dt,3)/3.0*Wbf_+dt*Wv_;
		W.block(3,9,3,3) = -pow(dt,2)/2.0*R_WI_last*Wbf_;
		W.block(9,0,3,3) = W.block(0,9,3,3).transpose();
		W.block(9,3,3,3) = W.block(3,9,3,3).transpose();
		W.block(9,9,3,3) = dt*Wbf_;
	} else {
		W.block(0,0,3,3) = pow(dt,3)/3.0*pManager_->Rf_+dt*Wr_;
		W.block(0,3,3,3) = pow(dt,2)/2.0*pManager_->Rf_;
		W.block(3,0,3,3) = W.block(0,3,3,3).transpose();
		W.block(3,3,3,3) = pow(dt,1)/1.0*pManager_->Rf_+dt*Wv_;
	}
	if(mbEstimateRotBias_){
		W.block(6,6,3,3) = pow(dt,1)/1.0*pManager_->Rw_+(G3+G3.transpose())*Wbw_+dt*Wq_;
		W.block(6,12,3,3) = -G2.transpose()*Wbw_;
		W.block(12,6,3,3) = W.block(6,12,3,3).transpose();
		W.block(12,12,3,3) = dt*Wbw_;
	} else {
		W.block(6,6,3,3) = pow(dt,1)/1.0*pManager_->Rw_+dt*Wq_;
	}
	for(int i=0;i<LSE_N_LEG;i++){
		W.block(15+i*3,15+i*3,3,3) = dt*R_WI_last*Wp_*R_IW_last;
	}

	x.P_ = F*x.P_*F.transpose()+W;

	// Update state and covariance
	Eigen::Matrix<double,15+3*LSE_N_LEG,1> dx;
	int contactNo = 0;
	for(int i=0;i < LSE_N_LEG;i++){
		if(x.CFC_[i]>0){
			contactNo++;
		}
	}
	if(contactNo>0){
		// Init matrices
		Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> G(3*contactNo,15+3*LSE_N_LEG);
		Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> R(3*contactNo,3*contactNo);
		Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> S(3*contactNo,3*contactNo);
		Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> K(15+3*LSE_N_LEG,3*contactNo);
		Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> y(3*contactNo,1);
		G.setZero();
		R.setZero();
		S.setZero();
		K.setZero();
		y.setZero();

		int j=0;
		for(int i=0;i < LSE_N_LEG;i++){
			if(x.CFC_[i]>0){
				G.block(3*j,0,3,3) = -R_IW;
				G.block(3*j,6,3,3) = R_IW*Rotations::vecToSqew(x.pLin_.col(i)-x.r_)*R_WI;
				G.block(3*j,15+3*i,3,3) = R_IW;
				if(mbAssumeFlatFloor_){
					G.block(3*j,17+3*i,3,1).setZero();
				}
				if(pManager_->legKinJac!=NULL){
					R.block(3*j,3*j,3,3) = pManager_->Rs_ + (*pManager_->legKinJac)(m.e_.col(i),i)*pManager_->Ra_*((*pManager_->legKinJac)(m.e_.col(i),i)).transpose();
				} else {
					R.block(3*j,3*j,3,3) = pManager_->Rs_;
				}
				y.block(3*j,0,3,1) = s.col(i)-R_IW*(x.p_.col(i)-x.r_);
				j++;
			}
		}

		S = G*x.P_*G.transpose() + R;
		K = x.P_*G.transpose()*S.inverse();
		// Correction vector
		dx = K*y;
		x.P_ = (Eigen::Matrix<double,15+3*LSE_N_LEG,15+3*LSE_N_LEG>::Identity() - K*G)*x.P_;

		// Get corrected state
		x.r_ = x.r_ + dx.block(0,0,3,1);
		x.v_ = x.v_ + dx.block(3,0,3,1);
		x.q_ = Rotations::quatL(Rotations::rotVecToQuat(-dx.block(6,0,3,1)))*x.q_;
		x.q_.normalize();
		if(mbEstimateAccBias_) x.bf_ = x.bf_ + dx.block(9,0,3,1);
		if(mbEstimateRotBias_) x.bw_ = x.bw_ + dx.block(12,0,3,1);
		for(int i=0;i < LSE_N_LEG;i++){
			if(x.CFC_[i]>0){
				x.p_.col(i) = x.p_.col(i) + dx.block(15+3*i,0,3,1);
				if(mbAssumeFlatFloor_){
					x.p_(2,i) = 0;
				}
			}
		}
	}
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
			pElem->QueryDoubleAttribute("x", &xInit_.r_(0));
			pElem->QueryDoubleAttribute("y", &xInit_.r_(1));
			pElem->QueryDoubleAttribute("z", &xInit_.r_(2));
		}
		pElem=hRoot.FirstChild("SyncSettings").FirstChild("Position").FirstChild("InitStd").Element();
		if (pElem){
			pElem->QueryDoubleAttribute("x", &xInit_.P_(0,0));
			pElem->QueryDoubleAttribute("y", &xInit_.P_(1,1));
			pElem->QueryDoubleAttribute("z", &xInit_.P_(2,2));
		}

		pElem=hRoot.FirstChild("SyncSettings").FirstChild("Velocity").FirstChild("Init").Element();
		if (pElem){
			pElem->QueryDoubleAttribute("x", &xInit_.v_(0));
			pElem->QueryDoubleAttribute("y", &xInit_.v_(1));
			pElem->QueryDoubleAttribute("z", &xInit_.v_(2));
		}
		pElem=hRoot.FirstChild("SyncSettings").FirstChild("Velocity").FirstChild("InitStd").Element();
		if (pElem){
			pElem->QueryDoubleAttribute("x", &xInit_.P_(3,3));
			pElem->QueryDoubleAttribute("y", &xInit_.P_(4,4));
			pElem->QueryDoubleAttribute("z", &xInit_.P_(5,5));
		}

		pElem=hRoot.FirstChild("SyncSettings").FirstChild("Attitude").FirstChild("Init").Element();
		if (pElem){
			pElem->QueryDoubleAttribute("x", &xInit_.q_(0));
			pElem->QueryDoubleAttribute("y", &xInit_.q_(1));
			pElem->QueryDoubleAttribute("z", &xInit_.q_(2));
			pElem->QueryDoubleAttribute("w", &xInit_.q_(3));
		}
		pElem=hRoot.FirstChild("SyncSettings").FirstChild("Attitude").FirstChild("InitStd").Element();
		if (pElem){
			pElem->QueryDoubleAttribute("x", &xInit_.P_(6,6));
			pElem->QueryDoubleAttribute("y", &xInit_.P_(7,7));
			pElem->QueryDoubleAttribute("z", &xInit_.P_(8,8));
		}

		pElem=hRoot.FirstChild("SyncSettings").FirstChild("AccelerometerBias").FirstChild("Init").Element();
		if (pElem){
			pElem->QueryDoubleAttribute("x", &xInit_.bf_(0));
			pElem->QueryDoubleAttribute("y", &xInit_.bf_(1));
			pElem->QueryDoubleAttribute("z", &xInit_.bf_(2));
		}
		pElem=hRoot.FirstChild("SyncSettings").FirstChild("AccelerometerBias").FirstChild("InitStd").Element();
		if (pElem){
			pElem->QueryDoubleAttribute("x", &xInit_.P_(9,9));
			pElem->QueryDoubleAttribute("y", &xInit_.P_(10,10));
			pElem->QueryDoubleAttribute("z", &xInit_.P_(11,11));
		}

		pElem=hRoot.FirstChild("SyncSettings").FirstChild("GyroscopeBias").FirstChild("Init").Element();
		if (pElem){
			pElem->QueryDoubleAttribute("x", &xInit_.bw_(0));
			pElem->QueryDoubleAttribute("y", &xInit_.bw_(1));
			pElem->QueryDoubleAttribute("z", &xInit_.bw_(2));
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
			pElem->QueryIntAttribute("useImu", &mInt);
			mbUseImu_ = mInt;
			pElem->QueryIntAttribute("useKin", &mInt);
			mbUseKin_ = mInt;
		}
		pElem=hRoot.FirstChild("SyncSettings").FirstChild("Foothold").Element();
		if (pElem){
			pElem->QueryIntAttribute("assumeFlatFloor", &mInt);
			mbAssumeFlatFloor_ = mInt;
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

	xInit_.rLast_ = xInit_.r_;
	xInit_.vLast_ = xInit_.v_;
	xInit_.qLast_ = xInit_.q_;
	xInit_.P_ = xInit_.P_*xInit_.P_;
	Wr_ = Wr_*Wr_;
	Wv_ = Wv_*Wv_;
	Wq_ = Wq_*Wq_;
	Wp_ = Wp_*Wp_;
	Wbf_ = Wbf_*Wbf_;
	Wbw_ = Wbw_*Wbw_;
}

}




