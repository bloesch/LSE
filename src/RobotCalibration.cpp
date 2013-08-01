/*!
* @file 	RobotCalibration.cpp
* @author 	Michael Blösch
* @date		10.10.2012
 */

#include "RobotCalibration.hpp"
#include "Manager.hpp"
#include "tinyxml.h"

#include "ceres/ceres.h"
#include "ceres/rotation.h"

using ceres::AutoDiffCostFunction;
using ceres::CostFunction;
using ceres::Problem;
using ceres::Solver;
using ceres::Solve;

struct CostFunctor {
  template <typename T> bool operator()(const T* const x, T* residual) const {
    residual[0] = T(10.0) - x[0];
    return true;
  }
};

struct GyrResidual {
	GyrResidual(const double *w,double sigma):w_(w),sigma_(sigma){}

	template <typename T>
	bool operator()(const T* const x, const T* const b, T* residual) const {
		residual[0] = (T(w_[0]) - x[0] - b[0])/T(sigma_);
		residual[1] = (T(w_[1]) - x[1] - b[1])/T(sigma_);
		residual[2] = (T(w_[2]) - x[2] - b[2])/T(sigma_);
		return true;
	}

private:
	const double* w_;
	const double sigma_;
};

struct AccResidual {
	AccResidual(const double *f,double sigma):f_(f),sigma_(sigma){}

	template <typename T>
	bool operator()(const T* const q, const T* const a, const T* const b, T* residual) const {
		// f = C*(a-g)
		T vec[3];
		vec[0] = a[0];
		vec[1] = a[1];
		vec[2] = a[2] + 9.81;
		T IfI[3];
		ceres::UnitQuaternionRotatePoint(q,vec,IfI);
		residual[0] = (T(f_[0]) - IfI[0] - b[0])/T(sigma_);
		residual[1] = (T(f_[1]) - IfI[1] - b[1])/T(sigma_);
		residual[2] = (T(f_[2]) - IfI[2] - b[2])/T(sigma_);
		return true;
	}

private:
	const double* f_;
	const double sigma_;
};

struct PosResidual {
	PosResidual(double t,double sigma):t_(t),sigma_(sigma){}

	template <typename T>
	bool operator()(const T* const r1, const T* const r2, const T* const v, T* residual) const {
		// r2 = r1 + dt*v
//		T WaI[3];
//		ceres::UnitQuaternionRotatePoint(x+9,x+6,WaI);
		residual[0] = (r2[0] - r1[0] - T(t_)*v[0])/T(sigma_);
		residual[1] = (r2[1] - r1[1] - T(t_)*v[1])/T(sigma_);
		residual[2] = (r2[2] - r1[2] - T(t_)*v[2])/T(sigma_);
		return true;
	}

private:
	const double t_;
	const double sigma_;
};

struct VelResidual {
	VelResidual(double t,double sigma):t_(t),sigma_(sigma){}

	template <typename T>
	bool operator()(const T* const v1, const T* const v2, const T* const a, T* residual) const {
		// v2 = v1 + dt*a
		residual[0] = (v2[0] - v1[0] - T(t_)*a[0])/T(sigma_);
		residual[1] = (v2[1] - v1[1] - T(t_)*a[1])/T(sigma_);
		residual[2] = (v2[2] - v1[2] - T(t_)*a[2])/T(sigma_);
		return true;
	}

private:
	const double t_;
	const double sigma_;
};

struct RotResidual {
	RotResidual(double t,double sigma):t_(t),sigma_(sigma){}

	template <typename T>
	bool operator()(const T* const q1, const T* const q2, const T* const w, T* residual) const {
		// q2 = exp(dt*w)*q1
		T theta[3];
		theta[0] = T(t_)*w[0];
		theta[1] = T(t_)*w[1];
		theta[2] = T(t_)*w[2];
		T dq[4];
		ceres::AngleAxisToQuaternion(theta,dq);
		T q2p[4];
		ceres::QuaternionProduct(dq,q1,q2p);
		T q2pi[4];
		q2pi[0] = q2p[0];
		q2pi[1] = -q2p[1];
		q2pi[2] = -q2p[2];
		q2pi[3] = -q2p[3];
		T qres[4];
		ceres::QuaternionProduct(q2,q2pi,qres);
		T vec[3];
		ceres::QuaternionToAngleAxis(qres,vec);
		residual[0] = (vec[0])/T(sigma_);
		residual[1] = (vec[1])/T(sigma_);
		residual[2] = (vec[2])/T(sigma_);
		return true;
	}

private:
	const double t_;
	const double sigma_;
};

struct ForwardKinematics {
	ForwardKinematics(const double* alpha, int legIndex):alpha_(alpha),legIndex_(legIndex){}

	template <typename T>
	bool operator()(const T* const par, T* residual) const {
		T corAlpha[3];
		corAlpha[0] = T(alpha_[3*legIndex_+0])+par[3*legIndex_+5];
		corAlpha[1] = T(alpha_[3*legIndex_+1])+par[3*legIndex_+6];
		corAlpha[2] = T(alpha_[3*legIndex_+2])+par[3*legIndex_+7];
		residual[0] = T((legIndex_<2)*2-1)*par[0]+par[3]*sin(corAlpha[1])+par[4]*sin(corAlpha[1]+corAlpha[2]);
		residual[1] = -T((legIndex_%2)*2-1)*par[1]-sin(corAlpha[0])*(par[2]+par[3]*cos(corAlpha[1])+par[4]*cos(corAlpha[1]+corAlpha[2]));
		residual[2] = cos(corAlpha[0])*(par[2]+par[3]*cos(corAlpha[1])+par[4]*cos(corAlpha[1]+corAlpha[2]));
		return true;
	}

private:
	const double* alpha_;
	const int legIndex_;
};

struct KinResidual {
	KinResidual(const double* alpha, int legIndex):alpha_(alpha),legIndex_(legIndex),forwardKinematics_(alpha,legIndex){}

	template <typename T>
	bool operator()(const T* const r, const T* const q, const T* const p, const T* const I_r_IB, const T* const q_BI, const T* const par, T* residual) const {
		T vec1[3];
		vec1[0] = p[0]-r[0];
		vec1[1] = p[1]-r[1];
		vec1[2] = p[2]-r[2];
		T vec2[3];
		ceres::UnitQuaternionRotatePoint(q,vec1,vec2);
		vec1[0] = vec2[0]-I_r_IB[0];
		vec1[1] = vec2[1]-I_r_IB[1];
		vec1[2] = vec2[2]-I_r_IB[2];
		ceres::UnitQuaternionRotatePoint(q_BI,vec1,vec2);
		forwardKinematics_(par,vec1);
		residual[0] = vec2[0]-vec1[0];
		residual[1] = vec2[1]-vec1[1];
		residual[2] = vec2[2]-vec1[2];
		return true;
	}

private:
	const double* alpha_;
	const int legIndex_;
	ForwardKinematics forwardKinematics_;
};

struct PosMeasResidual {
	PosMeasResidual(const LSE::PosMeas* m, const double sigma1, const double sigma2):sigma1_(sigma1),sigma2_(sigma2){
		VrVC_[0] = m->r_(0);
		VrVC_[1] = m->r_(1);
		VrVC_[2] = m->r_(2);
		qCV_[0] = m->q_(3);
		qCV_[1] = -m->q_(0);
		qCV_[2] = -m->q_(1);
		qCV_[3] = -m->q_(2);
	}

	template <typename T>
	bool operator()(const T* const r, const T* const q, const T* const W_r_WV, const T* const q_VW, const T* const I_r_IC, const T* const q_CI, T* residual) const {
		T qWV[4];
		T qIV[4];
		T qCV[4];
		T qVCm[4];
		T qE[4];
		T vecErrRot[3];
		qWV[0] = q_VW[0];
		qWV[1] = -q_VW[1];
		qWV[2] = -q_VW[2];
		qWV[3] = -q_VW[3];
		ceres::QuaternionProduct(q,qWV,qIV);
		ceres::QuaternionProduct(q_CI,qIV,qCV);
		qVCm[0] = T(qCV_[0]);
		qVCm[1] = -T(qCV_[1]);
		qVCm[2] = -T(qCV_[2]);
		qVCm[3] = -T(qCV_[3]);
		ceres::QuaternionProduct(qCV,qVCm,qE);
		ceres::QuaternionToAngleAxis(qE,vecErrRot);
		residual[0] = vecErrRot[0]/T(sigma1_);
		residual[1] = vecErrRot[1]/T(sigma1_);
		residual[2] = vecErrRot[2]/T(sigma1_);

		T qWI[4];
		T W_r_IC[3];
		T V_r_VC[3];
		T W_r_VC[3];
		qWI[0] = q[0];
		qWI[1] = -q[1];
		qWI[2] = -q[2];
		qWI[3] = -q[3];
		ceres::UnitQuaternionRotatePoint(qWI,I_r_IC,W_r_IC);
		V_r_VC[0] = T(VrVC_[0]);
		V_r_VC[1] = T(VrVC_[1]);
		V_r_VC[2] = T(VrVC_[2]);
		ceres::UnitQuaternionRotatePoint(qWV,V_r_VC,W_r_VC);
		residual[3] = (W_r_VC[0] + W_r_WV[0] - r[0] - W_r_IC[0])/T(sigma2_);
		residual[4] = (W_r_VC[1] + W_r_WV[1] - r[1] - W_r_IC[1])/T(sigma2_);
		residual[5] = (W_r_VC[2] + W_r_WV[2] - r[2] - W_r_IC[2])/T(sigma2_);
		return true;
	}

private:
	double VrVC_[3];
	double qCV_[4];
	const double sigma1_;
	const double sigma2_;
};

namespace LSE {

RobotCalibration::RobotCalibration(Manager* pManager,const char* pFilename){
	pManager_ = pManager;

	// Timesteps and stuff
	t1_ = 0;
	t2_ = 0;

	// Parameters
	dt_ = 0;
	mbUseImu_ = 0;
	mbUseEnc_ = 0;
	mbUsePos_ = 0;

	loadParam(pFilename);

	// Initialize parameter blocks
	N_=1;
	PB_states_ = new state;

	// This is the square root of the final value!!!
	Wr_ = 0.01*Eigen::Matrix3d::Identity();
	Wv_ = 0.01*Eigen::Matrix3d::Identity();
	Wq_ = 0.01*Eigen::Matrix3d::Identity();
}

RobotCalibration::~RobotCalibration(){
	delete[] PB_states_;
}

int RobotCalibration::calibrateRobot(const double& t,const double& T){
	// Check conditions

	// Initialize and recheck
	if(initialize(t,T)==0) return 0;

	// Count number of states
	itImu1_ = pManager_->imuMeasList_.lower_bound(t1_-pManager_->tImu_);
	itImu2_ = pManager_->imuMeasList_.lower_bound(t2_-pManager_->tImu_);
	itImu_ =  itImu1_;
	N_ = 1;
	while(itImu_ != itImu2_){
		itImu_++;
		N_++;
	}
	std::cout << "Number of States: " << N_ << std::endl;

	// Delete and initialize new parameters
	delete[] PB_states_;
	PB_states_ = new state[N_];
	itImu_ =  itImu1_;
	for(int i=0;i<N_;i++){
		PB_states_[i].t_ = itImu_->first;
		itImu_++;
	}
	PB_IrIB_[0] = -pManager_->B_r_BI_(0);
	PB_IrIB_[1] = -pManager_->B_r_BI_(1);
	PB_IrIB_[2] = -pManager_->B_r_BI_(2);
	PB_qBI_[0] = 1;
	PB_qBI_[1] = 0;
	PB_qBI_[2] = 0;
	PB_qBI_[3] = 0;
	for(int i=0;i<12;i++){
		PB_p_[i] = 0;
	}
	PB_pkin_[0] = 0.2525;
	PB_pkin_[1] = 0.185;
	PB_pkin_[2] = -0.0685;
	PB_pkin_[3] = -0.2;
	PB_pkin_[4] = -0.235;
	for(int i=0;i<12;i++){
		PB_pkin_[i+5] = 0;
	}
	PB_WrWV_[0] = 0;
	PB_WrWV_[1] = 0;
	PB_WrWV_[2] = 0;
	PB_qVW_[0] = 1;
	PB_qVW_[1] = 0;
	PB_qVW_[2] = 0;
	PB_qVW_[3] = 0;
	PB_IrIC_[0] = 0;
	PB_IrIC_[1] = 0;
	PB_IrIC_[2] = 0;
	PB_qCI_[0] = 1;
	PB_qCI_[1] = 0;
	PB_qCI_[2] = 0;
	PB_qCI_[3] = 0;
	PB_bw_[0] = 0;
	PB_bw_[1] = 0;
	PB_bw_[2] = 0;
	PB_bf_[0] = 0;
	PB_bf_[1] = 0;
	PB_bf_[2] = 0;

	// Build the problem.
	Problem problem;

	CostFunction* cost_function;
	itImu_ =  itImu1_;
	const double* mpDouble;
	double sigma = 1;
	double dt = 0;
	for(int i=0;i<N_-1;i++){
		// Calculate timestep
		dt = PB_states_[i+1].t_-PB_states_[i].t_;
		if(dt<0.001) dt = 0.001;

		mpDouble = itImu_->second.w_.data();
		sigma = std::sqrt(pManager_->Rw_(0,0)/dt);
		cost_function = new AutoDiffCostFunction<GyrResidual, 3, 3, 3>(new GyrResidual(mpDouble,sigma));
		problem.AddResidualBlock(cost_function, NULL, PB_states_[i].w_, PB_bw_);
		mpDouble = itImu_->second.f_.data();
		sigma = std::sqrt(pManager_->Rf_(0,0)/dt);
		cost_function = new AutoDiffCostFunction<AccResidual, 3, 4, 3, 3>(new AccResidual(mpDouble,sigma));
		problem.AddResidualBlock(cost_function, NULL, PB_states_[i].q_, PB_states_[i].a_, PB_bf_);
		itImu_++;

		sigma = std::sqrt(Wr_(0,0)*dt);
		cost_function = new AutoDiffCostFunction<PosResidual, 3, 3, 3, 3>(new PosResidual(dt,sigma));
		problem.AddResidualBlock(cost_function, NULL, PB_states_[i].r_, PB_states_[i+1].r_, PB_states_[i].v_);
		sigma = std::sqrt(Wv_(0,0)*dt);
		cost_function = new AutoDiffCostFunction<VelResidual, 3, 3, 3, 3>(new VelResidual(dt,sigma));
		problem.AddResidualBlock(cost_function, NULL, PB_states_[i].v_, PB_states_[i+1].v_, PB_states_[i].a_);
		sigma = std::sqrt(Wq_(0,0)*dt);
		cost_function = new AutoDiffCostFunction<RotResidual, 3, 4, 4, 3>(new RotResidual(dt,sigma));
		problem.AddResidualBlock(cost_function, NULL, PB_states_[i].q_, PB_states_[i+1].q_, PB_states_[i].w_);
	}

    ceres::LocalParameterization* quaternion_parameterization = new ceres::QuaternionParameterization;
	for(int i=0;i<N_;i++){
		problem.SetParameterization(PB_states_[i].q_,quaternion_parameterization);
	}
	problem.SetParameterBlockConstant(PB_bw_);
	problem.SetParameterBlockConstant(PB_bf_);

	//Add kinematics to problem
	itImu_ =  itImu1_;
	itEnc_ = pManager_->encMeasList_.lower_bound(t1_-pManager_->tEnc_);
	for(int i=0;i<N_;i++){
		while(itEnc_->first < itImu_->first){
			mpDouble = itEnc_->second.e_.data();
			for(int j=0;j<4;j++){
				sigma = 1.0;
				cost_function = new AutoDiffCostFunction<KinResidual, 3, 3, 4, 3, 3, 4, 17>(new KinResidual(mpDouble,j));
				problem.AddResidualBlock(cost_function, NULL, PB_states_[i].r_, PB_states_[i].q_, &PB_p_[j*3], PB_IrIB_, PB_qBI_, PB_pkin_);
			}
			itEnc_++;
			if(itEnc_==pManager_->encMeasList_.end())
				break;
		}
		itImu_++;
	}

    problem.SetParameterization(PB_qBI_,quaternion_parameterization);
//    problem.SetParameterBlockConstant(PB_pkin_);
//    problem.SetParameterBlockConstant(PB_IrIB_);
//    problem.SetParameterBlockConstant(PB_qBI_);

	//Add 6DOF pose meas to problem
	ceres::LossFunction* loss_function = new ceres::HuberLoss(1.0);
	double sigma1 = 1;
	double sigma2 = 1;
	itImu_ =  itImu1_;
	itPos_ = pManager_->posMeasList_.lower_bound(t1_-pManager_->tPos_);
	for(int i=0;i<N_;i++){
		while(itPos_->first < itImu_->first){
			sigma1 = std::sqrt(pManager_->Rposq_(0,0));
			sigma2 = std::sqrt(pManager_->Rposr_(0,0));
			cost_function = new AutoDiffCostFunction<PosMeasResidual, 6, 3, 4, 3, 4, 3, 4>(new PosMeasResidual(&itPos_->second,sigma1,sigma2));
			problem.AddResidualBlock(cost_function, NULL, PB_states_[i].r_, PB_states_[i].q_, PB_WrWV_, PB_qVW_, PB_IrIC_, PB_qCI_);
			itPos_++;
			if(itPos_==pManager_->posMeasList_.end())
				break;
		}
		itImu_++;
	}

    problem.SetParameterization(PB_qVW_,quaternion_parameterization);
    problem.SetParameterization(PB_qCI_,quaternion_parameterization);
//    problem.SetParameterBlockConstant(PB_WrWV_);
//    problem.SetParameterBlockConstant(PB_IrIC_);
//    problem.SetParameterBlockConstant(PB_qVW_);
//    problem.SetParameterBlockConstant(PB_qCI_);

	// Run the solver!
	Solver::Options options;
	options.minimizer_progress_to_stdout = true;
	options.max_num_iterations = 50;
	options.minimizer_type = ceres::TRUST_REGION;
	options.trust_region_strategy_type = ceres::LEVENBERG_MARQUARDT;
//	options.trust_region_strategy_type = ceres::DOGLEG;
	options.dogleg_type = ceres::TRADITIONAL_DOGLEG;
	options.num_threads = 7;
	options.linear_solver_type = ceres::SPARSE_NORMAL_CHOLESKY;
	options.num_linear_solver_threads = 7;
	Solver::Summary summary;
	Solve(options, &problem, &summary);

	std::cout << summary.BriefReport() << "\n";
	std::cout << PB_WrWV_[0] << "\t" << PB_WrWV_[1] << "\t" << PB_WrWV_[2] << "\t" << PB_qVW_[0] << "\t" << PB_qVW_[1] << "\t" << PB_qVW_[2] << "\t" << PB_qVW_[3] << std::endl;
	std::cout << PB_IrIC_[0] << "\t" << PB_IrIC_[1] << "\t" << PB_IrIC_[2] << "\t" << PB_qCI_[0] << "\t" << PB_qCI_[1] << "\t" << PB_qCI_[2] << "\t" << PB_qCI_[3] << std::endl;
	std::cout << PB_IrIB_[0] << "\t" << PB_IrIB_[1] << "\t" << PB_IrIB_[2] << "\t" << PB_qBI_[0] << "\t" << PB_qBI_[1] << "\t" << PB_qBI_[2] << "\t" << PB_qBI_[3] << std::endl;
	std::cout << PB_pkin_[0] << "\t" << PB_pkin_[1] << "\t" << PB_pkin_[2] << "\t" << PB_pkin_[3] << "\t" << PB_pkin_[4] << std::endl;
	return 1;
}

int RobotCalibration::initialize(const double& t,const double& T){
	t1_ = t-T;
	t2_ = t;

	// Abort if no measurements
	if(mbUseImu_ && pManager_->imuMeasList_.empty())
		return 0;
	if(mbUseEnc_ && pManager_->encMeasList_.empty())
		return 0;
	if(mbUsePos_ && pManager_->posMeasList_.empty())
		return 0;
	if(dt_ == 0 && ! mbUseImu_)
		return 0;


	// Maximal range -> should avoid any subsequent problem
	if(dt_ == 0){
		itImu1_ = pManager_->imuMeasList_.lower_bound(t1_-pManager_->tImu_);
		if(itImu1_ == pManager_->imuMeasList_.end())
			return 0;
		itImu2_ = pManager_->imuMeasList_.upper_bound(t2_-pManager_->tImu_);
		if(itImu2_ == pManager_->imuMeasList_.begin())
			return 0;
		itImu2_--;
		if(mbUseEnc_){
			itImu_ = pManager_->imuMeasList_.lower_bound(pManager_->encMeasList_.begin()->first+pManager_->tEnc_-pManager_->tImu_);
			if(itImu_ != pManager_->imuMeasList_.begin())
				itImu_--;
			if(itImu1_->first < itImu_->first)
				itImu1_ = itImu_;
			itImu_ = pManager_->imuMeasList_.upper_bound(pManager_->encMeasList_.rbegin()->first+pManager_->tEnc_-pManager_->tImu_);
			if(itImu_ == pManager_->imuMeasList_.end())
				itImu_--;
			if(itImu2_->first > itImu_->first)
				itImu2_ = itImu_;
			if(itImu1_ == itImu2_)
				return 0;
		}
		if(mbUsePos_){
			itImu_ = pManager_->imuMeasList_.lower_bound(pManager_->posMeasList_.begin()->first+pManager_->tPos_-pManager_->tImu_);
			if(itImu_ != pManager_->imuMeasList_.begin())
				itImu_--;
			if(itImu1_->first < itImu_->first)
				itImu1_ = itImu_;
			itImu_ = pManager_->imuMeasList_.upper_bound(pManager_->posMeasList_.rbegin()->first+pManager_->tPos_-pManager_->tImu_);
			if(itImu_ == pManager_->imuMeasList_.end())
				itImu_--;
			if(itImu2_->first > itImu_->first)
				itImu2_ = itImu_;
			if(itImu1_ == itImu2_)
				return 0;
		}
	} else {
		return 0; // TODO change
	}

	t1_ = itImu1_->first;
	t2_ = itImu2_->first;

	return 1;
}

int RobotCalibration::getN(){
	return N_;
}

const RobotCalibration::state* RobotCalibration::getBatch(){
	return PB_states_;
}


void RobotCalibration::loadParam(const char* pFilename){
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

		pElem=hRoot.FirstChild("RobotCalibrationSettings").Element();
		if (pElem){
			pElem->QueryDoubleAttribute("dt", &dt_);
		}

		pElem=hRoot.FirstChild("RobotCalibrationSettings").FirstChild("IMU").Element();
		if (pElem){
			mInt = 1;
			pElem->QueryIntAttribute("use", &mInt);
			mbUseImu_ = mInt;
		}

		pElem=hRoot.FirstChild("RobotCalibrationSettings").FirstChild("Kinematic").Element();
		if (pElem){
			mInt = 1;
			pElem->QueryIntAttribute("use", &mInt);
			mbUseEnc_ = mInt;
		}

		pElem=hRoot.FirstChild("RobotCalibrationSettings").FirstChild("PoseSensor").Element();
		if (pElem){
			mInt = 1;
			pElem->QueryIntAttribute("use", &mInt);
			mbUsePos_ = mInt;
		}

		pElem=hRoot.FirstChild("RobotCalibrationSettings").FirstChild("Position").FirstChild("PreNoiStd").Element();
		if (pElem){
			pElem->QueryDoubleAttribute("x", &Wr_(0,0));
			pElem->QueryDoubleAttribute("y", &Wr_(1,1));
			pElem->QueryDoubleAttribute("z", &Wr_(2,2));
		}

		pElem=hRoot.FirstChild("RobotCalibrationSettings").FirstChild("Velocity").FirstChild("PreNoiStd").Element();
		if (pElem){
			pElem->QueryDoubleAttribute("x", &Wv_(0,0));
			pElem->QueryDoubleAttribute("y", &Wv_(1,1));
			pElem->QueryDoubleAttribute("z", &Wv_(2,2));
		}

		pElem=hRoot.FirstChild("RobotCalibrationSettings").FirstChild("Attitude").FirstChild("PreNoiStd").Element();
		if (pElem){
			pElem->QueryDoubleAttribute("x", &Wq_(0,0));
			pElem->QueryDoubleAttribute("y", &Wq_(1,1));
			pElem->QueryDoubleAttribute("z", &Wq_(2,2));
		}
	}

	Wr_ = Wr_*Wr_;
	Wv_ = Wv_*Wv_;
	Wq_ = Wq_*Wq_;
}

}





