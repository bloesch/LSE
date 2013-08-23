/*!
* @file 	Manager.cpp
* @author 	Michael Blösch
* @date		10.10.2012
 */

#include "Manager.hpp"
#include "FilterOCEKF.hpp"
#include "FilterVUKF.hpp"
#include "FilterVUKF2.hpp"
#ifdef USE_CERES
#include "FilterInertialOF.hpp"
#endif
#include "DelayCalibration.hpp"
#include "tinyxml.h"
#include <iostream>
#include <vector>

namespace LSE {

Manager::Manager(const char* pFilename,Eigen::Vector3d (*f)(Eigen::Matrix<double,LSE_DOF_LEG,1>,int),Eigen::Matrix<double,3,LSE_DOF_LEG> (*J)(Eigen::Matrix<double,LSE_DOF_LEG,1>,int)):
legKin(f),legKinJac(J),g_(0.0,0.0,-9.81){
	tImu_ = 0;
	tEnc_ = 0;
	tPos_ = 0;

	// Init all parameters
	B_r_BI_.setZero();
	q_IB_ = Rotations::quatIdentity();
	B_r_BK_.setZero();
	q_KB_ = Rotations::quatIdentity();
	Rf_ = 0.002*Eigen::Matrix3d::Identity();
	Rw_ = 0.000873*Eigen::Matrix3d::Identity();
	Rs_ = 0.01*Eigen::Matrix3d::Identity();
	Ra_ = 0.001*Eigen::Matrix<double,LSE_DOF_LEG,LSE_DOF_LEG>::Identity();
	Rda_ = 0.01*Eigen::Matrix<double,LSE_DOF_LEG,LSE_DOF_LEG>::Identity();
	Rposr_ = 0.01*Eigen::Matrix<double,3,3>::Identity();
	Rposq_ = 0.01*Eigen::Matrix<double,3,3>::Identity();

	activeFilter_ = 0;
	loadParam(pFilename);

	// Initialize filter
	pFilterList_[0] = new FilterVUKF(this,pFilename);
	pFilterList_[1] = new FilterOCEKF(this,pFilename);
	pFilterList_[2] = new FilterVUKF2(this,pFilename);
	pDelayCalibration_ = new DelayCalibration(this,pFilename);
#ifdef USE_CERES
	pFilterList_[3] = new FilterInertialOF(this,pFilename);
	pRobotCalibration_ = new RobotCalibration(this,pFilename);
#endif

	// Logging Stuff
	isLogging_ = false;

	std::cout << "LSE Estimator ID: " << activeFilter_ << std::endl;

//	// Testing stuff
//	Eigen::Vector3d v;
//	v << 0.1, 0.2, 0.3;
//	OF::Expression<Eigen::Vector3d,double,double> x(v);
//	OF::Expression<Eigen::Vector3d,double,double> y(v);
//	OF::Expression<double,double,double> a(2);
//	OF::Expression<Eigen::Vector3d,double,Eigen::Vector3d> ax(&OF::ScalarVectorAddition,&a,&x);
//	OF::Expression<Eigen::Vector3d,Eigen::Vector3d,Eigen::Vector3d> z(&OF::VectorVectorAddition,&ax,&y);
//
//	z.fullEval();
//
//	OF::Expression<Eigen::Vector3d,double,Eigen::Vector3d> xx(&OF::ScalarVectorMult,&a,&z);
//	xx.fullEval();
//
//	OF::Expression<Eigen::Vector3d,Eigen::Vector3d,Eigen::Vector3d> xxx(&OF::CrossProduct,&xx,&x);
//	xxx.fullEval();
//
//	OF::Expression<double,Eigen::Vector3d,Eigen::Vector3d> xxxx(&OF::DotProduct,&xx,&x);
//	xxxx.fullEval();
//
//	std::cout << z.x_ << std::endl;
//	std::cout << xx.x_ << std::endl;
//	std::cout << xxx.x_ << std::endl;
//	std::cout << xxxx.x_ << std::endl;
//	std::cout << xxxx.op_->J1(xxxx.exp1_->x_,xxxx.exp2_->x_) << std::endl;
//	std::cout << xxxx.op_->J2(xxxx.exp1_->x_,xxxx.exp2_->x_) << std::endl;
//	std::cout << xxx.op_->J2(xxx.exp1_->x_,xxx.exp2_->x_) << std::endl;

//	// Test Ceres Stuff
////	  google::InitGoogleLogging(argv[0]);
//
//	  // The variable to solve for with its initial value. It will be
//	  // mutated in place by the solver.
//	  double x = 0.5;
//	  const double initial_x = x;
//
//	  // Build the problem.
//	  Problem problem;
//
//	  // Set up the only cost function (also known as residual). This uses
//	  // auto-differentiation to obtain the derivative (jacobian).
//	  CostFunction* cost_function = new AutoDiffCostFunction<CostFunctor, 1, 1>(new CostFunctor);
//	  problem.AddResidualBlock(cost_function, NULL, &x);
//
//	  // Run the solver!
//	  Solver::Options options;
//	  options.minimizer_progress_to_stdout = true;
//	  Solver::Summary summary;
//	  Solve(options, &problem, &summary);
//
//	  std::cout << summary.BriefReport() << "\n";
//	  std::cout << "x : " << initial_x << " -> " << x << "\n";

}

Manager::~Manager(){
	for(int i=0;i<NUM_FILTERS;i++){
		delete pFilterList_[i];
	}
	delete pDelayCalibration_;
#ifdef USE_CERES
	delete pRobotCalibration_;
#endif
}

void Manager::addImuMeas(const double& t,const ImuMeas& m){
	imuMeasList_.insert(std::pair<double,ImuMeas>(t,m));
	if(imuMeasList_.size() > LSE_MEAS_N){
		imuMeasList_.erase(imuMeasList_.begin());
	}
}

void Manager::addEncMeas(const double& t,const EncMeas& m){
	encMeasList_.insert(std::pair<double,EncMeas>(t,m));
	if(encMeasList_.size() > LSE_MEAS_N){
		encMeasList_.erase(encMeasList_.begin());
	}
}

void Manager::addPosMeas(const double& t,const PosMeas& m){
	posMeasList_.insert(std::pair<double,PosMeas>(t,m));
	if(posMeasList_.size() > LSE_MEAS_N){
		posMeasList_.erase(posMeasList_.begin());
	}
}

const ImuMeas* Manager::getImuMeas(double& t){
	std::map<double,ImuMeas>::iterator it;
	it = imuMeasList_.upper_bound(t);
	if(it == imuMeasList_.end()){
		return NULL;
	} else {
		t = (*it).first;
		return &(*it).second;
	}
}

const EncMeas* Manager::getEncMeas(double& t){
	std::map<double,EncMeas>::iterator it;
	it = encMeasList_.upper_bound(t);
	if(it == encMeasList_.end()){
		return NULL;
	} else {
		t = (*it).first;
		return &(*it).second;
	}
}

const PosMeas* Manager::getPosMeas(double& t){
	std::map<double,PosMeas>::iterator it;
	it = posMeasList_.upper_bound(t);
	if(it == posMeasList_.end()){
		return NULL;
	} else {
		t = (*it).first;
		return &(*it).second;
	}
}

void Manager::clearMeas(){
	imuMeasList_.clear();
	encMeasList_.clear();
	posMeasList_.clear();
}

void Manager::setImuTD(const double& TD){
	tImu_ = TD;
}
void Manager::setEncTD(const double& TD){
	tEnc_ = TD;
}
void Manager::setPosTD(const double& TD){
	tPos_ = TD;
}
double Manager::getImuTD(){
	return tImu_;
}
double Manager::getEncTD(){
	return tEnc_;
}
double Manager::getPosTD(){
	return tPos_;
}

void Manager::update(const double& t){
	pFilterList_[activeFilter_]->update(t);
}

void Manager::update(){
	pFilterList_[activeFilter_]->update();
}

State Manager::getEst(){
	return pFilterList_[activeFilter_]->getEst();
}

void Manager::resetEstimate(const double& t){
	pFilterList_[activeFilter_]->resetEstimate(t);
}

int Manager::delayIdentification(const double& t,const double& T){
	int res = pDelayCalibration_->calibrateDelay(t,T);
	return res;
}
#ifdef USE_CERES
int Manager::robotCalibration(const double& t,const double& T){
	int res = pRobotCalibration_->calibrateRobot(t,T);
	return res;
}
#endif

Eigen::Matrix3d Manager::gamma(const int& k,const Eigen::Vector3d& w,const double& dt){
	int b = k%2;
	int m = (k-b)/2;
	double wNorm = w.norm();
	double factor1 = 0;
	double factor2 = 0;
	Eigen::Matrix3d wk;
	Eigen::Matrix3d wk2;
	Eigen::Matrix3d G_k;

	// Get sqew matrices
	wk = Rotations::vecToSqew(w);
	wk2 = wk*wk;

	// Compute first factor
	if(wNorm*dt >= 1e-5*sqrt((2*m+3)*(2*m+4))){
		factor1 = cos(wNorm*dt);
		for(int i=0;i <= m;i++){
			factor1 += -pow(-1,i)*pow(wNorm*dt,2*i)/factorial(2*i);
		}
		factor1 *= pow(-1,m+1)/pow(wNorm,2+2*m);
	} else {
		factor1 = pow(dt,2*m+2)/factorial(2*m+2);
	}

	// Compute second factor
	if(wNorm*dt >= 1e-5*sqrt((2*m+2*b+2)*(2*m+2*b+3))){
		factor2 = sin(wNorm*dt);
		for(int i=0;i <= m+b-1;i++){
			factor2 += -pow(-1,i)*pow(wNorm*dt,2*i+1)/factorial(2*i+1);
		}
		factor2 *= pow(-1,m+b)/pow(wNorm,1+2*m+2*b);
	} else {
		factor2 = pow(dt,1+2*m+2*b)/factorial(1+2*m+2*b);
	}

	if(b==0){
		G_k = pow(dt,k)/factorial(k)*Eigen::Matrix3d::Identity()+factor1*wk2+factor2*wk;
	} else {
		G_k = pow(dt,k)/factorial(k)*Eigen::Matrix3d::Identity()+factor1*wk+factor2*wk2;
	}
	return G_k;
}

Eigen::Matrix3d Manager::gamma(const int& k,const Eigen::Vector3d& v){
	int b = k%2;
	int m = (k-b)/2;
	double vNorm = v.norm();
	double factor1 = 0;
	double factor2 = 0;
	Eigen::Matrix3d vk;
	Eigen::Matrix3d vk2;
	Eigen::Matrix3d G_k;

	// Get sqew matrices
	vk = Rotations::vecToSqew(v);
	vk2 = vk*vk;

	// Compute first factor
	if(vNorm >= 1e-5*sqrt((2*m+3)*(2*m+4))){
		factor1 = cos(vNorm);
		for(int i=0;i <= m;i++){
			factor1 += -pow(-1,i)*pow(vNorm,2*i)/factorial(2*i);
		}
		factor1 *= pow(-1,m+1)/pow(vNorm,2+2*m);
	} else {
		factor1 = 1/factorial(2*m+2);
	}

	// Compute second factor
	if(vNorm >= 1e-5*sqrt((2*m+2*b+2)*(2*m+2*b+3))){
		factor2 = sin(vNorm);
		for(int i=0;i <= m+b-1;i++){
			factor2 += -pow(-1,i)*pow(vNorm,2*i+1)/factorial(2*i+1);
		}
		factor2 *= pow(-1,m+b)/pow(vNorm,1+2*m+2*b);
	} else {
		factor2 = 1/factorial(1+2*m+2*b);
	}

	if(b==0){
		G_k = 1/factorial(k)*Eigen::Matrix3d::Identity()+factor1*vk2+factor2*vk;
	} else {
		G_k = 1/factorial(k)*Eigen::Matrix3d::Identity()+factor1*vk+factor2*vk2;
	}
	return G_k;
}

int Manager::factorial(const int& k){
	// Compute factorial of k
	int kFac = 1;
	for(int i=2;i<=k;i++){
		kFac = kFac*i;
	}
	return kFac;
}

void Manager::loadParam(const char* pFilename){
	// Open parameter file
	TiXmlDocument doc(pFilename);
	if (!doc.LoadFile()){
		std::cout << "No parameter file found at: " << pFilename << std::endl;
		return;
	}

	// Define handles and elements
	TiXmlHandle hDoc(&doc);
	TiXmlElement* pElem;
	TiXmlHandle hRoot(0);

	// Get root
	pElem=hDoc.FirstChildElement("LeggedStateEstimator").Element();
	if (pElem){
		pElem->QueryIntAttribute("activeFilter", &activeFilter_);
		hRoot=TiXmlHandle(pElem);

		pElem=hRoot.FirstChild("MeasurementSettings").FirstChild("Imu").FirstChild("AccelerometerStd").Element();
		if (pElem){
			pElem->QueryDoubleAttribute("x", &Rf_(0,0));
			pElem->QueryDoubleAttribute("y", &Rf_(1,1));
			pElem->QueryDoubleAttribute("z", &Rf_(2,2));
		}
		pElem=hRoot.FirstChild("MeasurementSettings").FirstChild("Imu").FirstChild("GyroscopeStd").Element();
		if (pElem){
			pElem->QueryDoubleAttribute("x", &Rw_(0,0));
			pElem->QueryDoubleAttribute("y", &Rw_(1,1));
			pElem->QueryDoubleAttribute("z", &Rw_(2,2));
		}
		pElem=hRoot.FirstChild("MeasurementSettings").FirstChild("Imu").FirstChild("TransOffset").Element();
		if (pElem){
			pElem->QueryDoubleAttribute("x", &B_r_BI_(0));
			pElem->QueryDoubleAttribute("y", &B_r_BI_(1));
			pElem->QueryDoubleAttribute("z", &B_r_BI_(2));
		}
		pElem=hRoot.FirstChild("MeasurementSettings").FirstChild("Imu").FirstChild("RotOffset").Element();
		if (pElem){
			pElem->QueryDoubleAttribute("x", &q_IB_(0));
			pElem->QueryDoubleAttribute("y", &q_IB_(1));
			pElem->QueryDoubleAttribute("z", &q_IB_(2));
			pElem->QueryDoubleAttribute("w", &q_IB_(3));
		}
		pElem=hRoot.FirstChild("MeasurementSettings").FirstChild("Kinematic").FirstChild("EncoderStd").Element();
		for(int i=0;i<LSE_DOF_LEG && pElem;i++){
			pElem->QueryDoubleAttribute("a", &Ra_(i,i));
			pElem->QueryDoubleAttribute("da", &Rda_(i,i));
			pElem = pElem->NextSiblingElement("EncoderStd");
		}
		pElem=hRoot.FirstChild("MeasurementSettings").FirstChild("Kinematic").FirstChild("ContactStd").Element();
		if (pElem){
			pElem->QueryDoubleAttribute("x", &Rs_(0,0));
			pElem->QueryDoubleAttribute("y", &Rs_(1,1));
			pElem->QueryDoubleAttribute("z", &Rs_(2,2));
		}
		pElem=hRoot.FirstChild("MeasurementSettings").FirstChild("Kinematic").FirstChild("TransOffset").Element();
		if (pElem){
			pElem->QueryDoubleAttribute("x", &B_r_BK_(0));
			pElem->QueryDoubleAttribute("y", &B_r_BK_(1));
			pElem->QueryDoubleAttribute("z", &B_r_BK_(2));
		}
		pElem=hRoot.FirstChild("MeasurementSettings").FirstChild("Kinematic").FirstChild("RotOffset").Element();
		if (pElem){
			pElem->QueryDoubleAttribute("x", &q_KB_(0));
			pElem->QueryDoubleAttribute("y", &q_KB_(1));
			pElem->QueryDoubleAttribute("z", &q_KB_(2));
			pElem->QueryDoubleAttribute("w", &q_KB_(3));
		}
		pElem=hRoot.FirstChild("MeasurementSettings").FirstChild("PoseSensor").FirstChild("PositionStd").Element();
		if (pElem){
			pElem->QueryDoubleAttribute("x", &Rposr_(0,0));
			pElem->QueryDoubleAttribute("y", &Rposr_(1,1));
			pElem->QueryDoubleAttribute("z", &Rposr_(2,2));
		}
		pElem=hRoot.FirstChild("MeasurementSettings").FirstChild("PoseSensor").FirstChild("AttituteStd").Element();
		if (pElem){
			pElem->QueryDoubleAttribute("x", &Rposq_(0,0));
			pElem->QueryDoubleAttribute("y", &Rposq_(1,1));
			pElem->QueryDoubleAttribute("z", &Rposq_(2,2));
		}
	}

	Rf_ = Rf_*Rf_;
	Rw_ = Rw_*Rw_;
	Rs_ = Rs_*Rs_;
	Ra_ = Ra_*Ra_;
	Rda_ = Rda_*Rda_;
	Rposr_ = Rposr_*Rposr_;
	Rposq_ = Rposq_*Rposq_;
}

void Manager::enableLogging(const char* pLogfile){
	if(isLogging_==false){
		isLogging_ = true;
//
//		std::string str;
//		str = pLogfile

		std::ostringstream oss (std::ostringstream::out);
		oss << pLogfile << "_F" << activeFilter_ << "_" << pFilterList_[activeFilter_]->getKeyString() << ".txt";

		ofsLog_.open(oss.str().c_str());
		ofsLog_.setf(std::ios_base::fixed);
		ofsLog_.precision(15);
	}
}

void Manager::disableLogging(){
	if(isLogging_==true){
		isLogging_ = false;
		ofsLog_.close();
	}
}

#ifdef USE_CERES
int Manager::getLengthOfBC(){
	return pRobotCalibration_->getN();
}
#endif

#ifdef USE_CERES
const RobotCalibration::state* Manager::getBCData(){
	return pRobotCalibration_->getBatch();
}
#endif

void Manager::addOflMeas(const double& t,const OflMeas& m){
	oflMeasList_.insert(std::pair<double,OflMeas>(t,m));
	if(oflMeasList_.size() > LSE_MEAS_N){
		oflMeasList_.erase(oflMeasList_.begin());
	}
}

const OflMeas* Manager::getOflMeas(double& t){
	std::map<double,OflMeas>::iterator it;
		it = oflMeasList_.upper_bound(t);
		if(it == oflMeasList_.end()){
			return NULL;
		} else {
			t = (*it).first;
			return &(*it).second;
		}
}

}

