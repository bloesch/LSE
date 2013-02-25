/*!
* @file 	PythonManager.cpp
* @author 	Michael Blösch
* @date		10.10.2012
 */

#include "PythonManager.hpp"
#include <Eigen/Dense>

// Kinematics of robot
Eigen::Vector3d legKin(Eigen::Matrix<double,LSE_DOF_LEG,1> a,int i){
	// Todo: include the parameters into the function call -> later used for robot calibration
	double bx = 0.2525;
	double by = 0.185;
	double lH = -0.0685;
	double lT = -0.2;
	double lS = -0.235;

	Eigen::Vector3d s;
	s(0) = ((i<2)*2-1)*bx+lT*sin(a(1))+lS*sin(a(1)+a(2));
	s(1) = -((i%2)*2-1)*by-sin(a(0))*(lH+lT*cos(a(1))+lS*cos(a(1)+a(2)));
	s(2) = cos(a(0))*(lH+lT*cos(a(1))+lS*cos(a(1)+a(2)));
	return s;
}

// Kinematic Jacobian of robot
Eigen::Matrix<double,3,LSE_DOF_LEG> legKinJac(Eigen::Matrix<double,LSE_DOF_LEG,1> a,int i){
	double lH = -0.0685;
	double lT = -0.2;
	double lS = -0.235;

	Eigen::Matrix<double,3,LSE_DOF_LEG> J;
	J.setZero();
	J(0,1) = lS*cos(a(1)+a(2))+lT*cos(a(1));
	J(0,2) = lS*cos(a(1)+a(2));
	J(1,0) = -cos(a(0))*(lH+lT*cos(a(1))+lS*cos(a(1)+a(2)));
	J(1,1) = sin(a(0))*(lT*sin(a(1))+lS*sin(a(1)+a(2)));
	J(1,2) = lS*sin(a(0))*sin(a(1)+a(2));
	J(2,0) = -sin(a(0))*(lH+lT*cos(a(1))+lS*cos(a(1)+a(2)));
	J(2,1) = -cos(a(0))*(lT*sin(a(1))+lS*sin(a(1)+a(2)));
	J(2,2) = -lS*cos(a(0))*sin(a(1)+a(2));
	return J;
}

namespace LSE {

PythonManager::PythonManager(){
	pManager_ = new LSE::Manager("Parameters.xml",&legKin,&legKinJac);
}

PythonManager::~PythonManager(){
	delete pManager_;
}

//void PythonManager::addImuMeas(const double& t,const ImuMeas& m){
//	imuMeasList_.insert(std::pair<double,ImuMeas>(t,m));
//	if(imuMeasList_.size() > LSE_MEAS_N){
//		imuMeasList_.erase(imuMeasList_.begin());
//	}
//}
//
//void PythonManager::addEncMeas(const double& t,const EncMeas& m){
//	encMeasList_.insert(std::pair<double,EncMeas>(t,m));
//	if(encMeasList_.size() > LSE_MEAS_N){
//		encMeasList_.erase(encMeasList_.begin());
//	}
//}
//
//void PythonManager::addPosMeas(const double& t,const PosMeas& m){
//	posMeasList_.insert(std::pair<double,PosMeas>(t,m));
//	if(posMeasList_.size() > LSE_MEAS_N){
//		posMeasList_.erase(posMeasList_.begin());
//	}
//}
//
//const ImuMeas* PythonManager::getImuMeas(double& t){
//	std::map<double,ImuMeas>::iterator it;
//	it = imuMeasList_.upper_bound(t);
//	if(it == imuMeasList_.end()){
//		return NULL;
//	} else {
//		t = (*it).first;
//		return &(*it).second;
//	}
//}
//
//const EncMeas* PythonManager::getEncMeas(double& t){
//	std::map<double,EncMeas>::iterator it;
//	it = encMeasList_.upper_bound(t);
//	if(it == encMeasList_.end()){
//		return NULL;
//	} else {
//		t = (*it).first;
//		return &(*it).second;
//	}
//}
//
//const PosMeas* PythonManager::getPosMeas(double& t){
//	std::map<double,PosMeas>::iterator it;
//	it = posMeasList_.upper_bound(t);
//	if(it == posMeasList_.end()){
//		return NULL;
//	} else {
//		t = (*it).first;
//		return &(*it).second;
//	}
//}
//
//void PythonManager::setEncTD(const double& TD){
//	tEnc_ = TD;
//}
//void PythonManager::setPosTD(const double& TD){
//	tPos_ = TD;
//}
//double PythonManager::getEncTD(){
//	return tEnc_;
//}
//double PythonManager::getPosTD(){
//	return tPos_;
//}
//
//void PythonManager::update(const double& t){
//	pFilterList_[activeFilter_]->update(t);
//}
//
//void PythonManager::update(){
//	pFilterList_[activeFilter_]->update();
//}
//
//State PythonManager::getEst(){
//	return pFilterList_[activeFilter_]->getEst();
//}
//
//void PythonManager::resetEstimate(const double& t){
//	pFilterList_[activeFilter_]->resetEstimate(t);
//}
//
//int PythonManager::delayIdentification(const double& t,const double& T){
//	int res = pDelayCalibration_->calibrateDelay(t,T);
//	return res;
//}
//
#if WRAP_PYTHON
void PythonManager::update_pythont(double t){
	return pManager_->update(t);
}
void PythonManager::update_python(){
	return pManager_->update();
}
void PythonManager::getEst_python(){
	State x;
	x = pManager_->getEst();
	return;
}
void PythonManager::resetEstimate_python(double t){
	return pManager_->resetEstimate(t);
}
int PythonManager::delayIdentification_python(double t, double T){
	return pManager_->delayIdentification(t,T);
}
void PythonManager::setImuTD_python(double TD){
	return pManager_->setImuTD(TD);
}
void PythonManager::setEncTD_python(double TD){
	return pManager_->setEncTD(TD);
}
void PythonManager::setPosTD_python(double TD){
	return pManager_->setPosTD(TD);
}
double PythonManager::getImuTD_python(){
	return pManager_->getImuTD();
}
double PythonManager::getEncTD_python(){
	return pManager_->getEncTD();
}
double PythonManager::getPosTD_python(){
	return pManager_->getPosTD();
}
using namespace boost::python;
BOOST_PYTHON_MODULE(_PythonManager)
{
    class_<PythonManager>("PythonManager")
        .def("update", &PythonManager::update_pythont)
        .def("update", &PythonManager::update_python)
        .def("getEst", &PythonManager::getEst_python)
        .def("resetEstimate", &PythonManager::resetEstimate_python)
        .def("delayIdentification", &PythonManager::delayIdentification_python)
        .def("setImuTD", &PythonManager::setImuTD_python)
        .def("setEncTD", &PythonManager::setEncTD_python)
        .def("setPosTD", &PythonManager::setPosTD_python)
        .def("getImuTD", &PythonManager::getImuTD_python)
        .def("getEncTD", &PythonManager::getEncTD_python)
        .def("getPosTD", &PythonManager::getPosTD_python)
    ;
}
#endif

}
