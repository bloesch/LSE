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

PythonManager::PythonManager(std::string filename){
	pManager_ = new LSE::Manager(filename.c_str(),&legKin,&legKinJac);
}

PythonManager::~PythonManager(){
	delete pManager_;
}

#if WRAP_PYTHON
void PythonManager::addImuMeas_python(double t, PyObject* pyf, PyObject* pyw){
	ImuMeas imuMeas;
	imuMeas.f_(0) = ((double*)PyArray_DATA(pyf))[0];
	imuMeas.f_(1) = ((double*)PyArray_DATA(pyf))[1];
	imuMeas.f_(2) = ((double*)PyArray_DATA(pyf))[2];
	imuMeas.w_(0) = ((double*)PyArray_DATA(pyw))[0];
	imuMeas.w_(1) = ((double*)PyArray_DATA(pyw))[1];
	imuMeas.w_(2) = ((double*)PyArray_DATA(pyw))[2];
	pManager_->addImuMeas(t,imuMeas);
}

int PythonManager::getImuMeas_python(double t, PyObject* pyf, PyObject* pyw){
	const ImuMeas* imuMeas;
	imuMeas = pManager_->getImuMeas(t);
	if(imuMeas != NULL){
		((double*)PyArray_DATA(pyf))[0] = imuMeas->f_(0);
		((double*)PyArray_DATA(pyf))[1] = imuMeas->f_(1);
		((double*)PyArray_DATA(pyf))[2] = imuMeas->f_(2);
		((double*)PyArray_DATA(pyw))[0] = imuMeas->w_(0);
		((double*)PyArray_DATA(pyw))[1] = imuMeas->w_(1);
		((double*)PyArray_DATA(pyw))[2] = imuMeas->w_(2);
		return 1;
	} else {
		return 0;
	}
}

void PythonManager::addEncMeas_python(double t, PyObject* pye, PyObject* pyv, PyObject* pyCF){
	EncMeas encMeas;
	for(int i=0;i<LSE_N_LEG*LSE_DOF_LEG;i++){
		encMeas.e_(i) = ((double*)PyArray_DATA(pye))[i];
		encMeas.v_(i) = ((double*)PyArray_DATA(pyv))[i];
	}
	for(int i=0;i<LSE_N_LEG;i++){
		encMeas.CF_[i] = (int)(0<((double*)PyArray_DATA(pyCF))[i]);
	}
	pManager_->addEncMeas(t,encMeas);
}

int PythonManager::getEncMeas_python(double t, PyObject* pye, PyObject* pyv, PyObject* pyCF){
	const EncMeas* encMeas;
	encMeas = pManager_->getEncMeas(t);
	if(encMeas != NULL){
		for(int i=0;i<LSE_N_LEG*LSE_DOF_LEG;i++){
			((double*)PyArray_DATA(pye))[i] = encMeas->e_(i);
			((double*)PyArray_DATA(pyv))[i] = encMeas->v_(i);
		}
		for(int i=0;i<LSE_N_LEG;i++){
			((double*)PyArray_DATA(pyCF))[i] = (double)encMeas->CF_[i];
		}
		return 1;
	} else {
		return 0;
	}
}

void PythonManager::addPosMeas_python(double t, PyObject* pyr, PyObject* pyq){
	PosMeas posMeas;
	posMeas.r_(0) = ((double*)PyArray_DATA(pyr))[0];
	posMeas.r_(1) = ((double*)PyArray_DATA(pyr))[1];
	posMeas.r_(2) = ((double*)PyArray_DATA(pyr))[2];
	posMeas.q_(0) = ((double*)PyArray_DATA(pyq))[0];
	posMeas.q_(1) = ((double*)PyArray_DATA(pyq))[1];
	posMeas.q_(2) = ((double*)PyArray_DATA(pyq))[2];
	posMeas.q_(3) = ((double*)PyArray_DATA(pyq))[3];
	pManager_->addPosMeas(t,posMeas);
}

int PythonManager::getPosMeas_python(double t, PyObject* pyr, PyObject* pyq){
	const PosMeas* posMeas;
	posMeas = pManager_->getPosMeas(t);
	if(posMeas != NULL){
		((double*)PyArray_DATA(pyr))[0] = posMeas->r_(0);
		((double*)PyArray_DATA(pyr))[1] = posMeas->r_(1);
		((double*)PyArray_DATA(pyr))[2] = posMeas->r_(2);
		((double*)PyArray_DATA(pyq))[0] = posMeas->q_(0);
		((double*)PyArray_DATA(pyq))[1] = posMeas->q_(1);
		((double*)PyArray_DATA(pyq))[2] = posMeas->q_(2);
		((double*)PyArray_DATA(pyq))[3] = posMeas->q_(3);
		return 1;
	} else {
		return 0;
	}
}

void PythonManager::addOflMeas_python(double t, PyObject* pyx, PyObject* pyu){
	OflMeas oflMeas;
	int N = PyArray_SIZE(pyx);
	N = floor(((double)N)/3.0);
	oflMeas.x_.reserve(N);
	oflMeas.u_.reserve(N);
	Eigen::Vector3d v;
	v.setZero();
	for(int i=0;i<N;i++){
		v(0) = ((double*)PyArray_DATA(pyx))[i*3+0];
		v(1) = ((double*)PyArray_DATA(pyx))[i*3+1];
		v(2) = ((double*)PyArray_DATA(pyx))[i*3+2];
		oflMeas.x_.push_back(v);
		v(0) = ((double*)PyArray_DATA(pyu))[i*3+0];
		v(1) = ((double*)PyArray_DATA(pyu))[i*3+1];
		v(2) = ((double*)PyArray_DATA(pyu))[i*3+2];
		oflMeas.u_.push_back(v);
	}
	pManager_->addOflMeas(t,oflMeas);
}

void PythonManager::clearMeas_python(){
	pManager_->clearMeas();
}

void PythonManager::update_pythont(double t){
	return pManager_->update(t);
}
void PythonManager::update_python(){
	return pManager_->update();
}
void PythonManager::getEst_python(PyObject* pyx){
	State x;
	x = pManager_->getEst();
	((double*)PyArray_DATA(pyx))[0] = x.r_(0);
	((double*)PyArray_DATA(pyx))[1] = x.r_(1);
	((double*)PyArray_DATA(pyx))[2] = x.r_(2);
	((double*)PyArray_DATA(pyx))[3] = x.v_(0);
	((double*)PyArray_DATA(pyx))[4] = x.v_(1);
	((double*)PyArray_DATA(pyx))[5] = x.v_(2);
	((double*)PyArray_DATA(pyx))[6] = x.q_(0);
	((double*)PyArray_DATA(pyx))[7] = x.q_(1);
	((double*)PyArray_DATA(pyx))[8] = x.q_(2);
	((double*)PyArray_DATA(pyx))[9] = x.q_(3);
	((double*)PyArray_DATA(pyx))[10] = x.w_(0);
	((double*)PyArray_DATA(pyx))[11] = x.w_(1);
	((double*)PyArray_DATA(pyx))[12] = x.w_(2);
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
#if USE_CERES
int PythonManager::robotCalibration_python(double t, double T){
	return pManager_->robotCalibration(t,T);
}
int PythonManager::getLengthOfBC_python(){
	return pManager_->getLengthOfBC();
}
void PythonManager::getBCData_python(PyObject* X){
	int N = pManager_->getLengthOfBC();
	const RobotCalibration::state* state = pManager_->getBCData();
	for(int i=0;i<N;i++){
		((double*)PyArray_DATA(X))[i*17] = state[i].t_;
		((double*)PyArray_DATA(X))[i*17+1] = state[i].r_[0];
		((double*)PyArray_DATA(X))[i*17+2] = state[i].r_[1];
		((double*)PyArray_DATA(X))[i*17+3] = state[i].r_[2];
		((double*)PyArray_DATA(X))[i*17+4] = state[i].v_[0];
		((double*)PyArray_DATA(X))[i*17+5] = state[i].v_[1];
		((double*)PyArray_DATA(X))[i*17+6] = state[i].v_[2];
		((double*)PyArray_DATA(X))[i*17+7] = state[i].a_[0];
		((double*)PyArray_DATA(X))[i*17+8] = state[i].a_[1];
		((double*)PyArray_DATA(X))[i*17+9] = state[i].a_[2];
		((double*)PyArray_DATA(X))[i*17+10] = state[i].q_[0];
		((double*)PyArray_DATA(X))[i*17+11] = state[i].q_[1];
		((double*)PyArray_DATA(X))[i*17+12] = state[i].q_[2];
		((double*)PyArray_DATA(X))[i*17+13] = state[i].q_[3];
		((double*)PyArray_DATA(X))[i*17+14] = state[i].w_[0];
		((double*)PyArray_DATA(X))[i*17+15] = state[i].w_[1];
		((double*)PyArray_DATA(X))[i*17+16] = state[i].w_[2];
	}

}
#endif
using namespace boost::python;
BOOST_PYTHON_MODULE(_PythonManager)
{
    class_<PythonManager>("PythonManager", init<std::string>())
        .def("addImuMeas", &PythonManager::addImuMeas_python)
        .def("getImuMeas", &PythonManager::getImuMeas_python)
        .def("addEncMeas", &PythonManager::addEncMeas_python)
        .def("getEncMeas", &PythonManager::getEncMeas_python)
        .def("addPosMeas", &PythonManager::addPosMeas_python)
        .def("getPosMeas", &PythonManager::getPosMeas_python)
        .def("clearMeas", &PythonManager::clearMeas_python)
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
#if USE_CERES
        .def("robotCalibration", &PythonManager::robotCalibration_python)
        .def("getLengthOfBC", &PythonManager::getLengthOfBC_python)
        .def("getBCData", &PythonManager::getBCData_python)
#endif
    ;
}
#endif

}
