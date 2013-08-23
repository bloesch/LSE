/*!
* @file 	Manager.hpp
* @author 	Michael Blösch
* @date		10.10.2012
* @brief	Main class of legged state estimator
* 			Handles measurements and filters
 */

#ifndef LSE_MANAGER_HPP_
#define LSE_MANAGER_HPP_

#ifdef USE_CERES
#define NUM_FILTERS 4
#else
#define NUM_FILTERS 3
#endif

#include "FilterBase.hpp"
#include "Common.hpp"
#include "Rotations.hpp"
#include "OptimizationFramework.hpp"
#ifdef USE_CERES
#include "RobotCalibration.hpp"
#endif
#include <Eigen/Dense>
#include <map>
#include <iostream>
#include <sstream>
#include <fstream>

namespace LSE {

class DelayCalibration;

/*! Manager */
class Manager{
public:
	/* -------------------- Constructor/Destructor --------------------- */
	/*! Constructor
	 * @param[in]	pFilename	filename of parameter-file
	 * @param[in] 	f			function pointer to forward kinematic
	 * @param[in] 	J			function pointer to forward kinematic Jacobian (default = NULL)
	 */
	Manager(const char* pFilename,Eigen::Vector3d (*f)(Eigen::Matrix<double,LSE_DOF_LEG,1>,int),Eigen::Matrix<double,3,LSE_DOF_LEG> (*J)(Eigen::Matrix<double,LSE_DOF_LEG,1>,int) = NULL);
	/*! Destructor */
	virtual ~Manager();

	/* -------------------- Measurement Handling --------------------- */
	/*! Adds an IMU measurement
	 * @param[in]	t	time of measurement
	 * @param[in] 	m	measurement data
	 */
	void addImuMeas(const double& t,const ImuMeas& m);
	/*! Adds encoders measurements
	 * @param[in]	t	time of measurement
	 * @param[in] 	m	measurement data
	 */
	void addEncMeas(const double& t,const EncMeas& m);
	/*! Adds pose sensor measurements
	 * @param[in]	t	time of measurement
	 * @param[in] 	m	measurement data
	 */
	void addPosMeas(const double& t,const PosMeas& m);
	/*! Searches for the first measurement after time t
	 * @return const pointer to measurement
	 * @param[in/out]	t	time of measurement, changed to precise measurement time
	 */
	const ImuMeas* getImuMeas(double& t);
	/*! Searches for the first measurement after time t
	 * @return const pointer to measurement
	 * @param[in/out]	t	time of measurement, changed to precise measurement time
	 */
	const EncMeas* getEncMeas(double& t);
	/*! Searches for the first measurement after time t
	 * @return const pointer to measurement
	 * @param[in/out]	t	time of measurement, changed to precise measurement time
	 */
	const PosMeas* getPosMeas(double& t);
	/*! Clears all measurement entries
	 */
	void clearMeas();

	void addOflMeas(const double& t,const OflMeas& m);
	const OflMeas* getOflMeas(double& t);

	/* -------------------- Filter and Calibration Handling --------------------- */
	/*! Updates the filter to time t
	 * @param[in]	t	desired update time
	 */
	void update(const double& t);
	/*! Updates the filter to the newest measurement time */
	void update();
	/*! Return current estimate of robot state (main body)
	 * @return	current robot state
	 */
	State getEst();
	/*! Resets the filter
	 * @param[in]	t	time used to initialize new state estimate
	 */
	void resetEstimate(const double& t);
	/*! Identifies the delay between the sensor modalities
	 * @return	results of identification (0:fail, 1:success)
	 * @param[in]	t	end of identification interval
	 * @param[in]	T	length of identification interval
	 */
	int delayIdentification(const double& t,const double& T);
	/*! Calibrates the robot
	 * @return	results of identification (0:fail, 1:success)
	 * @param[in]	t	end of identification interval
	 * @param[in]	T	length of identification interval
	 */
#ifdef USE_CERES
	int robotCalibration(const double& t,const double& T);
#endif

	/* -------------------- Time delay handling of modalities --------------------- */
	/*! Set the time delay parameter of the IMU
	 * @param[in]	TD	time delay
	 */
	void setImuTD(const double& TD);
	/*! Set the time delay parameter of the encoders
	 * @param[in]	TD	time delay
	 */
	void setEncTD(const double& TD);
	/*! Set the time delay parameter of the pose sensor
	 * @param[in]	TD	time delay
	 */
	void setPosTD(const double& TD);
	/*! Get the time delay parameter of the IMU
	 * @return	time delay
	 */
	double getImuTD();
	/*! Get the time delay parameter of the IMU
	 * @return	time delay
	 */
	double getEncTD();
	/*! Get the time delay parameter of the IMU
	 * @return	time delay
	 */
	double getPosTD();

#ifdef USE_CERES
	int getLengthOfBC();
	const RobotCalibration::state* getBCData();
#endif


	/* -------------------- Logging stuff (unclean) --------------------- */
	void enableLogging(const char* pLogfile);
	void disableLogging();

	/* -------------------- Friends --------------------- */
	friend class FilterOCEKF;
	friend class FilterVUKF;
	friend class FilterVUKF2;
	friend class DelayCalibration;
#ifdef USE_CERES
	friend class FilterInertialOF;
	friend class FilterFLS;
	friend class RobotCalibration;
#endif

private:
	/*! Loads overall parameters from parameter file
	 * @param[in]	pFilename	name of parameter file
	 */
	void loadParam(const char* pFilename);
	/*! Returns gamma 3x3 matrix
	 * @return	Gamma matrix
	 * @param[in]	k	Order
	 * @param[in]	w	Rotational rate
	 * @param[in]	dt	Time difference
	 */
	Eigen::Matrix3d gamma(const int& k,const Eigen::Vector3d& w,const double& dt);
	/*! Returns gamma 3x3 matrix
	 * @return	Gamma matrix
	 * @param[in]	k	Order
	 * @param[in]	v	Vector
	 */
	Eigen::Matrix3d gamma(const int& k,const Eigen::Vector3d& v);
	/*! Factorial
	 * @return	factorial of k
	 * @param[in]	k	Natural number
	 */
	int factorial(const int& k);

	/* -------------------- Different pointers --------------------- */
	/*! Pointer to Filter list */
	FilterBase* pFilterList_[NUM_FILTERS];
	/*! Index of active filter */
	int activeFilter_;
	/*! Pointer to time delay calibration routine */
	DelayCalibration* pDelayCalibration_;
	/*! Pointer to time robot calibration routine */
#ifdef USE_CERES
	RobotCalibration* pRobotCalibration_;
#endif
	/*! Function pointer to leg kinematics */
	Eigen::Vector3d (*legKin)(Eigen::Matrix<double,LSE_DOF_LEG,1>,int);
	/*! Function pointer to leg kinematics Jacobian */
	Eigen::Matrix<double,3,LSE_DOF_LEG> (*legKinJac)(Eigen::Matrix<double,LSE_DOF_LEG,1>,int);

	/* -------------------- Measurement Storage --------------------- */
	/*! Map storage of Imu Measurements */
	std::map<double,ImuMeas> imuMeasList_;
	/*! Map storage of encoder Measurements */
	std::map<double,EncMeas> encMeasList_;
	/*! Map storage of pose sensor Measurements */
	std::map<double,PosMeas> posMeasList_;
	/*! Map storage of pose sensor Measurements */
	std::map<double,OflMeas> oflMeasList_;

	/* -------------------- Parameters --------------------- */
	/*! Gravity vector in world coordinate frame */
	const Eigen::Vector3d g_;
	/*! Imu time offset, real time = timestamp + tImu_ */
	double tImu_;
	/*! Encoder time offset, real time = timestamp + tEnc_ */
	double tEnc_;
	/*! Pose sensor time offset, real time = timestamp + tPos_ */
	double tPos_;
	/*! Position of Imu frame w.r.t the body frame (expressed in body frame) */
	Eigen::Vector3d B_r_BI_;
	/*! Rotation from body frame to Imu frame */
	Rotations::Quat q_IB_;
	/*! Position of kinematic frame w.r.t the body frame (expressed in body frame) */
	Eigen::Vector3d B_r_BK_;
	/*! Rotation from body frame to kinematic frame */
	Rotations::Quat q_KB_;
	/*! Noise of accelerometer [m^2/s^3] (continuous form) */
	Eigen::Matrix3d Rf_;
	/*! Noise of gyroscope [rad^2/s] (continuous form) */
	Eigen::Matrix3d Rw_;
	/*! Noise of foothold measurements [m^2] (discrete form) */
	Eigen::Matrix3d Rs_;
	/*! Noise of pose measurements [m^2] (discrete form) */
	Eigen::Matrix3d Rposr_;
	/*! Noise of pose measurements [rad^2] (discrete form) */
	Eigen::Matrix3d Rposq_;
	/*! Noise of encoder measurement [rad^2] (position) (discrete form) */
	Eigen::Matrix<double,LSE_DOF_LEG,LSE_DOF_LEG> Ra_;
	/*! Noise of encoder measurement [rad^2/s^2] (velocity) (discrete form) */
	Eigen::Matrix<double,LSE_DOF_LEG,LSE_DOF_LEG> Rda_;


	/* -------------------- Logging stuff (unclean) --------------------- */
	std::ofstream ofsLog_;
	bool isLogging_;


};

}

#endif /* LSE_MANAGER_HPP_ */
