/*!
* @file 	Manager.hpp
* @author 	Michael Blösch
* @date		10.10.2012
* @brief	Main class of legged state estimator
* 			Handles measurements and filters
 */

#ifndef LSE_MANAGER_HPP_
#define LSE_MANAGER_HPP_

#include "Common.hpp"
#include "Rotations.hpp"
#include <Eigen/Dense>
#include <map>

namespace LSE {

class FilterOCEKF;
class FilterSync;
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
	~Manager();

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
	/*! Sets the Sampling time (only for fixed interval filters (SyncFilter)
	 * @param[in]	Ts	sampling time
	 */
	void setSamplingTime(double Ts);

	/* -------------------- Friends --------------------- */
	friend class FilterOCEKF;
	friend class FilterSync;
	friend class DelayCalibration;

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
	/*! Factorial
	 * @return	factorial of k
	 * @param[in]	k	Natural number
	 */
	int factorial(const int& k);

	/* -------------------- Different pointers --------------------- */
	/* Index of active Filter */
	int activeFilter_;
	/*! Pointer to OCEKF Filter */
	FilterOCEKF* pFilterOCEKF_;
	/*! Pointer to Sync Filter */
	FilterSync* pFilterSync_;
	/*! Pointer to time delay calibration routine */
	DelayCalibration* pDelayCalibration_;
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
	/*! Position of pose sensor frame w.r.t the body frame (expressed in body frame) */
	Eigen::Vector3d B_r_BK_;
	/*! Rotation from body frame to pose sensor frame */
	Rotations::Quat q_KB_;
	/*! Noise of accelerometer [m^2/s^3] (continuous form) */
	Eigen::Matrix3d Rf_;
	/*! Noise of gyroscope [rad^2/s] (continuous form) */
	Eigen::Matrix3d Rw_;
	/*! Noise of foothold measurements [m^2] (discrete form) */
	Eigen::Matrix3d Rs_;
	/*! Noise of encoder measurement [rad^2] (discrete form) */
	Eigen::Matrix<double,LSE_DOF_LEG,LSE_DOF_LEG> Ra_;

};

}

#endif /* LSE_MANAGER_HPP_ */
