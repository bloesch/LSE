/*!
* @file 	DelayCalibration.hpp
* @author 	Michael Blösch
* @date		10.10.2012
* @brief	Time delay calibration routine
 */

#ifndef DELAYCALIBRATION_HPP_
#define DELAYCALIBRATION_HPP_

#include "Common.hpp"
#include <Eigen/Dense>
#include "Rotations.hpp"
#include <map>
#include <vector>

namespace LSE {

class Manager;

/*! Time delay calibration routine */
class DelayCalibration{
public:
	/* -------------------- Constructor/Destructor and calibration call --------------------- */
	/*! Constructor
	 * @param[in]	pManager	pointer to main class Manager
	 * @param[in]	pFilename	filename of parameter-file
	 */
	DelayCalibration(Manager* pManager,const char* pFilename);
	/*! Destructor */
	~DelayCalibration();
	/*! Identifies the delay between the sensor modalities
	 * @return	results of identification (0:fail, 1:success)
	 * @param[in]	t	end of identification interval
	 * @param[in]	T	length of identification interval
	 */
	int calibrateDelay(const double& t,const double& T);


private:
	/* -------------------- Various private functions --------------------- */
	/*! Checks and initializes the calibration
	 * @return	results of check and initialization (0:fail, 1:success)
	 * @param[in]	t	end of identification interval
	 * @param[in]	T	length of identification interval
	 */
	int initialize(const double& t,const double& T);
	/*! Evaluates the rotational rate norms of the IMU using the current settings */
	void getImuNorms();
	/*! Evaluates the estimated rotational rate norms of the encoder data using the current settings */
	void getEncNorms();
	/*! Evaluates the estimated rotational rate norms of the pose sensor data using the current settings */
	void getPosNorms();
	/*! Computes the attitude of the main body given some fixed footholds and kinematic measurements
	 * @return	attitude of main body (quaternion)
	 * @param[in]	m	kinematic measurement
	 * @param[in]	p	location of footholds
	 */
	Rotations::Quat quatFromFootholds(const EncMeas& m,const Eigen::Matrix<double,3,LSE_N_LEG>& p);
	/*! Loads overall parameters from parameter file
	 * @param[in]	pFilename	name of parameter file
	 */
	void loadParam(const char* pFilename);

	/* -------------------- Pointers and timing stuff --------------------- */
	/*! Pointer to main class Manager */
	Manager* pManager_;
	/*! Start time of calibration interval */
	double t1_;
	/*! End time of calibration interval */
	double t2_;
	/*! Number of interpolation points */
	int N_;
	/*! Sampling time of IMU */
	double TsImu_;

	/* -------------------- Vectors and map iterators --------------------- */
	/*! Storage of IMU rotational rate norms */
	std::vector<double> imuRateNorm_;
	/*! Storage of kinematic rotational rate norms */
	std::vector<double> encRateNorm_;
	/*! Storage of pose sensor rotational rate norms */
	std::vector<double> posRateNorm_;
	/*! Iterator for the IMU measurements */
	std::map<double,ImuMeas>::iterator itImu_;
	/*! Iterator for the kinematic measurements (first) */
	std::map<double,EncMeas>::iterator itEnc_;
	/*! Iterator for the kinematic measurements (second) */
	std::map<double,EncMeas>::iterator itEnc2_;
	/*! Reverse iterator for the kinematic measurements (first) */
	std::map<double,EncMeas>::reverse_iterator ritEnc_;
	/*! Reverse iterator for the kinematic measurements (second) */
	std::map<double,EncMeas>::reverse_iterator ritEnc2_;
	/*! Iterator for the pose sensor measurements (first) */
	std::map<double,PosMeas>::iterator itPos_;
	/*! Iterator for the pose sensor measurements (second) */
	std::map<double,PosMeas>::iterator itPos2_;
	/*! Reverse iterator for the pose sensor measurements (first) */
	std::map<double,PosMeas>::reverse_iterator ritPos_;
	/*! Reverse iterator for the pose sensor measurements (second) */
	std::map<double,PosMeas>::reverse_iterator ritPos2_;

	/* -------------------- Temporary quantities used during differentiation and interpolation --------------------- */
	/*! Time of last evaluation */
	double lastTime_;
	/*! Time of new evaluation */
	double newTime_;
	/*! Quaternion describing earlier attitude */
	Rotations::Quat q1_;
	/*! Quaternion describing newer attitude */
	Rotations::Quat q2_;
	/*! Rotational rate norm of last time */
	double lastNorm_;
	/*! Rotational rate norm of new time */
	double newNorm_;
	/*! Current desired interpolation time */
	double interpolTime_;

	/* -------------------- Calibration routine parameters --------------------- */
	/*! Timesteps between interpolation points */
	double dt_;
	/*! Upper bound on the delay between the sensor modalities */
	double maxDelay_;
	/*! Discrete differentiation window for kinematic measurements */
	int difWindowEnc_;
	/*! Discrete differentiation window for pose sensor measurements */
	int difWindowPos_;
	/*! Flag whether to use the IMU measurements */
	bool mbUseImu_;
	/*! Flag whether to use the kinematic measurements */
	bool mbUseEnc_;
	/*! Flag whether to use the pose sensor measurements */
	bool mbUsePos_;
};

}

#endif /* DELAYCALIBRATION_HPP_ */
