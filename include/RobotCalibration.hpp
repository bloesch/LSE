/*!
* @file 	RobotCalibration.hpp
* @author 	Michael Blösch
* @date		10.10.2012
* @brief	Time delay calibration routine
 */

#ifndef RobotCalibration_HPP_
#define RobotCalibration_HPP_

#include "Common.hpp"
#include <Eigen/Dense>
#include "Rotations.hpp"
#include <map>
#include <vector>

namespace LSE {

class Manager;

/*! Time delay calibration routine */
class RobotCalibration{
public:
	struct state{
		state(){
			t_ = 0;
			r_[0] = 0;
			r_[1] = 0;
			r_[2] = 0;
			v_[0] = 0;
			v_[1] = 0;
			v_[2] = 0;
			a_[0] = 0;
			a_[1] = 0;
			a_[2] = 0;
			q_[0] = 1;
			q_[1] = 0;
			q_[2] = 0;
			q_[3] = 0;
			w_[0] = 0;
			w_[1] = 0;
			w_[2] = 0;
		}
		double t_;
		double r_[3];
		double v_[3];
		double a_[3];
		double q_[4];
		double w_[3];
	};
	state* PB_states_;

	/* -------------------- Constructor/Destructor and calibration call --------------------- */
	/*! Constructor
	 * @param[in]	pManager	pointer to main class Manager
	 * @param[in]	pFilename	filename of parameter-file
	 */
	RobotCalibration(Manager* pManager,const char* pFilename);
	/*! Destructor */
	~RobotCalibration();
	/*! Identifies the delay between the sensor modalities
	 * @return	results of identification (0:fail, 1:success)
	 * @param[in]	t	end of identification interval
	 * @param[in]	T	length of identification interval
	 */
	int calibrateRobot(const double& t,const double& T);

	int getN();
	const state* getBatch();


private:
	/* -------------------- Various private functions --------------------- */
	/*! Checks and initializes the calibration
	 * @return	results of check and initialization (0:fail, 1:success)
	 * @param[in]	t	end of identification interval
	 * @param[in]	T	length of identification interval
	 */
	int initialize(const double& t,const double& T);
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

	/* -------------------- Vectors and map iterators --------------------- */
	/*! Iterator for the IMU measurements */
	std::map<double,ImuMeas>::iterator itImu_;
	/*! Iterator for the IMU measurements (first) */
	std::map<double,ImuMeas>::iterator itImu1_;
	/*! Iterator for the IMU measurements (second) */
	std::map<double,ImuMeas>::iterator itImu2_;
	/*! Iterator for the kinematic measurements (first) */
	std::map<double,EncMeas>::iterator itEnc_;
//	/*! Iterator for the kinematic measurements (second) */
//	std::map<double,EncMeas>::iterator itEnc2_;
//	/*! Reverse iterator for the kinematic measurements (first) */
//	std::map<double,EncMeas>::reverse_iterator ritEnc_;
//	/*! Reverse iterator for the kinematic measurements (second) */
//	std::map<double,EncMeas>::reverse_iterator ritEnc2_;
	/*! Iterator for the pose sensor measurements (first) */
	std::map<double,PosMeas>::iterator itPos_;
//	/*! Iterator for the pose sensor measurements (second) */
//	std::map<double,PosMeas>::iterator itPos2_;
//	/*! Reverse iterator for the pose sensor measurements (first) */
//	std::map<double,PosMeas>::reverse_iterator ritPos_;
//	/*! Reverse iterator for the pose sensor measurements (second) */
//	std::map<double,PosMeas>::reverse_iterator ritPos2_;

	/* -------------------- Temporary quantities used during differentiation and interpolation --------------------- */
//	/*! Time of last evaluation */
//	double lastTime_;
//	/*! Time of new evaluation */
//	double newTime_;
//	/*! Quaternion describing earlier attitude */
//	Rotations::Quat q1_;
//	/*! Quaternion describing newer attitude */
//	Rotations::Quat q2_;
//	/*! Rotational rate norm of last time */
//	double lastNorm_;
//	/*! Rotational rate norm of new time */
//	double newNorm_;
//	/*! Current desired interpolation time */
//	double interpolTime_;

	/* -------------------- Calibration routine parameters --------------------- */
	/*! Timesteps between interpolation points, if 0 use IMU timesteps */
	double dt_;
	/*! Flag whether to use the IMU measurements */
	bool mbUseImu_;
	/*! Flag whether to use the kinematic measurements */
	bool mbUseEnc_;
	/*! Flag whether to use the pose sensor measurements */
	bool mbUsePos_;

	Eigen::Matrix3d Wr_;
	Eigen::Matrix3d Wv_;
	Eigen::Matrix3d Wq_;

	int N_;
	double PB_T_BI_[7];
	double PB_p_[12];
	double PB_IrIB_[3];
	double PB_qBI_[4];
	double PB_bf_[3];
	double PB_bw_[3];
	double PB_IrIC_[3];
	double PB_qCI_[4];
	double PB_WrWV_[3];
	double PB_qVW_[4];
	double PB_pkin_[17];
};

}

#endif /* RobotCalibration_HPP_ */
