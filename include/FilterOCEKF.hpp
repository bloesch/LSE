/*!
* @file 	FilterOCEKF.hpp
* @author 	Michael Blösch
* @date		10.10.2012
* @brief	OCEKF Filter for legged robots
 */

#ifndef FILTEROCEKF_HPP_
#define FILTEROCEKF_HPP_

#include "FilterBase.hpp"
#include "Common.hpp"
#include <Eigen/Dense>
#include "Rotations.hpp"

namespace LSE {

class Manager;

/*! Observability Constrained Extended Kalman Filter */
class FilterOCEKF: public FilterBase{
public:
	/* -------------------- Constructor/Destructor --------------------- */
	/*! Constructor
	 * @param[in]	pManager	pointer to main class Manager
	 * @param[in]	pFilename	filename of parameter-file
	 */
	FilterOCEKF(Manager* pManager,const char* pFilename);
	/*! Destructor */
	virtual ~FilterOCEKF();

	/* -------------------- Filter handling --------------------- */
	/*! Updates the filter to time t
	 * @param[in]	t	desired update time
	 */
	virtual void update(const double& t);
	/*! Updates the filter to the newest measurement time */
	virtual void update();
	/*! Return current estimate of robot state (main body)
	 * @return	current robot state
	 */
	virtual State getEst();
	/*! Return slippage detection
	 * @return	current slippage detection
	 */
	virtual SlippageDetection getSlippage();
	/*! Resets the filter
	 * @param[in]	t	time used to initialize new state estimate
	 */
	virtual void resetEstimate(const double& t);
	/*! Returns a string describing the main filter parameters
	 * @param[out] str	string characterize the parameter set of the filter
	 */
	virtual std::string getKeyString();
	/*! Makes and entry of the cuurent state into the log-file */
	void logState();


private:
	/*! Loads overall parameters from parameter file
	 * @param[in]	pFilename	name of parameter file
	 */
	void loadParam(const char* pFilename);
	/*! Structure of filter intern state */
	struct InternState{
		EIGEN_MAKE_ALIGNED_OPERATOR_NEW

		/*! Time of estimate */
		double t_;
		/*! Position estimate */
		Eigen::Vector3d r_;
		/*! Velocity estimate */
		Eigen::Vector3d v_;
		/*! Attitude estimate (quaternion) */
		Rotations::Quat q_;
		/*! Rotational rate estimate (bias corrected) */
		Eigen::Vector3d w_;
		/*! Foothold estimate */
		Eigen::Matrix<double,3,LSE_N_LEG> p_;
		/*! Contact flag counter */
		CF CFC_;
		/*! Estimate of accelerometer bias */
		Eigen::Vector3d bf_;
		/*! Estimate of gyroscope bias */
		Eigen::Vector3d bw_;
		/*! Estimate of covariance matrix */
		Eigen::Matrix<double,15+3*LSE_N_LEG,15+3*LSE_N_LEG> P_;
		/*! Last position estimate */
		Eigen::Vector3d rLast_;
		/*! Last velocity estimate */
		Eigen::Vector3d vLast_;
		/*! Last attitude estimate */
		Rotations::Quat qLast_;
		/*! Linearization point for footholds */
		Eigen::Matrix<double,3,LSE_N_LEG> pLin_;
		/*! Linearization point for rotational rate */
		Eigen::Vector3d wLin_;
		/*! Linearization point for accelerometer measurement for position Jacobian */
		Eigen::Vector3d f1Lin_;
		/*! Linearization point for accelerometer measurement for velocity Jacobian */
		Eigen::Vector3d f2Lin_;
		/*! Current corrected accelerometer measurement */
		Eigen::Vector3d f_;
		/*! Time of last update */
		double tLast_;
	};

	/* -------------------- Filtering/Predicting/Updating --------------------- */
	/*! Filters the referenced internal state up to the given time
	 * @param[in/out]	x	Filter state to be filtered
	 * @param[in]		t	Desired filter time
	 */
	void filterState(InternState& x,const double& t);
	/*! Predicts the referenced internal state using the given IMU measurement
	 * @param[in/out]	x	Filter state to be filtered
	 * @param[in]		t	Desired filter time
	 * @param[in]		m	IMU measurement
	 */
	void predictState(InternState& x,const double& t, const ImuMeas& m);
	/*! Updates the referenced internal state using the given Encoder measurement
	 * @param[in/out]	x	Filter state to be filtered
	 * @param[in]		m	Encoder measurement
	 */
	void encUpdateState(InternState& x,const EncMeas& m);

	/* -------------------- Pointers and filter states --------------------- */
	/*! Pointer to main class Manager */
	Manager* pManager_;
	/*! Safe state (where the chance is high that all measurements have arrived) */
	InternState xs_;
	/*! Predicted state */
	InternState xp_;

	/* -------------------- Parameters (can be specified in the parameter file) --------------------- */
	/*! Initialization state */
	InternState xInit_;
	/*! Predicition noise of position */
	Eigen::Matrix3d Wr_;
	/*! Predicition noise of velocity */
	Eigen::Matrix3d Wv_;
	/*! Predicition noise of attitude */
	Eigen::Matrix3d Wq_;
	/*! Predicition noise of footholds */
	Eigen::Matrix3d Wp_;
	/*! Predicition noise of accelerometer bias */
	Eigen::Matrix3d Wbf_;
	/*! Predicition noise of gyroscope bias */
	Eigen::Matrix3d Wbw_;

	/* -------------------- Flag for modes (can be specified in the parameter file) --------------------- */
	/*! True if accelerometer bias is co-estimated */
	bool mbEstimateAccBias_;
	/*! True if gyroscope bias is co-estimated */
	bool mbEstimateRotBias_;
	/*! True if imu measurements are used */
	bool mbUseImu_;
	/*! True if kinematic measurements are used */
	bool mbUseKin_;
	/*! True if the floor should be assumed to be flat */
	bool mbAssumeFlatFloor_;
	/*! True if using fixed timestepping (bug handling in SL)*/
	bool mbFixedTimeStepping_;
	/*! Timestep used if fixed time stepping*/
	double timeStep_;
};

}

#endif /* FILTEROCEKF_HPP_ */
