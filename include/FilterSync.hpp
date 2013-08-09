/*!
* @file 	FilterSync.hpp
* @author 	Michael Blösch
* @date		10.10.2012
* @brief	Sync Filter for legged robots
 */

#ifndef FILTERSYNC_HPP_
#define FILTERSYNC_HPP_

#define SF_state_dim (15+3*LSE_N_LEG)
#define SF_preNoise_dim (21+3*LSE_N_LEG)
#define SF_upNoise_dim (3*LSE_N_LEG)
#define g 9.81

#include "Common.hpp"
#include <Eigen/Dense>
#include "Rotations.hpp"

namespace LSE {

class Manager;

/*! Minimal Filter State */
class SyncFilterState{
public:
	/*! Position estimate */
	Eigen::Vector3d r_;
	/*! Velocity estimate */
	Eigen::Vector3d v_;
	/*! Attitude estimate (quaternion) */
	Rotations::Quat q_;
	/*! Foothold estimate */
	Eigen::Matrix<double,3,LSE_N_LEG> p_;
	/*! Estimate of accelerometer bias */
	Eigen::Vector3d bf_;
	/*! Estimate of gyroscope bias */
	Eigen::Vector3d bw_;

	SyncFilterState operator +(const Eigen::Matrix<double,SF_state_dim,1> &y) const;
	Eigen::Matrix<double,SF_state_dim,1> operator -(const SyncFilterState &y) const;
};

/*! Observability Constrained Extended Kalman Filter */
class FilterSync{
public:
	/* -------------------- Constructor/Destructor --------------------- */
	/*! Constructor
	 * @param[in]	pManager	pointer to main class Manager
	 * @param[in]	pFilename	filename of parameter-file
	 */
	FilterSync(Manager* pManager,const char* pFilename);
	/*! Destructor */
	~FilterSync();

	/* -------------------- Filter handling --------------------- */
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


private:
	typedef Eigen::Matrix<double,SF_state_dim,SF_state_dim> MatrixP;
	typedef Eigen::Matrix<double,SF_preNoise_dim,SF_preNoise_dim> MatrixPreCov;
	typedef Eigen::Matrix<double,SF_upNoise_dim,SF_upNoise_dim> MatrixUpCov;
	/*! Loads overall parameters from parameter file
	 * @param[in]	pFilename	name of parameter file
	 */
	void loadParam(const char* pFilename);
	/*! Structure of filter intern state */
	struct InternState{
		EIGEN_MAKE_ALIGNED_OPERATOR_NEW

		/*! Time of estimate */
		double t_;
		/*! Minimal Filter State */
		SyncFilterState x_;
		/*! Estimate of covariance matrix */
		MatrixP P_;
		/*! Contact flag counter */
		CF CFC_;
		/*! Rotational rate estimate (bias corrected) */
		Eigen::Vector3d w_;
		/*! Current corrected accelerometer measurement */
		Eigen::Vector3d f_;
	};


	/* -------------------- Pointers and filter states --------------------- */
	void predict(SyncFilterState& x,double Ts,ImuMeas imuMeas);
	void predict(SyncFilterState& x,double Ts,ImuMeas imuMeas,Eigen::Matrix<double,SF_preNoise_dim,1> n);

	/* -------------------- Pointers and filter states --------------------- */
	/*! Pointer to main class Manager */
	Manager* pManager_;
	/*! State */
	InternState x_;
	/*! Sigma Samples of Filter States*/
	SyncFilterState X_[1+2*(SF_state_dim+SF_preNoise_dim+SF_upNoise_dim)];
	/*! Sigma Samples of Prediction Noise*/
	Eigen::Matrix<double,SF_preNoise_dim,1+2*SF_preNoise_dim> PN_;
	/*! Sigma Samples of Update Noise*/
	Eigen::Matrix<double,SF_upNoise_dim,1+2*SF_upNoise_dim> UN_;

	/*! Prediction noise matrix */
	MatrixPreCov Npre_;
	/*! Cholesky of Prediction noise matrix */
	MatrixPreCov SNpre_;
	/*! Update noise matrix */
	MatrixUpCov Nup_;
	/*! Cholesky of update noise matrix */
	MatrixUpCov SNup_;
	/*! Cholesky covariance matrix */
	MatrixP SP_;

	/* -------------------- Parameters of unscented filter --------------------- */
	/*! Alpha */
	double UKFAlpha_;
	/*! Kappa */
	double UKFKappa_;
	/*! Beta */
	double UKFBeta_;
	/*! Weights */
	double UKFWs_,UKFWc_,UKFWi_;
	/*! Lambdas */
	double UKFLambda_,UKFGamma_;

	/* -------------------- Parameters (can be specified in the parameter file) --------------------- */
	/*! Initialization state */
	InternState xInit_;
	/*! Prediction noise of position */
	Eigen::Matrix3d Wr_;
	/*! Prediction noise of velocity */
	Eigen::Matrix3d Wv_;
	/*! Prediction noise of attitude */
	Eigen::Matrix3d Wq_;
	/*! Prediction noise of footholds */
	Eigen::Matrix3d Wp_;
	/*! Prediction noise of accelerometer bias */
	Eigen::Matrix3d Wbf_;
	/*! Prediction noise of gyroscope bias */
	Eigen::Matrix3d Wbw_;
	/*! Fixed Time Step Value */
	double Ts_;

	/* -------------------- Flag for modes (can be specified in the parameter file) --------------------- */
	/*! True if accelerometer bias is co-estimated */
	bool mbEstimateAccBias_;
	/*! True if gyroscope bias is co-estimated */
	bool mbEstimateRotBias_;
};

}

#endif /* FILTERSYNC_HPP_ */
