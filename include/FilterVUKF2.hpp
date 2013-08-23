/*!
* @file 	FilterVUKF2.hpp
* @author 	Michael Blösch
* @date		10.10.2012
* @brief	VUKF Filter for legged robots
 */

#ifndef FilterVUKF2_HPP_
#define FilterVUKF2_HPP_

#define VUKFF_state_dim (15)
#define VUKFF_preNoise_dim (21)
#define VUKFF_upNoise_dim (3*LSE_N_LEG)

#include "FilterBase.hpp"
#include "Common.hpp"
#include <Eigen/Dense>
#include "Rotations.hpp"

namespace LSE {

class Manager;
class StateVUKF;
class PreNoiseVUKF;
class MeasKinNoiseVUKF;

/*! Minimal Filter State */
class VUKFFilterState{
public:
	/*! Position estimate */
	Eigen::Vector3d r_;
	/*! Velocity estimate */
	Eigen::Vector3d v_;
	/*! Attitude estimate (quaternion) */
	Rotations::Quat q_;
	/*! Estimate of accelerometer bias */
	Eigen::Vector3d bf_;
	/*! Estimate of gyroscope bias */
	Eigen::Vector3d bw_;

	VUKFFilterState operator +(const Eigen::Matrix<double,VUKFF_state_dim,1> &y) const;
	Eigen::Matrix<double,VUKFF_state_dim,1> operator -(const VUKFFilterState &y) const;
};

/*! Observability Constrained Extended Kalman Filter */
class FilterVUKF2: public FilterBase{
public:public:
	/* -------------------- Constructor/Destructor --------------------- */
	/*! Constructor
	 * @param[in]	pManager	pointer to main class Manager
	 * @param[in]	pFilename	filename of parameter-file
	 */
	FilterVUKF2(Manager* pManager,const char* pFilename);
	/*! Destructor */
	virtual ~FilterVUKF2();

	/* -------------------- Filter handling --------------------- */
	/*! Updates the filter to time t (may include prediction)
	 * @param[in]	t	desired update time
	 */
	virtual void update(const double& t);
	/*! Updates the filter to the newest measurement time (may include prediction)*/
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

private:
	typedef Eigen::Matrix<double,VUKFF_state_dim,VUKFF_state_dim> MatrixP;
	typedef Eigen::Matrix<double,VUKFF_preNoise_dim,VUKFF_preNoise_dim> MatrixPreCov;
	typedef Eigen::Matrix<double,VUKFF_upNoise_dim,VUKFF_upNoise_dim> MatrixUpCov;
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
		VUKFFilterState x_;
		/*! Estimate of covariance matrix */
		VUKFFilterState X_[1+2*(VUKFF_state_dim+VUKFF_preNoise_dim+VUKFF_upNoise_dim)];
		/*! Estimate of covariance matrix */
		MatrixP P_;
		/*! Sigma Samples of Prediction Noise*/
		Eigen::Matrix<double,VUKFF_preNoise_dim,2*VUKFF_preNoise_dim> PN_;
		/*! Sigma Samples of Update Noise*/
		Eigen::Matrix<double,VUKFF_upNoise_dim,2*VUKFF_upNoise_dim> UN_;
		/*! Contact flag counter */
		CF CFC_;
		/*! Legs used for kinematic update */
		SlippageDetection slippageDetection_;
		/*! Rotational rate estimate (bias corrected) */
		Eigen::Vector3d w_;
		/*! Current corrected accelerometer measurement */
		Eigen::Vector3d f_;
		/*! Innovation of update setp */
		Eigen::Matrix<double,12,1> y_;
		/*! Flag if Sigma points samples */
		bool mbSigmaSampled_;
	};

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

	/* -------------------- Filtering/Predicting/Updating --------------------- */
	/*! Filters the state x up to the given time t
	 * @param[in/out]	x	Filter state to be filtered
	 * @param[in]		t	Desired filter time
	 */
	void filterState(InternState& x,const double& t);
	/*! Predicts the state x using the given IMU measurement
	 * @param[in/out]	x	Filter state to be filtered
	 * @param[in]		t	Desired filter time
	 * @param[in]		m	IMU measurement
	 */
	void predictState(InternState& x,const double& t, const ImuMeas& m);
	/*! Updates the state x using the given Encoder measurement
	 * @param[in/out]	x	Filter state to be filtered
	 * @param[in]		m	Encoder measurement
	 */
	void encUpdateState(InternState& x,const EncMeas& m);
	/*! Concervative outlier detection based on predicted innovation covariance
	 * @param[in/out]	x		State
	 * @param[in]		Pyinv	Innovation information matrix
	 */
	void outlierDetection(InternState& x,const Eigen::Matrix<double,12,12>& Pyinv);
	/*! Makes and entry of the cuurent state into the log-file */
	void logState();

	void samplePredictionNoise(InternState& x,double dt);
	void sampleUpdateNoise(InternState& x);

	/* -------------------- P=rediction function --------------------- */
	void predict(VUKFFilterState& x,double dt,ImuMeas imuMeas);
	void predict(VUKFFilterState& x,double dt,ImuMeas imuMeas,Eigen::Matrix<double,VUKFF_preNoise_dim,1> n);

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
	/*! Predicition noise of position [m^2/s] (continuous form) */
	Eigen::Matrix3d Wr_;
	/*! Prediction noise of velocity */
	Eigen::Matrix3d Wv_;
	/*! Prediction noise of attitude */
	Eigen::Matrix3d Wq_;
	/*! Predicition noise of accelerometer bias [m^2/s^5] (continuous form */
	Eigen::Matrix3d Wbf_;
	/*! Predicition noise of gyroscope bias [rad^2/s^3] (continuous form */
	Eigen::Matrix3d Wbw_;
	/*! Threshold for kinematic outliers (there is an underlying chi-square distribution, dof=3) */
	double kinOutTh_;
	/*! Factor used during outlier restoration (must be larger than 1)*/
	double restorationFactor_;

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

	/* -------------------- Flag for modes (can be specified in the parameter file) --------------------- */
	/*! True if accelerometer bias is co-estimated */
	bool mbEstimateAccBias_;
	/*! True if gyroscope bias is co-estimated */
	bool mbEstimateRotBias_;
	/*! True if kinematic measurements are used */
	bool mbUseKin_;
	/*! True if using fixed timestepping (bug handling in SL)*/
	bool mbFixedTimeStepping_;
	/*! Timestep used if fixed time stepping*/
	double timeStep_;
};

}


#endif /* FilterVUKF2_HPP_ */
