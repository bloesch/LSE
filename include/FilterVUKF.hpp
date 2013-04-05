/*!
* @file 	FilterVUKF.hpp
* @author 	Michael Blösch
* @date		10.10.2012
* @brief	VUKF Filter for legged robots
 */

#ifndef FILTERVUKF_HPP_
#define FILTERVUKF_HPP_

#include "FilterBase.hpp"
#include "Common.hpp"
#include <Eigen/Dense>
#include "Rotations.hpp"

namespace LSE {

class Manager;
class StateVUKF;

/*! Observability Constrained Extended Kalman Filter */
class FilterVUKF: public FilterBase{
public:public:
	/* -------------------- Constructor/Destructor --------------------- */
	/*! Constructor
	 * @param[in]	pManager	pointer to main class Manager
	 * @param[in]	pFilename	filename of parameter-file
	 */
	FilterVUKF(Manager* pManager,const char* pFilename);
	/*! Destructor */
	virtual ~FilterVUKF();

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
	/*! Resets the filter
	 * @param[in]	t	time used to initialize new state estimate
	 */
	virtual void resetEstimate(const double& t);
	/*! Returns a string describing the main filter parameters
	 * @param[out] str	string characterize the parameter set of the filter
	 */
	virtual std::string getKeyString();

private:
	typedef Eigen::Matrix<double,30,30> Matrix30d;
	/*! Loads overall parameters from parameter file
	 * @param[in]	pFilename	name of parameter file
	 */
	void loadParam(const char* pFilename);
	/*! Structure of filter intern state augmented with process noise*/
	struct AugmentedState{
		EIGEN_MAKE_ALIGNED_OPERATOR_NEW

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
		/*! Noise on position prediction */
		Eigen::Vector3d nr_;
		/*! Accelerometer measurement (noise affected, not bias corrected) */
		Eigen::Vector3d f_;
		/*! Rotational rate measurement (noise affected, not bias corrected) */
		Eigen::Vector3d w_;
		/*! Random walk of accelerometer bias */
		Eigen::Vector3d nbf_;
		/*! Random walk of gyroscope bias */
		Eigen::Vector3d nbw_;
	};
	/*! Structure of filter intern state */
	struct InternState{
		EIGEN_MAKE_ALIGNED_OPERATOR_NEW

		/*! Time of estimate */
		double t_;
		/*! State estimate */
		AugmentedState x_;
		/*! Sigma points of state */
		AugmentedState X_[1+4*(LSE_VUKF_N)];
		/*! Flag if Sigma points samples */
		bool mbSigmaSampled_;
		/*! Contact flag counter */
		CF CFC_;
		/*! Legs used for kinematic update */
		CF LegArray_;
		/*! Estimate of covariance matrix */
		Eigen::Matrix<double,15,15> P_;
		/*! Innovation of update setp */
		Eigen::Matrix<double,12,1> y_;

		/* -------------------- Operator overloading --------------------- */
		/*! Assignement operator overloading */
		InternState& operator= (const InternState& x) {
			t_ = x.t_;
			x_ = x.x_;
			for(int i=0;i<1+4*(LSE_VUKF_N);i++){
				X_[i] = x.X_[i];
			}
			CFC_ = x.CFC_;
			LegArray_ = x.LegArray_;
			P_ = x.P_;
			return *this;
		}
	};

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

class StateVUKF{
public:
	Eigen::Matrix<double,3,1>& r(){
		return x_[0];
	}
	Eigen::Matrix<double,3,1>& v(){
		return x_[1];
	}
	Eigen::Matrix<double,3,1>& b(){
		return x_[2];
	}
	Eigen::Matrix<double,3,1>& c(){
		return x_[3];
	}
	Rotations::Quat& q(){
		return x_(0);
	}


private:
	Rotations::LieG<0,4,1> x_;
};

class PreNoiseVUKF{
public:
	Eigen::Matrix<double,3,1>& nr(){
		return n_[0];
	}
	Eigen::Matrix<double,3,1>& nf(){
		return n_[1];
	}
	Eigen::Matrix<double,3,1>& nw(){
		return n_[2];
	}
	Eigen::Matrix<double,3,1>& nb(){
		return n_[3];
	}
	Eigen::Matrix<double,3,1>& nc(){
		return n_[4];
	}


private:
	Rotations::LieG<0,5,0> n_;

};

class MeasKinNoiseVUKF{
public:
	Eigen::Matrix<double,3,1>& n(int i){
		return n_[i];
	}


private:
	Rotations::LieG<0,4,0> n_;

};

}


#endif /* FILTERVUKF_HPP_ */
