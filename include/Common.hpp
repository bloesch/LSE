/*!
* @file 	Common.hpp
* @author 	Michael Blösch
* @date		10.10.2012
* @brief	Defines common stuff
 */

#ifndef LSE_COMMON_HPP_
#define LSE_COMMON_HPP_

#include "Rotations.hpp"
#include <vector>

#define LSE_MEAS_N 100000
#define LSE_DOF_LEG 3
#define LSE_N_LEG 4

#define LSE_VUKF_N 15

namespace LSE {

/*! Contact Flag Class */
class CF{
public:
	/*! Integer array containing contact flags */
	int F_[LSE_N_LEG];

	/*! Returns the number of contact points
	 * @return	number of contact points
	 */
	int sum(){
		int s = 0;
		for(int i=0;i<LSE_N_LEG;i++){
			if(F_[i]){
				s++;
			}
		}
		return s;
	}

	/* -------------------- Operator overloading --------------------- */
	/*! Assignement operator overloading */
	CF& operator= (const CF& cf) {
		for(int i=0;i<LSE_N_LEG;i++){
			F_[i] = cf.F_[i];
		}
		return *this;
	}
	/*! Element access operator overloading (const version) */
	const int& operator[](unsigned int i) const{ return F_[i];}
	/*! Element access operator overloading */
	int& operator[](unsigned int i) { return F_[i];}
};

/*! IMU measurement structure */
struct ImuMeas{
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW
	/*! Accelerometer measurement */
	Eigen::Vector3d f_;
	/*! Gyroscope measurement */
	Eigen::Vector3d w_;
};

/*! Encoder measurement structure */
struct EncMeas{
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW
	/*! Encoder position measurement */
	Eigen::Matrix<double,LSE_DOF_LEG,LSE_N_LEG> e_;
	/*! Encoder velocity measurement */
	Eigen::Matrix<double,LSE_DOF_LEG,LSE_N_LEG> v_;
	/*! Contact flag measurement */
	CF CF_;
};

/*! Pose sensor measurement structure */
struct PosMeas{
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW
	/*! Position measurement */
	Eigen::Vector3d r_;
	/*! Attitude measurement */
	Rotations::Quat q_;
};

/*! OF measurement structure */
struct OflMeas{
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW
	/*! Bearing vector */
	std::vector<Eigen::Vector3d> x_;
	/*! OF measurement */
	std::vector<Eigen::Vector3d> u_;
};

/*! State of robot main body */
struct State{
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW
	/*! Time */
	double t_;
	/*! Position */
	Eigen::Vector3d r_;
	/*! Velocity */
	Eigen::Vector3d v_;
	/*! Attitude */
	Rotations::Quat q_;
	/*! Rotational rate */
	Eigen::Vector3d w_;
	/*! Covariance matrix */
	Eigen::Matrix<double,12,12> P_;
};

}

#endif /* LSE_COMMON_HPP_ */
