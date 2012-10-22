/*!
* @file 	Rotations.hpp
* @author 	Michael Blösch
* @date		10.10.2012
* @brief	Rotation stuff... convention:
* 			roll-pitch-yaw: 	alias
* 			rotation matrix: 	alibi
* 			rotation vector: 	alibi
* 			quaternion:			alibi
 */

#ifndef ROTATION_HPP_
#define ROTATION_HPP_

#include <Eigen/Dense>

namespace Rotations {

typedef Eigen::Vector4d Quat;

/*! Converts vector to sqew matrix
 * @return	corresponding sqew-matrix
 * @param[in] 	v	vector
 */
inline Eigen::Matrix3d vecToSqew(const Eigen::Vector3d& v){
	Eigen::Matrix3d M;
	M << 0, -v(2), v(1), v(2), 0, -v(0), -v(1), v(0), 0;
	return M;
}

/*! Limits the norm of a rotation vector to pi
 * @return 	rescaled rotation vector
 * @param[in] 	v	rotation vector
 */
inline Eigen::Vector3d rangePi(const Eigen::Vector3d& v){
	const double a = v.norm();
	if (a<=M_PI){
		return v;
	} else {
		const double a2 = -2.0*M_PI*floor((a+M_PI)/(2*M_PI))+a;
		return v/a*a2;
	}
}

/*! Converts a quaternion to a rotation matrix
 * @return 	corresponding rotation matrix
 * @param[in]	q	quaternion
 */
inline Eigen::Matrix3d quatToRotMat(const Quat& q){
	Eigen::Vector3d v = q.block(0,0,3,1);
	return (2*q(3)*q(3)-1)*Eigen::Matrix3d::Identity() + 2*q(3)*vecToSqew(v) + 2*v*v.transpose();
}

/*! Converts a quaternion to a rotation vector
 * @return 	corresponding rotation vector
 * @param[in]	q	quaternion
 */
inline Eigen::Vector3d quatToRotVec(const Quat& q){
	Eigen::Vector3d v;
	const double c = q(3);
	v = q.block(0,0,3,1);
	const double s = v.norm();
	if(s >= 1e-10){
		const double a = 2*atan2(s,c);
		return v*a/s;
	} else {
		return v*2;
	}
}

/*! Converts a rotation vector to a quaternion
 * @return 	corresponding quaternion
 * @param[in]	v	rotation vector
 */
inline Quat rotVecToQuat(const Eigen::Vector3d& v){
	Quat q;
	const double a = v.norm();
	q(3) = cos(a/2);
	if(a >= 1e-10){
		q.block(0,0,3,1) = sin(a/2)/a*v;
	} else {
		q.block(0,0,3,1) = v;
	}
	q.normalize();
	return q;
}

/*! Computes the inverse of a quaternion
 * @return 	corresponding quaternion inverse
 * @param[in]	q	quaternion
 */
inline Quat quatInverse(const Quat& q){
	Quat q2;
	q2.block(0,0,3,1) = -q.block(0,0,3,1);
	q2(3) = q(3);
	return q2;
}

/*! Returns the identity quaternion
 * @return 	identity quaternion
 */
inline Quat quatIdentity(){
	Quat q;
	q.setZero();
	q(3) = 1;
	return q;
}

/*! Computes the left-hand multiplication matrix from a given quaternion
 * @return 	left-hand multiplication matrix
 * @param[in]	q	quaternion
 */
inline Eigen::Matrix<double,4,4> quatL(const Quat& q){
	Eigen::Matrix<double,4,4> M;
	M.setIdentity();
	M = M*q(3);
	M.block(0,0,3,3) += vecToSqew(q.block(0,0,3,1));
	M.block(3,0,1,4) = -q.transpose();
	M.block(0,3,4,1) = q;
	return M;
}

/*! Computes the right-hand multiplication matrix from a given quaternion
 * @return 	right-hand multiplication matrix
 * @param[in]	q	quaternion
 */
inline Eigen::Matrix<double,4,4> quatR(const Quat& q){
	Eigen::Matrix<double,4,4> M;
	M.setIdentity();
	M = M*q(3);
	M.block(0,0,3,3) -= vecToSqew(q.block(0,0,3,1));
	M.block(3,0,1,4) = -q.transpose();
	M.block(0,3,4,1) = q;
	return M;
}

//static void quatToRpy(const Eigen::Quaterniond& q, Eigen::Vector3d& rpy){
//	const double q0 = q.w();
//	const double q1 = q.x();
//	const double q2 = q.y();
//	const double q3 = q.z();
//	rpy(0) = -atan2(2*(q0*q1+q2*q3),1-2*(pow(q1,2)+pow(q2,2)));
//	rpy(1) = -asin(2*(q0*q2-q1*q3));
//	rpy(2) = -atan2(2*(q0*q3+q1*q2),1-2*(pow(q2,2)+pow(q3,2)));
//}
//
//static void rpyToQuat(const Eigen::Vector3d& rpy, Eigen::Quaterniond& q){
//	const double cy = cos(rpy(2)/2);
//	const double sy = sin(rpy(2)/2);
//	const double cc = cos(rpy(1)/2)*cos(rpy(0)/2);
//	const double cs = cos(rpy(1)/2)*sin(rpy(0)/2);
//	const double sc = sin(rpy(1)/2)*cos(rpy(0)/2);
//	const double ss = sin(rpy(1)/2)*sin(rpy(0)/2);
//
//	q.w() = cy*cc-sy*ss;
//	q.x() = -cy*cs-sy*sc;
//	q.y() = -cy*sc+sy*cs;
//	q.z() = -cy*ss-sy*cc;
//	q.normalize();
//}
//
//static void rpyToEar(const Eigen::Vector3d& rpy, Eigen::Matrix3d& ear){
//	const double cp = cos(rpy(1));
//	const double sp = sin(rpy(1));
//	const double cy = cos(rpy(2));
//	const double sy = sin(rpy(2));
//	ear << cp*cy, sy, 0, -cp*sy, cy, 0, sp, 0, 1;
//}
//
//static void rpyToEarInv(const Eigen::Vector3d& rpy, Eigen::Matrix3d& earInv){
//	const double t2 = cos(rpy(1));
//	const double t3 = 1.0/t2;
//	const double t4 = sin(rpy(2));
//	const double t5 = cos(rpy(2));
//	const double t6 = tan(rpy(1));
//	earInv(0,0) = t3*t5;
//	earInv(0,1) = -t3*t4;
//	earInv(0,2) = 0;
//	earInv(1,0) = t4;
//	earInv(1,1) = t5;
//	earInv(1,2) = 0;
//	earInv(2,0) = -t5*t6;
//	earInv(2,1) = t4*t6;
//	earInv(2,2) = 1;
//}
//
//static void rotMatToQuat(const Eigen::Matrix3d& mat, Eigen::Quaterniond& q){
//	// Bad precision!
//	double w;
//	double x;
//	double y;
//	double z;
//	w = sqrt(1+mat(0,0)+mat(1,1)+mat(2,2))/2;
//	if(w>0.02){
//		x = 1/4/w*(mat(2,1)-mat(1,2));
//		y = 1/4/w*(mat(0,2)-mat(2,0));
//		z = 1/4/w*(mat(1,0)-mat(0,1));
//	} else {
//		x = sqrt(1+mat(0,0)-mat(1,1)-mat(2,2))/2;
//		y = 1/4/x*(mat(0,1)+mat(1,0));
//		z = 1/4/x*(mat(0,2)+mat(2,0));
//		w = 1/4/x*(mat(2,1)-mat(1,2));
//	}
//	q.w() = w;
//	q.x() = x;
//	q.y() = y;
//	q.z() = z;
//	q.normalize();
//}
//
//static void rotMatToRotVec(const Eigen::Matrix3d& mat, Eigen::Vector3d& vec){
//	double cosTheta = 0.5*(mat(0,0)+mat(1,1)+mat(2,2)-1);
//	if(cosTheta > 1.0){
//		cosTheta = 1;
//	} else if(cosTheta < -1.0){
//		cosTheta = -1.0;
//	}
//	double theta = acos(cosTheta);
//	vec(0) = (mat(2,1)-mat(1,2));
//	vec(1) = (mat(0,2)-mat(2,0));
//	vec(2) = (mat(1,0)-mat(0,1));
//	if(std::abs(theta) > 1e-10){
//		vec *= 0.5/sin(theta)*theta;
//	} else {
//		vec *= 0.5;
//	}
//}
//
//
//static void rpyToRotMat(const Eigen::Vector3d& rpy, Eigen::Matrix3d& mat){
//	const double t1 = cos(rpy(2));
//	const double t2 = sin(rpy(0));
//	const double t3 = sin(rpy(2));
//	const double t4 = cos(rpy(0));
//	const double t5 = sin(rpy(1));
//	const double t6 = cos(rpy(1));
//	mat(0,0) = t1*t6;
//	mat(0,1) = t3*t4+t1*t2*t5;
//	mat(0,2) = t2*t3-t1*t4*t5;
//	mat(1,0) = -t3*t6;
//	mat(1,1) = t1*t4-t2*t3*t5;
//	mat(1,2) = t1*t2+t3*t4*t5;
//	mat(2,0) = t5;
//	mat(2,1) = -t2*t6;
//	mat(2,2) = t4*t6;
//}

}

#endif /* ROTATION_HPP_ */
