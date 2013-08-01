/*!
* @file 	OptimizationFramework.hpp
* @author 	Michael Blösch
* @date		10.10.2012
* @brief	Rotation stuff... convention:
* 			roll-pitch-yaw: 	alias
* 			rotation matrix: 	alibi
* 			rotation vector: 	alibi
* 			quaternion:			alibi
 */

#ifndef LSE_OPTIMIZATIONFRAMEWORK_HPP_
#define LSE_OPTIMIZATIONFRAMEWORK_HPP_

#include <Eigen/Dense>
#include "Rotations.hpp"

namespace LSE {
namespace OF {

//---------------------------------------------------------------------------------------------------------------------------------------------------------//
//---------------------------------------------------------------------- Operations -----------------------------------------------------------------------//
template<class T,class T1,class T2>
class Operation{
public:
	Operation(){};
	~Operation(){};
	virtual T f(const T1& x1, const T2& x2) = 0;
};

template<>
class Operation<double,double,double>{
public:
	Operation(){};
	~Operation(){};
	virtual double f(const double& x1, const double& x2) = 0;
	virtual double J1(const double& x1, const double& x2) = 0;
	virtual double J2(const double& x1, const double& x2) = 0;
};

template<>
class Operation<double,Eigen::Vector3d,double>{
public:
	Operation(){};
	~Operation(){};
	virtual double f(const Eigen::Vector3d& x1, const double& x2) = 0;
	virtual Eigen::Matrix<double,1,3> J1(const Eigen::Vector3d& x1, const double& x2) = 0;
	virtual double J2(const Eigen::Vector3d& x1, const double& x2) = 0;
};

template<>
class Operation<double,double,Eigen::Vector3d>{
public:
	Operation(){};
	~Operation(){};
	virtual double f(const double& x1, const Eigen::Vector3d& x2) = 0;
	virtual double J1(const double& x1, const Eigen::Vector3d& x2) = 0;
	virtual Eigen::Matrix<double,1,3> J2(const double& x1, const Eigen::Vector3d& x2) = 0;
};

template<>
class Operation<double,Eigen::Vector3d,Eigen::Vector3d>{
public:
	Operation(){};
	~Operation(){};
	virtual double f(const Eigen::Vector3d& x1, const Eigen::Vector3d& x2) = 0;
	virtual Eigen::Matrix<double,1,3> J1(const Eigen::Vector3d& x1, const Eigen::Vector3d& x2) = 0;
	virtual Eigen::Matrix<double,1,3> J2(const Eigen::Vector3d& x1, const Eigen::Vector3d& x2) = 0;
};

template<>
class Operation<Eigen::Vector3d,double,double>{
public:
	Operation(){};
	~Operation(){};
	virtual Eigen::Vector3d f(const double& x1, const double& x2) = 0;
	virtual Eigen::Vector3d J1(const double& x1, const double& x2) = 0;
	virtual Eigen::Vector3d J2(const double& x1, const double& x2) = 0;
};

template<>
class Operation<Eigen::Vector3d,Eigen::Vector3d,double>{
public:
	Operation(){};
	~Operation(){};
	virtual Eigen::Vector3d f(const Eigen::Vector3d& x1, const double& x2) = 0;
	virtual Eigen::Matrix<double,3,3> J1(const Eigen::Vector3d& x1, const double& x2) = 0;
	virtual Eigen::Vector3d J2(const Eigen::Vector3d& x1, const double& x2) = 0;
};

template<>
class Operation<Eigen::Vector3d,double,Eigen::Vector3d>{
public:
	Operation(){};
	~Operation(){};
	virtual Eigen::Vector3d f(const double& x1, const Eigen::Vector3d& x2) = 0;
	virtual Eigen::Vector3d J1(const double& x1, const Eigen::Vector3d& x2) = 0;
	virtual Eigen::Matrix<double,3,3> J2(const double& x1, const Eigen::Vector3d& x2) = 0;
};

template<>
class Operation<Eigen::Vector3d,Eigen::Vector3d,Eigen::Vector3d>{
public:
	Operation(){};
	~Operation(){};
	virtual Eigen::Vector3d f(const Eigen::Vector3d& x1, const Eigen::Vector3d& x2) = 0;
	virtual Eigen::Matrix<double,3,3> J1(const Eigen::Vector3d& x1, const Eigen::Vector3d& x2) = 0;
	virtual Eigen::Matrix<double,3,3> J2(const Eigen::Vector3d& x1, const Eigen::Vector3d& x2) = 0;
};

static class ScalarScalarAddition: public Operation<double,double,double>{
public:
	ScalarScalarAddition(){};
	~ScalarScalarAddition(){};
	double f(const double& x1, const double& x2){
		return x1 + x2;
	}
	double J1(const double& x1, const double& x2){
		return 1.0;
	}
	double J2(const double& x1, const double& x2){
		return 1.0;
	}
} ScalarScalarAddition;

static class ScalarVectorAddition: public Operation<Eigen::Vector3d,double,Eigen::Vector3d>{
public:
	ScalarVectorAddition(){};
	~ScalarVectorAddition(){};
	Eigen::Vector3d f(const double& x1, const Eigen::Vector3d& x2){
		return x1*Eigen::Vector3d::Ones() + x2;
	}
	Eigen::Vector3d J1(const double& x1, const Eigen::Vector3d& x2){
		return Eigen::Vector3d::Ones();
	}
	Eigen::Matrix<double,3,3> J2(const double& x1, const Eigen::Vector3d& x2){
		return Eigen::Matrix<double,3,3>::Identity();
	}
} ScalarVectorAddition;

static class VectorVectorAddition: public Operation<Eigen::Vector3d,Eigen::Vector3d,Eigen::Vector3d>{
public:
	VectorVectorAddition(){};
	~VectorVectorAddition(){};
	Eigen::Vector3d f(const Eigen::Vector3d& x1, const Eigen::Vector3d& x2){
		return x1 + x2;
	}
	Eigen::Matrix<double,3,3> J1(const Eigen::Vector3d& x1, const Eigen::Vector3d& x2){
		return Eigen::Matrix<double,3,3>::Identity();
	}
	Eigen::Matrix<double,3,3> J2(const Eigen::Vector3d& x1, const Eigen::Vector3d& x2){
		return Eigen::Matrix<double,3,3>::Identity();
	}
} VectorVectorAddition;

static class ScalarScalarMult: public Operation<double,double,double>{
public:
	ScalarScalarMult(){};
	~ScalarScalarMult(){};
	double f(const double& x1, const double& x2){
		return x1*x2;
	}
	double J1(const double& x1, const double& x2){
		return x2;
	}
	double J2(const double& x1, const double& x2){
		return x1;
	}
} ScalarScalarMult;

static class ScalarVectorMult: public Operation<Eigen::Vector3d,double,Eigen::Vector3d>{
public:
	ScalarVectorMult(){};
	~ScalarVectorMult(){};
	Eigen::Vector3d f(const double& x1, const Eigen::Vector3d& x2){
		return x1*x2;
	}
	Eigen::Vector3d J1(const double& x1, const Eigen::Vector3d& x2){
		return x2;
	}
	Eigen::Matrix<double,3,3> J2(const double& x1, const Eigen::Vector3d& x2){
		return x1*Eigen::Matrix<double,3,3>::Identity();
	}
} ScalarVectorMult;

static class CrossProduct: public Operation<Eigen::Vector3d,Eigen::Vector3d,Eigen::Vector3d>{
public:
	CrossProduct(){};
	~CrossProduct(){};
	Eigen::Vector3d f(const Eigen::Vector3d& x1, const Eigen::Vector3d& x2){
		return x1.cross(x2); //TODO Check!
	}
	Eigen::Matrix<double,3,3> J1(const Eigen::Vector3d& x1, const Eigen::Vector3d& x2){
		return -Rotations::vecToSqew(x2);
	}
	Eigen::Matrix<double,3,3> J2(const Eigen::Vector3d& x1, const Eigen::Vector3d& x2){
		return Rotations::vecToSqew(x1);
	}
} CrossProduct;

static class DotProduct: public Operation<double,Eigen::Vector3d,Eigen::Vector3d>{
public:
	DotProduct(){};
	~DotProduct(){};
	double f(const Eigen::Vector3d& x1, const Eigen::Vector3d& x2){
		return x1.dot(x2); //TODO Check!
	}
	Eigen::Matrix<double,1,3> J1(const Eigen::Vector3d& x1, const Eigen::Vector3d& x2){
		return x2.transpose();
	}
	Eigen::Matrix<double,1,3> J2(const Eigen::Vector3d& x1, const Eigen::Vector3d& x2){
		return x1.transpose();
	}
} DotProduct;

//---------------------------------------------------------------------------------------------------------------------------------------------------------//
//---------------------------------------------------------------------- Expressions ----------------------------------------------------------------------//
template<class T>
class ExpressionBase{
public:
	ExpressionBase(){};
	~ExpressionBase(){};
	T x_;
	virtual void fullEval() = 0;
};

template<class T,class T1,class T2>
class Expression:public ExpressionBase<T>{
public:
	Expression(Operation<T,T1,T2>* op, ExpressionBase<T1>* exp1, ExpressionBase<T2>* exp2){
		op_ = op;
		exp1_ = exp1;
		exp2_ = exp2;
	}
	Expression(){
		op_ = NULL;
		exp1_ = NULL;
		exp2_ = NULL;
	}
	Expression(T x){
		op_ = NULL;
		exp1_ = NULL;
		exp2_ = NULL;
		this->x_ = x;
	}
	~Expression(){};
	void eval(){
		if(op_ != NULL){
			this->x_ = op_->f(exp1_->x_,exp2_->x_);
		}
	}
	void fullEval(){
		if(op_ != NULL){
			if(exp1_ != NULL) exp1_->fullEval();
			if(exp2_ != NULL) exp2_->fullEval();
			this->x_ = op_->f(exp1_->x_,exp2_->x_);
		};
	}

	Operation<T,T1,T2>* op_;
	ExpressionBase<T1>* exp1_;
	ExpressionBase<T2>* exp2_;
};


}
}

#endif /* LSE_OPTIMIZATIONFRAMEWORK_HPP_ */
