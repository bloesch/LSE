/*!
* @file 	OptimizationFramework.cpp
* @author 	Michael Blösch
* @date		10.10.2012
 */

#include "OptimizationFramework.hpp"

namespace LSE {
namespace OF {

//---------------------------------------------------------------------------------------------------------------------------------------------------------//
//---------------------------------------------------------------------- Operations -----------------------------------------------------------------------//


double ScalarScalarAddition::f(const double& x1, const double& x2){
	return x1 + x2;
}

//static class ScalarVectorAddition: public Operation<Eigen::Vector3d,double,Eigen::Vector3d>{
//public:
//	ScalarVectorAddition();
//	~ScalarVectorAddition();
//	Eigen::Vector3d f(const double& x1, const Eigen::Vector3d& x2){
//		return x1*Eigen::Vector3d::Ones() + x2;
//	}
//} ScalarVectorAddition;
//
//static class VectorVectorAddition: public Operation<Eigen::Vector3d,Eigen::Vector3d,Eigen::Vector3d>{
//public:
//	VectorVectorAddition();
//	~VectorVectorAddition();
//	Eigen::Vector3d f(const Eigen::Vector3d& x1, const Eigen::Vector3d& x2){
//		return x1 + x2;
//	}
//} VectorVectorAddition;
//
////---------------------------------------------------------------------------------------------------------------------------------------------------------//
////---------------------------------------------------------------------- Expressions ----------------------------------------------------------------------//
//template<class T>
//class ExpressionBase{
//public:
//	ExpressionBase();
//	virtual ~ExpressionBase();
//	T x_;
//	virtual void fullEval();
//};
//
//template<class T,class T1,class T2>
//class Expression:public ExpressionBase<T>{
//public:
//	Expression(Operation<T,T1,T2>* op, ExpressionBase<T1>* exp1, ExpressionBase<T2>* exp2){
//		op_ = op;
//		exp1_ = exp1;
//		exp2_ = exp2;
//	}
//	Expression(){
//		op_ = NULL;
//		exp1_ = NULL;
//		exp2_ = NULL;
//	}
//	Expression(T x){
//		this->x_ = x;
//	}
//	~Expression();
//	void eval(){
//		this->x_ = op_->f(exp1_->x_,exp2_->x_);
//	}
//	void fullEval(){
//		exp1_->fullEval();
//		exp2_->fullEval();
//		this->x_ = op_->f(exp1_->x_,exp2_->x_);
//	}
//
//	Operation<T,T1,T2>* op_;
//	ExpressionBase<T1>* exp1_;
//	ExpressionBase<T2>* exp2_;
//};

}
}

