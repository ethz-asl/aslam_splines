/*
* UnitQuaternionBSpline.hpp
*
*  Created on: May 10, 2012
*      Author: hannes
*/


#include <sm/kinematics/quaternion_algebra.hpp>
#include <sm/kinematics/rotations.hpp>

#include "bsplines/manifolds/UnitQuaternionManifold.hpp"


namespace manifolds {

#define _TEMPLATE template <typename TScalar, typename TConfigurationDerived>
#define _CLASS DiffManifold<UnitQuaternionManifoldConf<TScalar>, TConfigurationDerived>

	using namespace sm::kinematics;

	_TEMPLATE
	const TScalar _CLASS::epsilon6thRoot=pow(std::numeric_limits<TScalar>::epsilon(), 1.0/6);

	_TEMPLATE
	inline void _CLASS::multInto(const point_t & a, const point_t & b, point_t & result)
	{
		result = qplus(a, b);
	}

	_TEMPLATE
	inline void _CLASS::logAtIdInto(const point_t & to, tangent_vector_t & result)
	{
		result = quat2AxisAngle(to);
	}


	_TEMPLATE
	inline void _CLASS::logInto(const point_t & from, const point_t & to, tangent_vector_t & result)
	{
		logAtIdInto(qplus(quatInv(from), to), result);
	}

	_TEMPLATE
	inline typename _CLASS::point_t _CLASS::expAtId(const tangent_vector_t & vec)
	{
		return axisAngle2quat(vec);
	}

	_TEMPLATE
	inline void _CLASS::expAtIdInto(const tangent_vector_t & vec, point_t & result)
	{
		result = expAtId(vec);
	}

	_TEMPLATE
	inline void _CLASS::expInto(const point_t & point, const tangent_vector_t & vec, point_t & result)
	{
		result = qplus(point, expAtId(vec));
	}

	_TEMPLATE
	inline void _CLASS::getIdentityInto(point_t & result) {
		result = quatIdentity();
	}

	_TEMPLATE
	inline const typename _CLASS::dmatrix_t & _CLASS::V(){
		static const dmatrix_t V = 0.5 * dmatrix_t::Identity();
		return V;
	}

	_TEMPLATE
	inline Eigen::Matrix3d _CLASS::S(const tangent_vector_t & vec){
		TScalar phi = vec.norm();

		if(phi == 0){
			return Eigen::Matrix3d::Identity();
		}

		Eigen::Matrix3d vecCross = crossMx(vec);

		TScalar phiAbs = fabs(phi);
		TScalar phiSquare = phi * phi;

		TScalar a;
		if(phiAbs > epsilon6thRoot){
			TScalar siPhiHalf = sin(phi / 2);
			a = (2 * siPhiHalf * siPhiHalf / phiSquare);
		}
		else{
			a = (1.0/2) * (1 - (1.0 / (24 / 2)) * phiSquare * (1 - (1.0 / (720 / 24)) * phiSquare));
		}

		TScalar b;
		if(phiAbs > epsilon6thRoot){
			b = ((1 - sin(phi) / phi)/phiSquare);
		}
		else{
			b = (1.0/6) * (1 - (1.0 / (120 / 6)) * phiSquare * (1 - (1.0 / (5040 / 120)) * phiSquare));
		}

		return Eigen::Matrix3d::Identity() - a * vecCross + b * vecCross * vecCross;
	}

	_TEMPLATE
	inline Eigen::Matrix3d _CLASS::LByVec(const tangent_vector_t & vec){
		TScalar phi = vec.norm();
		if(phi == 0){
			return Eigen::Matrix3d::Identity();
		}

		TScalar phiAbs = fabs(phi);
		Eigen::Matrix3d vecCross = crossMx(vec);

		TScalar a;
		if(phiAbs > epsilon6thRoot){
			TScalar phiHalf = 0.5 * phi;
			a = ((1 - phiHalf / tan(phiHalf))/phi/phi);
		}
		else{
			TScalar phiSquare = phi * phi;
			a = 1.0 / 12 * (1 + 1.0 / 60 * phiSquare * (1 + 1.0/(30240 / 720) * phiSquare));
		}
		return Eigen::Matrix3d::Identity() + 0.5 * vecCross + a * vecCross * vecCross;
	}

//TODO implement: void dlog(const point_t & point, const tangent_vector_t & vec, dmatrix_t & result){
//		result = quatPlus(point) * quatOPlus(axisAngle2quat(vec)) * V() * S(vec);
//	}

	_TEMPLATE
	inline void _CLASS::dexpAtIdInto(const tangent_vector_t & vec, dmatrix_t & result) const
	{
		result = quatOPlus(axisAngle2quat(vec)) * _CLASS::V() * _CLASS::S(vec);
	}

	_TEMPLATE
	inline void _CLASS::dexpInto(const point_t & point, const tangent_vector_t & vec, dmatrix_t & result) const
	{
		result = quatPlus(point) * dexpAtId(vec);
	}

	_TEMPLATE
	bool _CLASS::isInManifold(const point_t & pt)
	{
		return fabs(pt.norm() - 1) < fabs(1E-9);
	}

	_TEMPLATE
	void _CLASS::projectIntoManifold(point_t & pt)
	{
		pt /= pt.norm();
	}

	_TEMPLATE
	void _CLASS::randomizePoint(point_t & pt)
	{
		double norm;
		do{ //TODO improve : realize uniform distribution on S3
			pt = 5*(Eigen::Vector4d::Random() * 2 - Eigen::Vector4d::Ones());
		}while((norm = pt.norm()) < 1e-10);
		pt /= norm;
	}

#undef _TEMPLATE
#undef _CLASS
}
