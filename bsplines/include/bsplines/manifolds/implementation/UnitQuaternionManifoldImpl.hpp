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
		_CLASS::DERIVED::logAtIdInto(qplus(quatInv(from), to), result);
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
		return sm::kinematics::quatV<TScalar>();
	}

	_TEMPLATE
	inline Eigen::Matrix3d _CLASS::S(const tangent_vector_t & vec){
		return expDiffMat(vec);
	}

	_TEMPLATE
	inline Eigen::Matrix3d _CLASS::LByVec(const tangent_vector_t & vec){
		return logDiffMat(vec);
	}

	_TEMPLATE
	void _CLASS::dlogInto(const point_t & from, const point_t & to, dmatrix_transposed_t & result) const {
		auto fromInv = quatInv(from);
		result = quatLogJacobian2(qplus(fromInv, to)) * quatPlus(fromInv);
	}

	_TEMPLATE
	void _CLASS::dlogAtIdInto(const point_t & to, dmatrix_transposed_t & result) const {
		result = quatLogJacobian2(to);
	}

	_TEMPLATE
	inline void _CLASS::dexpAtIdInto(const tangent_vector_t & vec, dmatrix_t & result) const
	{
		result = quatExpJacobian(vec);
	}

	_TEMPLATE
	inline void _CLASS::dexpInto(const point_t & point, const tangent_vector_t & vec, dmatrix_t & result) const
	{
		result = quatPlus(point) * this->dexpAtId(vec);
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
