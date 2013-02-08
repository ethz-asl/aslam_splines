/*
 * EuclideanSpaceImpl.hpp
 *
 *  Created on: Jul 26, 2012
 *      Author: hannes
 */


#include "bsplines/manifolds/EuclideanSpace.hpp"

namespace manifolds {

	template <int IDimension, typename TScalar, typename TConfigurationDerived>
	inline void DiffManifold<EuclideanSpaceConf<IDimension, TScalar>, TConfigurationDerived>::
	getIdentityInto(point_t & pt) const {
		pt = point_t::Zero(this->getPointSize());
	}

	template <int IDimension, typename TScalar, typename TConfigurationDerived>
	inline void DiffManifold<EuclideanSpaceConf<IDimension, TScalar>, TConfigurationDerived>::
	multInto(const point_t & a, const point_t & b, point_t & result)
	{
		result = a + b;
	}

	template <int IDimension, typename TScalar, typename TConfigurationDerived>
	inline void DiffManifold<EuclideanSpaceConf<IDimension, TScalar>, TConfigurationDerived>::
	logInto(const point_t & from, const point_t & to, tangent_vector_t & result) {
		result = to - from;
	}

	template <int IDimension, typename TScalar, typename TConfigurationDerived>
	inline void DiffManifold<EuclideanSpaceConf<IDimension, TScalar>, TConfigurationDerived>::
	expInto(const point_t & point, const tangent_vector_t & vec, point_t & result) {
		result = point + vec;
	}

	template <int IDimension, typename TScalar, typename TConfigurationDerived>
	inline void DiffManifold<EuclideanSpaceConf<IDimension, TScalar>, TConfigurationDerived>::
	expAtIdInto(const tangent_vector_t & vec, point_t & p) {
		p = vec;
	}

	template <int IDimension, typename TScalar, typename TConfigurationDerived>
	inline void DiffManifold<EuclideanSpaceConf<IDimension, TScalar>, TConfigurationDerived>::
	dexpAtIdInto(const tangent_vector_t & vec, dmatrix_t & result) const {
		result = dmatrix_t::Identity(result.rows(), result.cols());
	}

	template <int IDimension, typename TScalar, typename TConfigurationDerived>
	inline void DiffManifold<EuclideanSpaceConf<IDimension, TScalar>, TConfigurationDerived>::
	dexpInto(const point_t & from, const tangent_vector_t & vec, dmatrix_t & result) const {
		result = dmatrix_t::Identity(result.rows(), result.cols());
	}

	template <int IDimension, typename TScalar, typename TConfigurationDerived>
	inline void DiffManifold<EuclideanSpaceConf<IDimension, TScalar>, TConfigurationDerived>::
	randomizePoint(point_t & pt) const {
		pt = 2 * point_t::Random(this->getDimension()) - point_t::Ones(this->getDimension());
	}
}
