/*
 * LieGroup.hpp
 *
 *  Created on: 28.07.2012
 *      Author: hannes
 */

#ifndef LIEGROUP_HPP_
#define LIEGROUP_HPP_

#include "DiffManifold.hpp"

namespace manifolds {

	template <int IDimension, int IPointSize, typename TScalar>
	struct LieGroupConf : public DiffManifoldConfiguration<IDimension, IPointSize, TScalar> {
		typedef LieGroupConf<IDimension, IPointSize, TScalar> Conf;
		typedef DiffManifoldConfiguration<IDimension, IPointSize, TScalar> ParentConf;
		typedef DiffManifold<Conf> Manifold;
	};

	namespace internal{
		template <typename TConfigurationDerived, int IDimension, int IPointSize, typename TScalar>
		struct DiffManifoldPointUpdateTraits<LieGroupConf<IDimension, IPointSize, TScalar>, TConfigurationDerived> {
			typedef manifolds::internal::DiffManifoldConfigurationTypeTrait<TConfigurationDerived> Types;
			typedef typename Types::point_t point_t;
			typedef typename Types::tangent_vector_t tangent_vector_t;
			typedef typename TConfigurationDerived::Manifold Manifold;
			static inline void update(const Manifold & manifold, point_t & point, const tangent_vector_t & vec){
				point_t p (point.rows());
				Manifold::expAtIdInto(vec, p);
				point = Manifold::mult(p, point);
				manifold.projectIntoManifold(point);
			}
		};
	}


	template <int IDimension, int IPointSize, typename TScalar, typename TConfigurationDerived>
	class DiffManifold<LieGroupConf<IDimension, IPointSize, TScalar>, TConfigurationDerived> : public DiffManifold<typename LieGroupConf<IDimension, IPointSize, TScalar>::ParentConf, TConfigurationDerived> {
	public:
		typedef DiffManifold<typename LieGroupConf<IDimension, IPointSize, TScalar>::ParentConf, TConfigurationDerived> parent_t;
		typedef TConfigurationDerived configuration_t;
		typedef internal::DiffManifoldConfigurationTypeTrait<configuration_t> Types;
		typedef typename Types::scalar_t scalar_t;
		typedef typename Types::point_t point_t;
		typedef typename Types::tangent_vector_t tangent_vector_t;
		typedef typename Types::dmatrix_t dmatrix_t;
		typedef typename configuration_t::Manifold base_t;
		typedef typename Eigen::Matrix<scalar_t, Types::Dimension, Types::PointSize> dlog_matrix_t;

		DiffManifold(TConfigurationDerived confiuration) : parent_t (confiuration) {}

		void getIdentityInto(point_t & p) const;
		point_t getIdentity() const;
		point_t getDefaultPoint() const;

		static void multInto(const point_t & a, const point_t & b, point_t & result);
		static point_t mult(const point_t & a, const point_t & b);

		static void expAtIdInto(const tangent_vector_t & vec, point_t & result);
		static void expInto(const point_t & point, const tangent_vector_t & vec, point_t & result);
		void dexpAtIdInto(const tangent_vector_t & vec, dmatrix_t & result) const;
		void dexpInto(const point_t & point, const tangent_vector_t & vec, dmatrix_t & result) const;

		point_t expAtId(const tangent_vector_t & vec) const;
		point_t exp(const point_t & point, const tangent_vector_t & vec) const;
		dmatrix_t dexpAtId(const tangent_vector_t & vec) const;
		dmatrix_t dexp(const point_t & point, const tangent_vector_t & vec) const;

		static void logAtIdInto(const point_t & to, tangent_vector_t & result);
		static void logInto(const point_t & from, const point_t & to, tangent_vector_t & result);
		void dlogInto(const point_t & point, const tangent_vector_t & vec, dmatrix_t & result) const;
		void dlogAtIdInto(const tangent_vector_t & vec, dmatrix_t & result) const;

		tangent_vector_t logAtId(const point_t & to) const;
		tangent_vector_t log(const point_t & from, const point_t & to) const;
		dlog_matrix_t dlogAtId(const point_t & to) const;
		dlog_matrix_t dlog(const point_t & from, const point_t & to) const;
	};
}

#include "implementation/LieGroupImpl.hpp"

#endif /* LIEGROUP_HPP_ */
