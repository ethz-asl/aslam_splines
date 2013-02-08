
#include "bsplines/manifolds/LieGroup.hpp"

namespace manifolds {
	template<int IDimension, int IPointSize, typename TScalar, typename TConfigurationDerived>
	inline typename DiffManifold< LieGroupConf<IDimension, IPointSize, TScalar>, TConfigurationDerived>::point_t DiffManifold< LieGroupConf<IDimension, IPointSize, TScalar>, TConfigurationDerived>::
	getIdentity() const
	{
		point_t p(this->getPointSize());
		((DiffManifold< TConfigurationDerived, TConfigurationDerived> * )this)->getIdentityInto(p);
		return p;
	}

	template<int IDimension, int IPointSize, typename TScalar, typename TConfigurationDerived>
	inline typename DiffManifold< LieGroupConf<IDimension, IPointSize, TScalar>, TConfigurationDerived>::point_t DiffManifold< LieGroupConf<IDimension, IPointSize, TScalar>, TConfigurationDerived>::
	getDefaultPoint() const
	{
		return getIdentity();
	}

	template<int IDimension, int IPointSize, typename TScalar, typename TConfigurationDerived>
	inline typename DiffManifold< LieGroupConf<IDimension, IPointSize, TScalar>, TConfigurationDerived>::point_t DiffManifold< LieGroupConf<IDimension, IPointSize, TScalar>, TConfigurationDerived>::
	mult(const point_t & a, const point_t & b)
	{
		point_t p(a.rows());
		DiffManifold< TConfigurationDerived, TConfigurationDerived>::multInto(a, b, p);
		return p;
	}

	template<int IDimension, int IPointSize, typename TScalar, typename TConfigurationDerived>
	inline typename DiffManifold< LieGroupConf<IDimension, IPointSize, TScalar>, TConfigurationDerived>::point_t DiffManifold< LieGroupConf<IDimension, IPointSize, TScalar>, TConfigurationDerived>::
	exp(const point_t & point, const tangent_vector_t & vec) const
	{
		point_t p(point.rows());
		((DiffManifold< TConfigurationDerived, TConfigurationDerived> * )this)->expInto(point, vec, p);
		return p;
	}

	template<int IDimension, int IPointSize, typename TScalar, typename TConfigurationDerived>
	inline typename DiffManifold< LieGroupConf<IDimension, IPointSize, TScalar>, TConfigurationDerived>::point_t DiffManifold< LieGroupConf<IDimension, IPointSize, TScalar>, TConfigurationDerived>::
	expAtId(const tangent_vector_t & vec) const
	{
		point_t result(this->getPointSize());
		((DiffManifold< TConfigurationDerived, TConfigurationDerived> * )this)->expAtIdInto(vec, result);
		return result;
	}

	template<int IDimension, int IPointSize, typename TScalar, typename TConfigurationDerived>
	inline typename DiffManifold< LieGroupConf<IDimension, IPointSize, TScalar>, TConfigurationDerived>::dmatrix_t DiffManifold< LieGroupConf<IDimension, IPointSize, TScalar>, TConfigurationDerived>::
	dexpAtId(const tangent_vector_t & vec) const
	{
		dmatrix_t result(this->getPointSize(), this->getDimension());
		((DiffManifold< TConfigurationDerived, TConfigurationDerived> * )this)->dexpAtIdInto(vec, result);
		return result;
	}

	template<int IDimension, int IPointSize, typename TScalar, typename TConfigurationDerived>
	inline typename DiffManifold< LieGroupConf<IDimension, IPointSize, TScalar>, TConfigurationDerived>::tangent_vector_t DiffManifold< LieGroupConf<IDimension, IPointSize, TScalar>, TConfigurationDerived>::
	log(const point_t & from, const point_t & to) const
	{
		tangent_vector_t v(this->getDimension());
		((DiffManifold< TConfigurationDerived, TConfigurationDerived> * )this)->logInto(from, to, v);
		return v;
	}

	template<int IDimension, int IPointSize, typename TScalar, typename TConfigurationDerived>
	inline typename DiffManifold< LieGroupConf<IDimension, IPointSize, TScalar>, TConfigurationDerived>::tangent_vector_t DiffManifold< LieGroupConf<IDimension, IPointSize, TScalar>, TConfigurationDerived>::
	logAtId(const point_t & to) const
	{
		tangent_vector_t v(this->getDimension());
		((DiffManifold< TConfigurationDerived, TConfigurationDerived> * )this)->logAtIdInto(to, v);
		return v;
	}

	template<int IDimension, int IPointSize, typename TScalar, typename TConfigurationDerived>
	inline typename DiffManifold< LieGroupConf<IDimension, IPointSize, TScalar>, TConfigurationDerived>::dlog_matrix_t DiffManifold< LieGroupConf<IDimension, IPointSize, TScalar>, TConfigurationDerived>::
	dlogAtId(const point_t & to) const
	{
		dlog_matrix_t result(this->getDimension(), this->getDimension());
		((DiffManifold< TConfigurationDerived, TConfigurationDerived> * )this)->dlogAtIdInto(to, result);
		return result;
	}
}
