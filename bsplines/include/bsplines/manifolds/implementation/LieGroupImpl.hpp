
#include "bsplines/manifolds/LieGroup.hpp"

namespace manifolds {
#define _TEMPLATE template<int IDimension, int IPointSize, typename TScalar, typename TConfigurationDerived>
#define _CLASS DiffManifold< LieGroupConf<IDimension, IPointSize, TScalar>, TConfigurationDerived>

_TEMPLATE
inline typename _CLASS::point_t _CLASS::getIdentity() const
{
	point_t p(this->getPointSize());
	this->getDerived().getIdentityInto(p);
	return p;
}

_TEMPLATE
inline typename _CLASS::point_t _CLASS::getDefaultPoint() const
{
	return getIdentity();
}

_TEMPLATE
inline typename _CLASS::point_t _CLASS::mult(const point_t & a, const point_t & b)
{
	point_t p(a.rows());
	_CLASS::DERIVED::multInto(a, b, p);
	return p;
}

_TEMPLATE
inline typename _CLASS::point_t _CLASS::exp(const point_t & point, const tangent_vector_t & vec) const
{
	point_t p(point.rows());
	this->getDerived().expInto(point, vec, p);
	return p;
}

_TEMPLATE
inline typename _CLASS::point_t _CLASS::expAtId(const tangent_vector_t & vec) const
{
	point_t result(this->getPointSize());
	this->getDerived().expAtIdInto(vec, result);
	return result;
}

_TEMPLATE
inline typename _CLASS::dmatrix_t _CLASS::dexpAtId(const tangent_vector_t & vec) const
{
	dmatrix_t result((int)this->getPointSize(), (int)this->getDimension());
	this->getDerived().dexpAtIdInto(vec, result);
	return result;
}

_TEMPLATE
inline typename _CLASS::dmatrix_t _CLASS::dexp(const point_t & point, const tangent_vector_t & vec) const
{
	dmatrix_t result((int)this->getPointSize(), (int)this->getDimension());
	this->getDerived().dexpInto(point, vec, result);
	return result;
}

_TEMPLATE
inline typename _CLASS::tangent_vector_t _CLASS::log(const point_t & from, const point_t & to) const
{
	tangent_vector_t v(this->getDimension());
	this->getDerived().logInto(from, to, v);
	return v;
}

_TEMPLATE
inline typename _CLASS::tangent_vector_t _CLASS::logAtId(const point_t & to) const
{
	tangent_vector_t v(this->getDimension());
	this->getDerived().logAtIdInto(to, v);
	return v;
}

_TEMPLATE
inline typename _CLASS::dlog_matrix_t _CLASS::dlogAtId(const point_t & to) const
{
	dlog_matrix_t result(this->getDimension(), this->getDimension());
	this->getDerived().dlogAtIdInto(to, result);
	return result;
}

_TEMPLATE
inline typename _CLASS::dlog_matrix_t _CLASS::dlog(const point_t & from, const point_t & to) const
{
	dlog_matrix_t result(this->getDimension(), this->getDimension());
	this->getDerived().dlogInto(from, to, result);
	return result;
}

#undef _CLASS
#undef _TEMPLATE
}
