/*
 * OPTBSpline.hpp
 *
 *  Created on: 05.08.2012
 *      Author: hannes
 */

#ifndef OPTUNITQUATERNIONBSPLINE_HPP_
#define OPTUNITQUATERNIONBSPLINE_HPP_

#include "OPTBSpline.hpp"
#include "bsplines/UnitQuaternionBSpline.hpp"

namespace aslam {
namespace splines {

template <typename TDiffManifoldConfiguration, int ISplineOrder, typename TTimePolicy, typename TDerivedConf>
class OPTBSpline< ::bsplines::UnitQuaternionBSplineConfiguration<TDiffManifoldConfiguration, ISplineOrder, TTimePolicy >, TDerivedConf > : public OPTBSpline< typename ::bsplines::UnitQuaternionBSplineConfiguration<TDiffManifoldConfiguration, ISplineOrder, TTimePolicy>::ParentConf, TDerivedConf > {
protected:
	typedef typename ::bsplines::UnitQuaternionBSpline<ISplineOrder>::TYPE TBSpline;
	typedef OPTBSpline< typename ::bsplines::UnitQuaternionBSplineConfiguration<TDiffManifoldConfiguration, ISplineOrder, TTimePolicy>::ParentConf, TDerivedConf > parent_t;
	typedef typename parent_t::configuration_t configuration_t;
	typedef typename parent_t::CONF CONF;
public:
	typedef typename parent_t::spline_t spline_t;
	typedef typename parent_t::TimePolicy TimePolicy;
	typedef typename parent_t::time_t time_t;
	typedef typename parent_t::point_t point_t;
	typedef typename parent_t::SegmentIterator SegmentIterator;
	typedef typename parent_t::SegmentConstIterator SegmentConstIterator;
	typedef typename CONF::BSpline::angular_jacobian_t angular_jacobian_t;

	enum {
		PointSize = TBSpline::PointSize,
		Dimension = TBSpline::Dimension
	};

	typedef aslam::backend::VectorExpression<3> angular_derivative_expression_t;

	OPTBSpline(const CONF & config) : parent_t(configuration_t(config)){}

	template <int IMaxDerivativeOrder>
	class ExpressionFactory : public parent_t::template ExpressionFactory<IMaxDerivativeOrder> {
	public:
		typedef typename parent_t::template ExpressionFactory<IMaxDerivativeOrder>::eval_ptr_t eval_ptr_t;
		ExpressionFactory(spline_t & spline, const time_t & t) : parent_t::template ExpressionFactory<IMaxDerivativeOrder>(spline, t){}

		/// \brief get an expression
		angular_derivative_expression_t getAngularVelocityExpression();
		angular_derivative_expression_t getAngularAccelerationExpression();
	};

	template <int IMaxDerivativeOrder> inline ExpressionFactory<IMaxDerivativeOrder> getExpressionFactoryAt(const time_t & t){ return ExpressionFactory<IMaxDerivativeOrder>(*this, t); }
};

} // namespace splines
} // namespace aslam

#include "implementation/OPTUnitQuaternionBSplineImpl.hpp"

#endif /* OPTUNITQUATERNIONBSPLINE_HPP_ */
