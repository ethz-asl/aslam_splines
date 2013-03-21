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

namespace bsplines {

template <typename TDiffManifoldConfiguration, int ISplineOrder, typename TTimePolicy, typename TModifiedDerivedConf>
class DiffManifoldBSpline<aslam::splines::DesignVariableSegmentBSplineConf<UnitQuaternionBSplineConfiguration<TDiffManifoldConfiguration, ISplineOrder, TTimePolicy>, TModifiedDerivedConf>, aslam::splines::DesignVariableSegmentBSplineConf<TModifiedDerivedConf> >
	: public DiffManifoldBSpline<aslam::splines::DesignVariableSegmentBSplineConf<typename UnitQuaternionBSplineConfiguration<TDiffManifoldConfiguration, ISplineOrder, TTimePolicy>::ParentConf, TModifiedDerivedConf>, aslam::splines::DesignVariableSegmentBSplineConf<TModifiedDerivedConf> >{
protected:
	typedef DiffManifoldBSpline<aslam::splines::DesignVariableSegmentBSplineConf<typename UnitQuaternionBSplineConfiguration<TDiffManifoldConfiguration, ISplineOrder, TTimePolicy>::ParentConf, TModifiedDerivedConf>, aslam::splines::DesignVariableSegmentBSplineConf<TModifiedDerivedConf> > parent_t;
	typedef aslam::splines::DesignVariableSegmentBSplineConf<TModifiedDerivedConf> CONF;
private:
	typedef typename ::bsplines::UnitQuaternionBSpline<ISplineOrder>::TYPE TBSpline;
public:
	typedef typename parent_t::spline_t spline_t;
	typedef typename parent_t::TimePolicy TimePolicy;
	typedef typename parent_t::time_t time_t;
	typedef typename parent_t::point_t point_t;
	typedef typename parent_t::SegmentIterator SegmentIterator;
	typedef typename parent_t::SegmentConstIterator SegmentConstIterator;
	typedef typename parent_t::angular_jacobian_t angular_jacobian_t;

	enum {
		PointSize = TBSpline::PointSize,
		Dimension = TBSpline::Dimension
	};

	typedef aslam::backend::VectorExpression<3> angular_derivative_expression_t;

	DiffManifoldBSpline(const CONF & config = CONF()) : parent_t(config){}

	template <int IMaxDerivativeOrder>
	class ExpressionFactory : public parent_t::template ExpressionFactory<IMaxDerivativeOrder> {
	public:
		typedef typename parent_t::template ExpressionFactory<IMaxDerivativeOrder> ParentExpressionFactory;
		typedef typename ParentExpressionFactory::eval_ptr_t  eval_ptr_t;
		ExpressionFactory(const spline_t & spline, const time_t & t) : ParentExpressionFactory(spline, t){}

		/// \brief get an expression
		angular_derivative_expression_t getAngularVelocityExpression() const;
		angular_derivative_expression_t getAngularAccelerationExpression() const;
	};

	template <int IMaxDerivativeOrder> inline ExpressionFactory<IMaxDerivativeOrder> getExpressionFactoryAt(const time_t & t) const { return ExpressionFactory<IMaxDerivativeOrder>(*this, t); }
};

} // namespace bsplines

#include "implementation/OPTUnitQuaternionBSplineImpl.hpp"

#endif /* OPTUNITQUATERNIONBSPLINE_HPP_ */
