/*
 * OPTBSplineImpl.hpp
 *
 *  Created on: 05.08.2012
 *      Author: hannes
 */

#ifndef OPTBSPLINEIMPL_HPP_
#define OPTBSPLINEIMPL_HPP_

#include "aslam/backend/JacobianContainer.hpp"
#include "boost/typeof/typeof.hpp"
#include "bsplines/DiffManifoldBSpline.hpp"
#include "bsplines/manifolds/LieGroup.hpp"

using namespace aslam::backend;

namespace bsplines {

#define _TEMPLATE template<typename TModifiedConf, typename TModifiedDerivedConf>
#define _CLASS DiffManifoldBSpline<aslam::splines::DesignVariableSegmentBSplineConf<TModifiedConf, TModifiedDerivedConf>, aslam::splines::DesignVariableSegmentBSplineConf<TModifiedDerivedConf> >

_TEMPLATE
void _CLASS::init() {
	parent_t::init();
	// Create all of the design variables as maps into the vector of spline coefficients.

	_designVariables.clear();
	for(SegmentIterator i = this->firstRelevantSegment(), end = this->end(); i != end; i++)
	{
		_designVariables.push_back((dv_t* )& i->getDesignVariable());
	}
}

_TEMPLATE
size_t _CLASS::numDesignVariables() {
	return _designVariables.size();
}

_TEMPLATE
typename _CLASS::dv_t * _CLASS::designVariable(
		size_t i) {
	SM_ASSERT_LT(aslam::Exception, i,
			_designVariables.size(), "Index out of bounds");
	return _designVariables[i];
}

_TEMPLATE
const std::vector<typename _CLASS::dv_t *> & _CLASS::getDesignVariables() {
	return _designVariables;
}

_TEMPLATE
template <int IMaxDerivativeOrder>
_CLASS::ExpressionFactory<IMaxDerivativeOrder>::ExpressionFactory(const spline_t & spline, const time_t & t) : _evalPtr(new eval_t(spline, t)) {
}

namespace internal {
	template <typename TDiffManifoldBSplineConfiguration, typename TDiffManifoldBSplineConfigurationDerived>
	inline void SegmentData< ::aslam::splines::DesignVariableSegmentBSplineConf<TDiffManifoldBSplineConfiguration, TDiffManifoldBSplineConfigurationDerived> >::minimalDifferenceImplementation(const Eigen::MatrixXd& xHat, Eigen::VectorXd& outDifference) const {
		UpdateTraits::minimalDifference(_manifold, xHat, this->getControlVertex(), outDifference);
	}

	template <typename TDiffManifoldBSplineConfiguration, typename TDiffManifoldBSplineConfigurationDerived>
	inline void SegmentData< ::aslam::splines::DesignVariableSegmentBSplineConf<TDiffManifoldBSplineConfiguration, TDiffManifoldBSplineConfigurationDerived> >::minimalDifferenceAndJacobianImplementation(const Eigen::MatrixXd& xHat, Eigen::VectorXd& outDifference, Eigen::MatrixXd& outJacobian) const {
		UpdateTraits::minimalDifferenceAndJacobian(_manifold, xHat, this->getControlVertex(), outDifference, outJacobian);
	}

	template <typename TDerived>
	inline void addJacImpl(const DesignVariable * designVariable, JacobianContainer & outJacobians, const Eigen::MatrixXd * applyChainRule, const Eigen::MatrixBase<TDerived> & block){
		if(applyChainRule){
			outJacobians.add(const_cast<DesignVariable *>(designVariable), (*applyChainRule) * block);
		}
		else{
			outJacobians.add(const_cast<DesignVariable *>(designVariable), block);
		}
	}

	template <typename TJacobian, int IPointSize, int IDimension, bool Dynamic = TJacobian::SizeAtCompileTime == Eigen::Dynamic>
	struct AddJac {
		static inline void addJac(TJacobian & J, const int col, const DesignVariable * designVariable, JacobianContainer & outJacobians, const Eigen::MatrixXd * applyChainRule, int pointSize, int dimension){
			addJacImpl(designVariable, outJacobians, applyChainRule, J.template block<IPointSize,IDimension>(0,col));
		}
	};

	template <typename TJacobian, int IPointSize, int IDimension>
	struct AddJac<TJacobian, IPointSize, IDimension, true> {
		static inline void addJac(TJacobian & J, const int col, const DesignVariable * designVariable, JacobianContainer & outJacobians, const Eigen::MatrixXd * applyChainRule, int pointSize, int dimension){
			addJacImpl(designVariable, outJacobians, applyChainRule, J.block(0,col, pointSize, dimension));
		}
	};
}

_TEMPLATE
template <int IMaxDerivativeOrder>
typename _CLASS::expression_t _CLASS::ExpressionFactory<IMaxDerivativeOrder>::getValueExpression(const int derivativeOrder) const {
	typedef aslam::backend::VectorExpressionNode<PointSize> node_t;

	class ExpressionNode : public node_t {
		typedef typename node_t::vector_t vector_t;
		eval_ptr_t _evalPtr;
		const int _derivativeOrder;
	public:
		ExpressionNode(const eval_ptr_t & evalPtr, int derivativeOrder) : _evalPtr(evalPtr), _derivativeOrder(derivativeOrder){}
		virtual ~ExpressionNode(){}
	protected:
		inline void evaluateJacobiansImplementation(JacobianContainer & outJacobians, const Eigen::MatrixXd * applyChainRule) const {
			const int dimension=_evalPtr->getSpline().getDimension(), pointSize = _evalPtr->getSpline().getPointSize(), splineOrder = _evalPtr->getSpline().getSplineOrder();
			typename _CLASS::full_jacobian_t J(pointSize, dimension * splineOrder);
			_evalPtr->evalJacobian(_derivativeOrder, J);

			int col = 0;
			for(SegmentConstIterator i = _evalPtr->begin(), end = _evalPtr->end(); i != end; ++i)
			{
				internal::AddJac<typename _CLASS::full_jacobian_t, PointSize, Dimension>::addJac(J, col, &i->getDesignVariable(), outJacobians, applyChainRule, pointSize, dimension);
				col+=dimension;
			}
		}

		virtual vector_t evaluateImplementation() const {
			return _evalPtr->evalD(_derivativeOrder);
		}

		virtual void evaluateJacobiansImplementation(JacobianContainer & outJacobians) const {
			evaluateJacobiansImplementation(outJacobians, NULL);
		}

		virtual void evaluateJacobiansImplementation(JacobianContainer & outJacobians, const Eigen::MatrixXd & applyChainRule) const {
			evaluateJacobiansImplementation(outJacobians, &applyChainRule);
		}

		virtual void getDesignVariablesImplementation(DesignVariable::set_t & designVariables) const {
			for(SegmentConstIterator i = _evalPtr->begin(), end = _evalPtr->end(); i != end; ++i)
			{
		//TODO discuss why is this method declared const in the parent, while the set is not of const *DesignVariables? - requires the following const_cast?
				designVariables.insert(const_cast<DesignVariable *>(&i->getDesignVariable()));
			}
		}
	};
	return expression_t(boost::shared_ptr<node_t>(static_cast<node_t *> (new ExpressionNode(_evalPtr, derivativeOrder))));
}

_TEMPLATE
std::vector<typename _CLASS::dv_t *> _CLASS::getDesignVariables(time_t tk) {
	std::vector<dv_t *> dvs(this->getSplineOrder());
	using namespace std;

	int j = 0;
	for (SegmentIterator back = this->getSegmentIterator(tk), i = this->getFirstRelevantSegmentByLast(back), end=++back; i != end; ++i) {
		dvs[j++] = &i->getDesignVariable();
	}
	return dvs;
}


_TEMPLATE
typename _CLASS::time_t _CLASS::appendSegments(KnotGenerator<time_t> & knotGenerator, int numSegments, const point_t * value) {
	SegmentIterator i = this->end();

	time_t t = parent_t::appendSegments(knotGenerator, numSegments, value);
	for(SegmentIterator end = this->end(); i != end; i++)
	{
		_designVariables.push_back((dv_t* )& i->getDesignVariable());
	}
	return t;
}

_TEMPLATE
void _CLASS::removeSegment() {
//TODO	_bspline.removeCurveSegment();
//	_designVariables.erase(_designVariables.begin());
}

#undef _CLASS
#undef _TEMPLATE

} // namespace bsplines

#ifdef USE_EXTERN_TEMPLATES
#ifndef USE_EXTERN_TEMPLATES_INSTANCIATION_FILE
#define OPTBSPLINEIMPL_HPP_EXTERN_KEYWORD extern
#else
#define OPTBSPLINEIMPL_HPP_EXTERN_KEYWORD
#endif
#include <bsplines/EuclideanBSpline.hpp>
#include <bsplines/UnitQuaternionBSpline.hpp>
#include <aslam/splines/OPTUnitQuaternionBSpline.hpp>

#ifdef USE_EXTERN_TEMPLATES_INSTANCIATION_FILE
OPTBSPLINEIMPL_HPP_EXTERN_KEYWORD template class ::bsplines::DiffManifoldBSpline< ::bsplines::EuclideanBSplineConfiguration< ::manifolds::EuclideanSpaceConf<> > >;
OPTBSPLINEIMPL_HPP_EXTERN_KEYWORD template class ::aslam::splines::OPTBSpline< ::bsplines::EuclideanBSplineConfiguration< ::manifolds::EuclideanSpaceConf<> > >;

OPTBSPLINEIMPL_HPP_EXTERN_KEYWORD template class ::bsplines::DiffManifoldBSpline< ::bsplines::EuclideanBSplineConfiguration< ::manifolds::EuclideanSpaceConf<>, 1> >;
OPTBSPLINEIMPL_HPP_EXTERN_KEYWORD template class ::aslam::splines::OPTBSpline< ::bsplines::EuclideanBSplineConfiguration< ::manifolds::EuclideanSpaceConf<>, 1> >;

OPTBSPLINEIMPL_HPP_EXTERN_KEYWORD template class ::bsplines::DiffManifoldBSpline< ::bsplines::EuclideanBSplineConfiguration< ::manifolds::EuclideanSpaceConf<>, 2> >;
OPTBSPLINEIMPL_HPP_EXTERN_KEYWORD template class ::aslam::splines::OPTBSpline< ::bsplines::EuclideanBSplineConfiguration< ::manifolds::EuclideanSpaceConf<>, 2> >;


OPTBSPLINEIMPL_HPP_EXTERN_KEYWORD template class ::bsplines::DiffManifoldBSpline< ::bsplines::UnitQuaternionBSplineConfiguration< > >;
OPTBSPLINEIMPL_HPP_EXTERN_KEYWORD template class ::aslam::splines::DesignVariableSegmentBSplineConf< ::bsplines::UnitQuaternionBSplineConfiguration< > >;
OPTBSPLINEIMPL_HPP_EXTERN_KEYWORD template class ::aslam::splines::OPTBSpline< ::bsplines::UnitQuaternionBSplineConfiguration< > >;

OPTBSPLINEIMPL_HPP_EXTERN_KEYWORD template class ::bsplines::DiffManifoldBSpline< ::bsplines::UnitQuaternionBSplineConfiguration< ::manifolds::UnitQuaternionManifoldConf<>, 1> >;
OPTBSPLINEIMPL_HPP_EXTERN_KEYWORD template class ::aslam::splines::DesignVariableSegmentBSplineConf< ::bsplines::UnitQuaternionBSplineConfiguration< ::manifolds::UnitQuaternionManifoldConf<>, 1> >;
OPTBSPLINEIMPL_HPP_EXTERN_KEYWORD template class ::aslam::splines::OPTBSpline< ::bsplines::UnitQuaternionBSplineConfiguration< ::manifolds::UnitQuaternionManifoldConf<>, 1> >;

OPTBSPLINEIMPL_HPP_EXTERN_KEYWORD template class ::bsplines::DiffManifoldBSpline< ::bsplines::UnitQuaternionBSplineConfiguration< ::manifolds::UnitQuaternionManifoldConf<>, 2> >;
OPTBSPLINEIMPL_HPP_EXTERN_KEYWORD template class ::aslam::splines::DesignVariableSegmentBSplineConf< ::bsplines::UnitQuaternionBSplineConfiguration< ::manifolds::UnitQuaternionManifoldConf<>, 2> >;
OPTBSPLINEIMPL_HPP_EXTERN_KEYWORD template class ::aslam::splines::OPTBSpline< ::bsplines::UnitQuaternionBSplineConfiguration< ::manifolds::UnitQuaternionManifoldConf<>, 2> >;

OPTBSPLINEIMPL_HPP_EXTERN_KEYWORD template class ::bsplines::DiffManifoldBSpline< ::bsplines::UnitQuaternionBSplineConfiguration< ::manifolds::UnitQuaternionManifoldConf<>, 3> >;
OPTBSPLINEIMPL_HPP_EXTERN_KEYWORD template class ::aslam::splines::DesignVariableSegmentBSplineConf< ::bsplines::UnitQuaternionBSplineConfiguration< ::manifolds::UnitQuaternionManifoldConf<>, 3> >;
OPTBSPLINEIMPL_HPP_EXTERN_KEYWORD template class ::aslam::splines::OPTBSpline< ::bsplines::UnitQuaternionBSplineConfiguration< ::manifolds::UnitQuaternionManifoldConf<>, 3> >;

OPTBSPLINEIMPL_HPP_EXTERN_KEYWORD template class ::bsplines::DiffManifoldBSpline< ::bsplines::UnitQuaternionBSplineConfiguration< ::manifolds::UnitQuaternionManifoldConf<>, 4> >;
OPTBSPLINEIMPL_HPP_EXTERN_KEYWORD template class ::aslam::splines::DesignVariableSegmentBSplineConf< ::bsplines::UnitQuaternionBSplineConfiguration< ::manifolds::UnitQuaternionManifoldConf<>, 4> >;
OPTBSPLINEIMPL_HPP_EXTERN_KEYWORD template class ::aslam::splines::OPTBSpline< ::bsplines::UnitQuaternionBSplineConfiguration< ::manifolds::UnitQuaternionManifoldConf<>, 4> >;

OPTBSPLINEIMPL_HPP_EXTERN_KEYWORD template class ::aslam::splines::OPTBSpline<typename ::bsplines::EuclideanBSpline<3>::CONF>;
OPTBSPLINEIMPL_HPP_EXTERN_KEYWORD template class ::aslam::splines::OPTBSpline<typename ::bsplines::EuclideanBSpline<4>::CONF>;
OPTBSPLINEIMPL_HPP_EXTERN_KEYWORD template class ::aslam::splines::OPTBSpline<typename ::bsplines::EuclideanBSpline<Eigen::Dynamic, 1>::CONF>;
OPTBSPLINEIMPL_HPP_EXTERN_KEYWORD template class ::aslam::splines::OPTBSpline<typename ::bsplines::EuclideanBSpline<2, 1>::CONF>;
OPTBSPLINEIMPL_HPP_EXTERN_KEYWORD template class ::aslam::splines::OPTBSpline<typename ::bsplines::EuclideanBSpline<3, 1>::CONF>;
OPTBSPLINEIMPL_HPP_EXTERN_KEYWORD template class ::aslam::splines::OPTBSpline<typename ::bsplines::EuclideanBSpline<4, 1>::CONF>;
OPTBSPLINEIMPL_HPP_EXTERN_KEYWORD template class ::aslam::splines::OPTBSpline<typename ::bsplines::EuclideanBSpline<Eigen::Dynamic, 2>::CONF>;
OPTBSPLINEIMPL_HPP_EXTERN_KEYWORD template class ::aslam::splines::OPTBSpline<typename ::bsplines::EuclideanBSpline<2, 2>::CONF>;
OPTBSPLINEIMPL_HPP_EXTERN_KEYWORD template class ::aslam::splines::OPTBSpline<typename ::bsplines::EuclideanBSpline<3, 2>::CONF>;
OPTBSPLINEIMPL_HPP_EXTERN_KEYWORD template class ::aslam::splines::OPTBSpline<typename ::bsplines::EuclideanBSpline<4, 2>::CONF>;
OPTBSPLINEIMPL_HPP_EXTERN_KEYWORD template class ::aslam::splines::OPTBSpline<typename ::bsplines::EuclideanBSpline<Eigen::Dynamic, 3>::CONF>;
OPTBSPLINEIMPL_HPP_EXTERN_KEYWORD template class ::aslam::splines::OPTBSpline<typename ::bsplines::EuclideanBSpline<2, 3>::CONF>;
OPTBSPLINEIMPL_HPP_EXTERN_KEYWORD template class ::aslam::splines::OPTBSpline<typename ::bsplines::EuclideanBSpline<3, 3>::CONF>;
OPTBSPLINEIMPL_HPP_EXTERN_KEYWORD template class ::aslam::splines::OPTBSpline<typename ::bsplines::EuclideanBSpline<4, 3>::CONF>;
#endif

#undef OPTBSPLINEIMPL_HPP_EXTERN_KEYWORD

#endif

#endif /* OPTBSPLINEIMPL_HPP_ */
