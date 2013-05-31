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

using namespace aslam::backend;

namespace bsplines {

#define _TEMPLATE template<typename TModifiedConf, typename TModifiedDerivedConf>
#define _CLASS DiffManifoldBSpline<aslam::splines::DesignVariableSegmentBSplineConf<TModifiedConf, TModifiedDerivedConf>, aslam::splines::DesignVariableSegmentBSplineConf<TModifiedDerivedConf> >

//_TEMPLATE
//_CLASS::~OPTBSpline() {
//
//}

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
template <int IMaxDerivativeOrder>
_CLASS::ExpressionFactory<IMaxDerivativeOrder>::ExpressionFactory(const spline_t & spline, const time_t & t) : _evalPtr(new eval_t(spline, t)) {
}

namespace internal {
	template <typename TDerived>
	inline void addJacImpl(const DesignVariable * designVariable, JacobianContainer & outJacobians, const Eigen::MatrixXd * applyChainRule, const Eigen::MatrixBase<TDerived> & block){
		if(applyChainRule){
			//TODO discuss: why not const DesignVariable * in add ?
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
	struct AddJac<TJacobian, IPointSize, IDimension, true>  {
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
void _CLASS::addSegment(time_t t, const point_t & p) {
//TODO	_bspline.addCurveSegment(t, p);
//	_designVariables.push_back(
//			new _CLASS::dv_t(
//					_bspline.fixedSizeVvCoefficientVector<D>(
//							_bspline.numVvCoefficients()
//									- 1)));
}

_TEMPLATE
void _CLASS::addSegment2(time_t t, const point_t & p, double lambda) {
//TODO	_bspline.addCurveSegment2(t, p, lambda);
//	_designVariables.push_back(
//			new _CLASS::dv_t(
//					_bspline.fixedSizeVvCoefficientVector<D>(
//							_bspline.numVvCoefficients()
//									- 1)));
}

_TEMPLATE
void _CLASS::removeSegment() {
//TODO	_bspline.removeCurveSegment();
//	_designVariables.erase(_designVariables.begin());
}

#undef _CLASS
#undef _TEMPLATE

} // namespace bsplines

#endif /* OPTBSPLINEIMPL_HPP_ */
