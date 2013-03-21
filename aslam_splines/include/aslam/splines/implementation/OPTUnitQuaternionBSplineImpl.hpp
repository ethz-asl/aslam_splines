/*
 * OPTBSplineImpl.hpp
 *
 *  Created on: 05.08.2012
 *      Author: hannes
 */

#ifndef OPTUNITQUATERNIONBSPLINEIMPL_HPP_
#define OPTUNITQUATERNIONBSPLINEIMPL_HPP_

#include "aslam/backend/JacobianContainer.hpp"
#include "boost/typeof/typeof.hpp"

using namespace aslam::backend;

namespace bsplines {

#define _TEMPLATE template <typename TDiffManifoldConfiguration, int ISplineOrder, typename TTimePolicy, typename TModifiedDerivedConf>
#define _CLASS DiffManifoldBSpline<aslam::splines::DesignVariableSegmentBSplineConf<UnitQuaternionBSplineConfiguration<TDiffManifoldConfiguration, ISplineOrder, TTimePolicy>, TModifiedDerivedConf>, aslam::splines::DesignVariableSegmentBSplineConf<TModifiedDerivedConf> >

_TEMPLATE
template <int IMaxDerivativeOrder>
typename _CLASS::angular_derivative_expression_t _CLASS::ExpressionFactory<IMaxDerivativeOrder>::getAngularVelocityExpression() const {
	typedef aslam::backend::VectorExpressionNode<3> node_t;

	class ExpressionNode : public node_t {
	public:
		ExpressionNode(const eval_ptr_t & evalPtr) : _evalPtr(evalPtr){}
	private:
		typedef typename node_t::vector_t vector_t;
		const eval_ptr_t _evalPtr;
		inline void evaluateJacobiansImplementation(JacobianContainer & outJacobians, const Eigen::MatrixXd * applyChainRule) const {
			const int dimension=_evalPtr->getSpline().getDimension(), pointSize = dimension, splineOrder = _evalPtr->getSpline().getSplineOrder();
			typename _CLASS::angular_jacobian_t J(pointSize, dimension * splineOrder);
			_evalPtr->evalAngularVelocityJacobian(J);

			int col = 0;
			for(SegmentConstIterator i = _evalPtr->begin(), end = _evalPtr->end(); i != end; ++i)
			{
				internal::AddJac<typename _CLASS::angular_jacobian_t, Dimension, Dimension>::addJac(J, col, &i->getDesignVariable(), outJacobians, applyChainRule, pointSize, dimension);
				col+=dimension;
			}
		}

	protected:
		virtual vector_t evaluateImplementation() const {
			return _evalPtr->evalAngularVelocity();
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
				designVariables.insert(const_cast<DesignVariable *>(&i->getDesignVariable()));
			}
		}
	};
	return angular_derivative_expression_t(boost::shared_ptr<node_t>(static_cast<node_t *> (new ExpressionNode(this->_evalPtr))));
}

_TEMPLATE
template <int IMaxDerivativeOrder>
typename _CLASS::angular_derivative_expression_t _CLASS::ExpressionFactory<IMaxDerivativeOrder>::getAngularAccelerationExpression() const {
	typedef aslam::backend::VectorExpressionNode<3> node_t;

	class ExpressionNode : public node_t {
	private:
		typedef typename node_t::vector_t vector_t;
		eval_ptr_t _evalPtr;
	public:
		ExpressionNode(const eval_ptr_t & evalPtr) : _evalPtr(evalPtr){}
		virtual ~ExpressionNode(){};
	private:
		inline void evaluateJacobiansImplementation(JacobianContainer & outJacobians, const Eigen::MatrixXd * applyChainRule) const {
			const int dimension=_evalPtr->getSpline().getDimension(), pointSize = dimension, splineOrder = _evalPtr->getSpline().getSplineOrder();
			typename _CLASS::angular_jacobian_t J(pointSize, dimension * splineOrder);
			_evalPtr->evalAngularAccelerationJacobian(J);

			int col = 0;
			for(SegmentConstIterator i = _evalPtr->begin(), end = _evalPtr->end(); i != end; ++i)
			{
				internal::AddJac<typename _CLASS::angular_jacobian_t, Dimension, Dimension>::addJac(J, col, &i->getDesignVariable(), outJacobians, applyChainRule, pointSize, dimension);
				col+=dimension;
			}
		}

	protected:
		virtual vector_t evaluateImplementation() const {
			return _evalPtr->evalAngularAcceleration();
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
				designVariables.insert(const_cast<DesignVariable *>(&i->getDesignVariable()));
			}
		}
	};
	return angular_derivative_expression_t(boost::shared_ptr<node_t>(static_cast<node_t *> (new ExpressionNode(this->_evalPtr))));
}

#undef _CLASS
#undef _TEMPLATE

} // namespace bsplines

#endif /* OPTUNITQUATERNIONBSPLINEIMPL_HPP_ */
