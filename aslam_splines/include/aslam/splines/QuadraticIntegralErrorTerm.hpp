/*
 * IntegrationErrorTerm.hpp
 *
 *  Created on: Jan 23, 2013
 *      Author: hannes
 */

#ifndef QINTEGRALERRORTERM_HPP_
#define QINTEGRALERRORTERM_HPP_


#include <iostream>
#include <boost/shared_ptr.hpp>
#include <boost/ptr_container/ptr_vector.hpp>
#include <sm/kinematics/rotations.hpp>
#include <sm/random.hpp>
#include <aslam/backend/OptimizationProblem.hpp>
#include <aslam/backend/Optimizer.hpp>
#include <aslam/backend/ErrorTerm.hpp>
#include <aslam/backend/ExpressionErrorTerm.hpp>

#include <aslam/splines/OPTBSpline.hpp>
#include <bsplines/BSplinePose.hpp>
#include <bsplines/UnitQuaternionBSpline.hpp>
#include <bsplines/NumericIntegrator.hpp>

namespace aslam {
namespace splines {
namespace integration {

namespace algorithms = ::numeric_integrator::algorithms;
typedef algorithms::Default DefaultAlgorithm;

namespace internal {
	inline void addErrorTermToProblem(OptimizationProblem & problem, ErrorTerm *et, bool problemOwnsErrorTerms){
		problem.addErrorTerm(et, problemOwnsErrorTerms);
	}
	inline void addErrorTermToProblem(OptimizationProblem & problem, boost::shared_ptr<ErrorTerm> et, bool problemOwnsErrorTerms){
		problem.addErrorTerm(et);
	}
}


/**
To add an error term of the form

	e = Int E(t)^T W(t) E(t) dt

to a problem use this function when E(t) is implemented as an ErrorTerm. If f is an expression use the next function.

The ErrorTermFactory functor must have an 

	ErrorTerm operator()(TTime t, double factor) const

member that returns E(t) with its invR = W * factor.
So don't forget to take the square root of factor when you use the setSqrtInvR method.
 */
template <typename Algorithm = DefaultAlgorithm, typename TTime, typename ErrorTermFactory>
void addQuadraticIntegralErrorTerms(OptimizationProblem & problem, const TTime & a, const TTime & b, int numberOfPoints, const ErrorTermFactory & errorTermFactory, bool problemOwnsErrorTerms = true)
{
	if(a == b) return;

	auto integrator = Algorithm().template getIntegrator<double>(a, b, numberOfPoints);
	SM_ASSERT_TRUE_DBG(std::runtime_error, !integrator.isAtEnd(), "too few integration points given : " << numberOfPoints);

	const double commonFactor = integrator.getCommonFactor();
	for(; !integrator.isAtEnd(); integrator.next()){
		double valueFactor = integrator.getValueFactor();
		internal::addErrorTermToProblem(problem, errorTermFactory(integrator.getIntegrationScalar(), commonFactor * valueFactor), problemOwnsErrorTerms);
	}
}

/**
To add an error term of the form

	e = \Int E(t)^T sqrtInvR^T sqrtInvR E(t) dt

to a problem use this function when E(t) is implemented as an VectorExpression.

The ExpressionFactory functor must have an

	VectorExpression operator()(TTime t) const
or
	GenericMatrixExpression<n, 1> operator()(TTime t) const

member that returns E(t).
 */
template <typename Algorithm = DefaultAlgorithm, typename TTime, typename ExpressionFactory, typename DerivedMatrix>
void addQuadraticIntegralExpressionErrorTerms(OptimizationProblem & problem, const TTime & a, const TTime & b, int numberOfPoints, const ExpressionFactory & expressionFactory, const Eigen::MatrixBase<DerivedMatrix> & sqrtInvR)
{
	addQuadraticIntegralErrorTerms<Algorithm>(problem, a, b, numberOfPoints, [&expressionFactory, &sqrtInvR](TTime t, double f){ return toErrorTermSqrt(expressionFactory(t), sqrtInvR * sqrt(f)); }, true);
}

}
}
}

#endif /* QINTEGRALERRORTERM_HPP_ */
