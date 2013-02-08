/*
 * TestOPTBSpline.cpp
 *
 *  Created on: 05.08.2012
 *      Author: hannes
 */


#include <sm/eigen/gtest.hpp>
#include <sm/eigen/NumericalDiff.hpp>
#include <aslam/splines/OPTBSpline.hpp>
#include <aslam/splines/OPTUnitQuaternionBSpline.hpp>
#include <bsplines/EuclideanBSpline.hpp>
#include <bsplines/UnitQuaternionBSpline.hpp>

using namespace aslam::backend;
using namespace sm::kinematics;
using namespace aslam::splines;
using namespace bsplines;
using namespace std;

const int numSegments = 10, numberOfTimesToProbe = 10;
const double tmin = 0;
const double tmax = numSegments;

template <typename TConf, int ISplineOrder, int IDim, bool BDimRequired> struct ConfCreator {
	static inline TConf create(){
		return TConf(typename TConf::ManifoldConf(IDim), ISplineOrder);
	}
};

template <typename TConf, int ISplineOrder, int IDim> struct ConfCreator<TConf, ISplineOrder, IDim, false> {
	static inline TConf create(){
		BOOST_STATIC_ASSERT_MSG(IDim == TConf::Dimension::VALUE, "impossible dimension selected!");
		return TConf(typename TConf::ManifoldConf(), ISplineOrder);
	}
};

template <typename TConf, int ISplineOrder, int IDim> inline TConf createConf(){
	return ConfCreator<TConf, ISplineOrder, IDim, TConf::Dimension::IS_DYNAMIC>::create();
}


template <typename TSplineMap, int ISplineOrder, int IDim>
struct OPTSplineSpecializationTester
{
	typedef OPTBSpline<typename TSplineMap::CONF> TestBSpline;
	static void test(TestBSpline & spline, double t){}
};


template <typename TSplineMap, int ISplineOrder, int IDim>
struct OPTSplineTester{
	typedef OPTBSpline<typename TSplineMap::CONF> TestBSpline;
	typedef typename TestBSpline::point_t point_t;
	typedef typename TestBSpline::tangent_vector_t tangent_vector_t;


	void static testCompilationAndExpressions(){

		TestBSpline bspline(createConf<typename TSplineMap::CONF, ISplineOrder, IDim>());

		const int pointSize = bspline.getPointSize();

		typename TestBSpline::point_t initPoint(pointSize);

		bspline.getManifold().randomizePoint(initPoint);
		bspline.initConstantUniformSpline(tmin, tmax, numSegments, initPoint);
		double update[IDim];
		update[0] = 1;
		for (typename TestBSpline::SegmentIterator i = bspline.begin(), end =
				bspline.end(); i != end; i++) {
			point_t & p = i->getControlVertex();
			bspline.getManifold().randomizePoint(p);
			point_t op = p;

			typename TestBSpline::dv_t & dv = *bspline.getDesignVariables(i.getKnot())[ISplineOrder - 1];
			SM_ASSERT_EQ(std::runtime_error, &dv, &i->getDesignVariable(),
					"");

			dv.update(update, IDim);

			point_t opu = op;
			::manifolds::internal::DiffManifoldPointUpdateTraits<typename TSplineMap::CONF::ManifoldConf>::update(opu, Eigen::Map<const tangent_vector_t>(update, IDim));

			sm::eigen::assertEqual(opu, p, SM_SOURCE_FILE_POS);

			dv.revertUpdate();
			sm::eigen::assertEqual(op, p, SM_SOURCE_FILE_POS);
		}


		for(int k = 0; k < numberOfTimesToProbe; k++){
			double t = tmin + ((double) rand() / RAND_MAX) * (tmax - tmin);
			typename TestBSpline::template ExpressionFactory<1> fact = bspline.template getExpressionFactoryAt < 1 > (t);

			BOOST_AUTO(expression, fact.toExpression(1));

			sm::eigen::assertEqual(fact.getEvaluator().evalD(1), expression.evaluate(), SM_SOURCE_FILE_POS);


			typename TestBSpline::template Evaluator<1> eval = bspline.template getEvaluatorAt < 1 > (t);
			sm::eigen::assertEqual(eval.evalD(1), expression.evaluate(), SM_SOURCE_FILE_POS);

			JacobianContainer jac(pointSize);
			DesignVariable::set_t set;
			expression.getDesignVariables(set);

			std::vector<DesignVariable *> varVec = bspline.getDesignVariables(t);
			SM_ASSERT_EQ(std::runtime_error, set.size(), varVec.size(), "");

			int c = 0;
			for(std::vector<DesignVariable *>::iterator i = varVec.begin(), end = varVec.end(); i != end; i++){
				(*i)->setActive(true);
				(*i)->setBlockIndex(c++);
				SM_ASSERT_EQ(std::runtime_error, set.count(*i), 1, "");
			}

			expression.evaluateJacobians(jac);

			typename TestBSpline::full_jacobian_t J;
			eval.evalJacobian(1, J);
			sm::eigen::assertEqual(J, jac.asDenseMatrix(), SM_SOURCE_FILE_POS);

			OPTSplineSpecializationTester<TSplineMap, ISplineOrder, IDim>::test(bspline, t);

			for(std::vector<DesignVariable *>::iterator i = varVec.begin(), end = varVec.end(); i != end; i++){
				(*i)->setActive(false);
				(*i)->setBlockIndex(-1);
			}
		}
	}
};



template <int IEigenSplineOrder, int ISplineOrder, int IDim>
struct OPTSplineSpecializationTester<UnitQuaternionBSpline<IEigenSplineOrder>, ISplineOrder, IDim>
{
	typedef OPTBSpline<typename UnitQuaternionBSpline<IEigenSplineOrder>::CONF> TestBSpline;

	static void test(TestBSpline & bspline, double t){
		typename TestBSpline::template ExpressionFactory<2> fact = bspline.template getExpressionFactoryAt < 2 > (t);
		BOOST_AUTO(avexpression, fact.getAngularVelocityExpression());
		BOOST_AUTO(aaexpression, fact.getAngularAccelerationExpression());

		typename TestBSpline::template Evaluator<2> eval = bspline.template getEvaluatorAt < 2 > (t);
		sm::eigen::assertEqual(eval.evalAngularVelocity(), avexpression.evaluate(), SM_SOURCE_FILE_POS);
		sm::eigen::assertEqual(eval.evalAngularAcceleration(), aaexpression.evaluate(), SM_SOURCE_FILE_POS);

		JacobianContainer jac(3);
		typename TestBSpline::angular_jacobian_t J;

		avexpression.evaluateJacobians(jac);
		eval.evalAngularVelocityJacobian(J);
		sm::eigen::assertEqual(J, jac.asDenseMatrix(), SM_SOURCE_FILE_POS);

		jac.clear();
		aaexpression.evaluateJacobians(jac);
		eval.evalAngularAccelerationJacobian(J);
		sm::eigen::assertEqual(J, jac.asDenseMatrix(), SM_SOURCE_FILE_POS);
	}
};

TEST(OPTBSplineTestSuite, testCompilationAndExpressions)
{
	OPTSplineTester<EuclideanBSpline<2, Eigen::Dynamic>, 2, 3>::testCompilationAndExpressions();
	OPTSplineTester<EuclideanBSpline<3, 2>, 3, 2>::testCompilationAndExpressions();
	OPTSplineTester<EuclideanBSpline<Eigen::Dynamic, 2>, 3, 2>::testCompilationAndExpressions();
	OPTSplineTester<UnitQuaternionBSpline<2>, 2, 3>::testCompilationAndExpressions();
	OPTSplineTester<UnitQuaternionBSpline<Eigen::Dynamic>, 4, 3>::testCompilationAndExpressions();
}
