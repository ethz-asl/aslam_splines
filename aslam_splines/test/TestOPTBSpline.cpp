/*
 * TestOPTBSpline.cpp
 *
 *  Created on: 05.08.2012
 *      Author: hannes
 */

#include <sm/eigen/gtest.hpp>
#include <sm/eigen/NumericalDiff.hpp>
#include <sm/kinematics/Transformation.hpp>
#include <aslam/splines/OPTBSpline.hpp>
#include <aslam/splines/OPTUnitQuaternionBSpline.hpp>
#include <bsplines/EuclideanBSpline.hpp>
#include <bsplines/UnitQuaternionBSpline.hpp>
#include <aslam/backend/TransformationBasic.hpp>
#include <aslam/backend/EuclideanPoint.hpp>
#include <aslam/backend/TransformationExpression.hpp>
#include <aslam/backend/Vector2RotationQuaternionExpressionAdapter.hpp>
#include <aslam/backend/ErrorTermTransformation.hpp>
#include <aslam/backend/test/ErrorTermTestHarness.hpp>
#include <aslam/backend/test/ExpressionTests.hpp>
#include <aslam/backend/test/RotationExpressionTests.hpp>

using namespace aslam::backend;
using namespace sm::kinematics;
using namespace aslam::splines;
using namespace bsplines;
using namespace std;

const int numSegments = 3, numberOfTimesToProbe = 10;
const double tmin = 0;
const double tmax = numSegments / 2.0;
const bool randomizeControlVertices = true;

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
	typedef typename OPTBSpline<typename TSplineMap::CONF>::BSpline TestBSpline;
	static void test(TestBSpline & spline, double t){}
};


template <typename TSplineMap, int ISplineOrder>
struct MaxDerivative {
	enum { VALUE = ISplineOrder };
};


template <typename TSplineMap, int ISplineOrder, int IDim>
struct OPTSplineTester{
	typedef typename OPTBSpline<typename TSplineMap::CONF>::BSpline TestBSpline;
	typedef typename TestBSpline::point_t point_t;
	typedef typename TestBSpline::tangent_vector_t tangent_vector_t;


	void static testCompilationAndExpressions(){

		TestBSpline bspline(createConf<typename TSplineMap::CONF, ISplineOrder, IDim>()), bsplineCopyAssign(createConf<typename TSplineMap::CONF, ISplineOrder, IDim>()), bsplineCopyConstruct(bspline);
		bsplineCopyAssign = bspline;

		const int pointSize = bspline.getPointSize();

		SM_ASSERT_EQ(std::runtime_error, ISplineOrder, bspline.getSplineOrder(), "");
		SM_ASSERT_EQ(std::runtime_error, IDim, bspline.getDimension(), "");

		bspline.initConstantUniformSpline(tmin, tmax, numSegments, bspline.getManifold().getDefaultPoint());
		double update[IDim];
		update[0] = 1;
		for (typename TestBSpline::SegmentIterator i = bspline.begin(), end = bspline.end(); i != end; i++) {
			point_t & p = i->getControlVertex();
			if(randomizeControlVertices){
				bspline.getManifold().randomizePoint(p);
			}
			point_t op = p;

			typename TestBSpline::dv_t & dv = *bspline.getDesignVariables(i.getKnot())[ISplineOrder - 1];
			SM_ASSERT_EQ(std::runtime_error, &dv, &i->getDesignVariable(), "");

			dv.update(update, IDim);

			point_t opu = op;
			::manifolds::internal::DiffManifoldPointUpdateTraits<typename TSplineMap::CONF::ManifoldConf>::update(opu, Eigen::Map<const tangent_vector_t>(update, IDim));

			sm::eigen::assertEqual(opu, p, SM_SOURCE_FILE_POS);

			dv.revertUpdate();
			sm::eigen::assertEqual(op, p, SM_SOURCE_FILE_POS);
		}


		for(int k = 0; k < numberOfTimesToProbe; k++){
			double t = tmin + ((k % 2 == 0) ? ((double) (k / 2) / ((numberOfTimesToProbe - 1)/ 2)) : ((double) rand() / RAND_MAX)) * (tmax - tmin);
			const int maxDerivativeOrder = MaxDerivative<TSplineMap, ISplineOrder>::VALUE;
			typename TestBSpline::template ExpressionFactory<maxDerivativeOrder> fact = bspline.template getExpressionFactoryAt <maxDerivativeOrder> (t);
			typename TestBSpline::template Evaluator<maxDerivativeOrder> eval = bspline.template getEvaluatorAt <maxDerivativeOrder > (t);
			for(int derivativeOrder = 0; derivativeOrder <= maxDerivativeOrder; derivativeOrder++){
				BOOST_AUTO(expression, fact.getValueExpression(derivativeOrder));

				sm::eigen::assertEqual(fact.getEvaluator().evalD(derivativeOrder), expression.evaluate(), SM_SOURCE_FILE_POS);
				sm::eigen::assertEqual(eval.evalD(derivativeOrder), expression.evaluate(), SM_SOURCE_FILE_POS);

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
				eval.evalJacobian(derivativeOrder, J);
				sm::eigen::assertEqual(J, jac.asDenseMatrix(), SM_SOURCE_FILE_POS);
				fact.getEvaluator().evalJacobian(derivativeOrder, J);
				sm::eigen::assertEqual(J, jac.asDenseMatrix(), SM_SOURCE_FILE_POS);

				{
					SCOPED_TRACE("");
					testJacobian(expression, false, 1E-3, 1E-6);
				}

				OPTSplineSpecializationTester<TSplineMap, ISplineOrder, IDim>::test(bspline, t);

				for(std::vector<DesignVariable *>::iterator i = varVec.begin(), end = varVec.end(); i != end; i++){
					(*i)->setActive(false);
					(*i)->setBlockIndex(-1);
				}
			}
		}
	}
};

template <int IEigenSplineOrder, int ISplineOrder>
struct MaxDerivative<UnitQuaternionBSpline<IEigenSplineOrder>, ISplineOrder> {
	enum { LIMIT=UnitQuaternionBSpline<IEigenSplineOrder>::TYPE::MaxSupportedDerivativeOrderJacobian, VALUE = LIMIT < ISplineOrder ? LIMIT : ISplineOrder };
};

template <int IEigenSplineOrder, int ISplineOrder, int IDim>
struct OPTSplineSpecializationTester<UnitQuaternionBSpline<IEigenSplineOrder>, ISplineOrder, IDim>
{
	typedef typename OPTBSpline<typename UnitQuaternionBSpline<IEigenSplineOrder>::CONF>::BSpline TestBSpline;

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

		{
			SCOPED_TRACE("");
			testJacobian(avexpression, false, 1E-3, 1E-6);
		}
		{
			SCOPED_TRACE("");
			testJacobian(aaexpression, false, 1E-3, 1E-6);
		}
	}
};

TEST(OPTBSplineTestSuite, testCompilationAndExpressionsEuclidean)
{
	OPTSplineTester<EuclideanBSpline<2, Eigen::Dynamic>, 2, 3>::testCompilationAndExpressions();
	OPTSplineTester<EuclideanBSpline<3, 2>, 3, 2>::testCompilationAndExpressions();
	OPTSplineTester<EuclideanBSpline<Eigen::Dynamic, 2>, 3, 2>::testCompilationAndExpressions();
	OPTSplineTester<EuclideanBSpline<Eigen::Dynamic, 2>, 4, 2>::testCompilationAndExpressions();
	OPTSplineTester<EuclideanBSpline<Eigen::Dynamic, 2>, 5, 2>::testCompilationAndExpressions();
	OPTSplineTester<EuclideanBSpline<Eigen::Dynamic, 2>, 6, 2>::testCompilationAndExpressions();
}

TEST(OPTBSplineTestSuite, testCompilationAndExpressionsUnitQuaternionStaticOrder2)
{
	OPTSplineTester<UnitQuaternionBSpline<2>, 2, 3>::testCompilationAndExpressions();
}

TEST(OPTBSplineTestSuite, testCompilationAndExpressionsUnitQuaternionDynamicOrder2)
{
	OPTSplineTester<UnitQuaternionBSpline<Eigen::Dynamic>, 2, 3>::testCompilationAndExpressions();
}

TEST(OPTBSplineTestSuite, testCompilationAndExpressionsUnitQuaternionDynamicOrder3)
{
	OPTSplineTester<UnitQuaternionBSpline<Eigen::Dynamic>, 3, 3>::testCompilationAndExpressions();
}

TEST(OPTBSplineTestSuite, testCompilationAndExpressionsUnitQuaternionDynamicOrder4)
{
	OPTSplineTester<UnitQuaternionBSpline<Eigen::Dynamic>, 4, 3>::testCompilationAndExpressions();
}

TEST(OPTBSplineTestSuite, testCompilationAndExpressionsUnitQuaternionDynamicOrder5)
{
	OPTSplineTester<UnitQuaternionBSpline<Eigen::Dynamic>, 5, 3>::testCompilationAndExpressions();
}

TEST(OPTBSplineTestSuite, testPoseErrorWithOPTSplines)
{
	try {
		using namespace aslam::backend;
		sm::kinematics::Transformation T_random;
		T_random.setRandom(0.05, 0.01);
		sm::kinematics::Transformation T_prior;

		OPTBSpline<typename UnitQuaternionBSpline<2>::CONF>::BSpline testSpline;

		testSpline.initConstantUniformSpline(0, 10, 2, quatRandom());

		double t = 5;
		int blockIndex = 0;
		for(DesignVariable * pDV : testSpline.getDesignVariables()){
			pDV->setActive(true);
			pDV->setBlockIndex(blockIndex++);
		}
		auto valueExpression = testSpline.getExpressionFactoryAt<0>(t).getValueExpression();
		{
			SCOPED_TRACE("");
			testJacobian(valueExpression);
		}

		RotationExpression rexp(Vector2RotationQuaternionExpressionAdapter::adapt(VectorExpression<4>(valueExpression)));
		for(int i = 0 ; i < 3 ; i++)
		{
			SCOPED_TRACE("");
			Eigen::Vector3d v;
			v.setZero();
			v(i) = 1;
			testJacobian(rexp * EuclideanExpression(v));
		}

		{
			SCOPED_TRACE("");
			testJacobian(rexp);
		}

		EuclideanPoint ep(T_prior.t());
		ep.setActive(true);
		ep.setBlockIndex(blockIndex++);

		EuclideanExpression eexp(&ep);


		TransformationBasic Tb(rexp, eexp);
		TransformationExpression T(&Tb);

		HomogeneousExpression he(Eigen::Vector4d::Random().eval());
		{
			SCOPED_TRACE("");
			testJacobian(T * he);
		}

		Eigen::MatrixXd N = Eigen::MatrixXd::Zero(6,6);
		N(0,0) = 1e-3;
		N(1,1) = 1e-3;
		N(2,2) = 1e-3;
		N(3,3) = 1e0;
		N(4,4) = 1e0;
		N(5,5) = 1e0;

		// Create the ErrorTerm
		ErrorTermTransformation ett(T, T_random, N);
		// Create the test harness
		aslam::backend::ErrorTermTestHarness<6> harness(&ett);

		// Run the unit tests.
		{
			SCOPED_TRACE("");
			harness.testAll(1e-5);
		}
	}
	catch(const std::exception & e)
	{
		FAIL() << e.what();
	}
}

