#include "DiffManifoldBSplineTests.hpp"

namespace bsplines {

TEST(UnitQuaternionBSplineTestSuite, testQuaternionBSplineCompilation)
{
	UnitQuaternionBSpline<splineOrder>::TYPE rbspline;
	const UnitQuaternionBSpline<splineOrder>::TYPE::point_t p = UnitQuaternionBSpline<splineOrder>::TYPE::point_t(1, 0, 0, 0);
	rbspline.initConstantUniformSpline(minTime, maxTime, numberOfSegments, p);
	sm::eigen::assertEqual(rbspline.getEvaluatorAt<0>(minTime).eval(), p, SM_SOURCE_FILE_POS);
}

TEST(UnitQuaternionBSplineTestSuite, differentEvalMethodsEvalTheSame)
{
	UQTestSpline spline;
	initMinimalSpline(spline);

	UQTestSpline::point_t p;
	for(int i = 0, n = knot_arithmetics::getNumControlVerticesRequired(2, splineOrder) ; i < n; i ++){
		spline.getManifold().randomizePoint(p);
		spline.addControlVertex(i, p);
	}

	UQTestSpline::full_jacobian_t jac1, jac2;

	for(int i = 0, n = 10; i< n; i ++){
		double t = (spline.getMaxTime() - spline.getMinTime()) / (n - 1) * i + spline.getMinTime();
		UQTestSpline::Evaluator<2> eval = spline.getEvaluatorAt<2>(t);

//		std::cout << std::endl << "t = "<< t << std::endl;
		sm::eigen::assertNear(eval.evalDRecursive(0), eval.evalGeneric(), 1E-9, SM_SOURCE_FILE_POS);
		sm::eigen::assertNear(eval.evalDRecursive(1), eval.evalD1Special(), 1E-9, SM_SOURCE_FILE_POS);

		//compare generic recursive computation of Jacobian with the optimized one.s
		eval.evalJacobianDRecursive(0, jac1);
		eval.evalJacobian(jac2);
		sm::eigen::assertNear(jac1, jac2, 1E-9, SM_SOURCE_FILE_POS);
	}
}

TEST(UnitQuaternionBSplineTestSuite, testDExp)
{
	DExpTester<UnitQuaternionManifoldConf<>::Manifold >::testFunc(100, 10);
}

TEST(UnitQuaternionBSplineTestSuite, evalRiD1)
{
	SplineEvalRiDTester<UQTestSpline, 1>::testFunc(10, 10);
}
TEST(UnitQuaternionBSplineTestSuite, evalRiD2)
{
	SplineEvalRiDTester<UQTestSpline, 2>::testFunc(10, 10);
}
TEST(UnitQuaternionBSplineTestSuite, evalRiD3)
{
	SplineEvalRiDTester<UQTestSpline, 3>::testFunc(10, 10);
}

TEST(UnitQuaternionBSplineTestSuite, evalD1)
{
	SplineEvalDTester<UQTestSpline, 1>::testFunc(10, 10);
}

TEST(UnitQuaternionBSplineTestSuite, evalD2)
{
	SplineEvalDTester<UQTestSpline, 2>::testFunc(10, 10);
}

TEST(UnitQuaternionBSplineTestSuite, evalD3)
{
	SplineEvalDTester<UQTestSpline, 3>::testFunc(10, 10);
}

// TODO implement generic angular velocity tests
TEST(UnitQuaternionBSplineTestSuite, evalAngularVelocityAndAcceleration)
{
	UQTestSpline spline;

	std::vector<UQTestSpline::point_t> points;
	std::vector<double> times;

	Eigen::Vector3d direction(0.5, 0.5, 0);

	const int numPoints = 4;

	Eigen::MatrixXd interpolationPointsE = Eigen::MatrixXd::Zero(4, numPoints);
	Eigen::VectorXd timesE = Eigen::VectorXd::Zero(numPoints);

	for(int i = -splineOrder; i < numPoints; i++){
		UQTestSpline::point_t p = spline.getManifold().expAtId(direction * i);
		spline.addControlVertex(i, p);
	}

	spline.init();

//	BSplineFitter<UQTestSpline> f;
//
//	const int numSegments = numPoints;
//	for(int i = 0; i < numPoints; i++){
//		times.push_back(i);
//		timesE(i, 0) = i;
//		UQTestSpline::point_t p = spline.getManifold().expAtId(direction * i);
//		points.push_back(p);
//		interpolationPointsE.col(i) = p;
//		spline.addControlVertex(i, p);
//	}

//	double lambda = 1E-5;

//	f.initUniformSplineDense(spline, times, points, numSegments, lambda);

//
//	for(UQTestSpline::SegmentIterator i = spline.getAbsoluteBegin(); i != spline.getAbsoluteEnd(); i++){
//		std::cout << "t=" << i.getTime() <<" -> " << i->getControlVertex() << std::endl; // XXX: debug output of i
//		if(i.getTime() >= numPoints -1) break;
//	}

	Eigen::VectorXd accel = Eigen::VectorXd::Zero(3);

	const int numTimes = numPoints;
	for(int i = 0; i < numTimes; i ++){
		double t = spline.getMinTime() + (spline.getMaxTime() - spline.getMinTime()) * (double) i / (numTimes - 1);
		UQTestSpline::Evaluator<2> eval = spline.getEvaluatorAt<2>(t);

		sm::eigen::assertNear(eval.evalAngularVelocity() , direction, 1E-9, SM_SOURCE_FILE_POS);
		sm::eigen::assertNear(eval.evalAngularAcceleration() , accel, 1E-9, SM_SOURCE_FILE_POS);
	}
}


TEST(UnitQuaternionBSplineTestSuite, SplineEvalRiDJacobianTesterD0)
{
	SplineEvalRiDJacobianTester<UQTestSpline, 0>::testFunc(10, 1);
}
TEST(UnitQuaternionBSplineTestSuite, SplineEvalRiDJacobianTesterD1)
{
	SplineEvalRiDJacobianTester<UQTestSpline, 1>::testFunc(10, 1);
}
TEST(UnitQuaternionBSplineTestSuite, SplineEvalRiDJacobianTesterD2)
{
	SplineEvalRiDJacobianTester<UQTestSpline, 2>::testFunc(10, 1);
}

// Check the Jacobian calculation.
TEST(UnitQuaternionBSplineTestSuite, testBSplineJacobianD0)
{
	BSplineJacobianTester<UQTestSpline, 0>::testFunc(10, 1);
}
TEST(UnitQuaternionBSplineTestSuite, testBSplineJacobianD1)
{
	BSplineJacobianTester<UQTestSpline, 1>::testFunc(10, 1);
}
TEST(UnitQuaternionBSplineTestSuite, testBSplineJacobianD2)
{
	BSplineJacobianTester<UQTestSpline, 2>::testFunc(10, 1);
}

TEST(UnitQuaternionBSplineTestSuite, testAngularVelocitiyJacobian)
{
	AngularDerivativeJacobianTestser<UQTestSpline, 1>::testFunc(10, 1);
}
TEST(UnitQuaternionBSplineTestSuite, testAngularAccelerationJacobian)
{
	AngularDerivativeJacobianTestser<UQTestSpline, 2>::testFunc(10, 1);
}

} // namespace bsplines
