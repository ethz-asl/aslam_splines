#include "bsplines/manifolds/EuclideanSpace.hpp"
#include "bsplines/manifolds/UnitQuaternionManifold.hpp"
#include "gtest/gtest.h"
#include <sm/eigen/gtest.hpp>

#define TEST_SPLINES

#include "NodeDistributedCacheTests.cpp"
#include "NumericIntegratorTests.cpp"


namespace bsplines {

	template <int IDim, bool dynamic, typename TScalar >
	void testManifolds() {
		using namespace manifolds;
		typedef Eigen::Matrix<TScalar, IDim, 1> Vec;

		Vec rand = Vec::Random();

		typename EuclideanSpaceConf<dynamic? Eigen::Dynamic : IDim, TScalar>::Manifold m(IDim);
		sm::eigen::assertNear(m.getIdentity(), Vec::Zero(), 1E-9, SM_SOURCE_FILE_POS);
		sm::eigen::assertNear(m.exp(Vec::Ones(), Vec::Ones()), 2 * Vec::Ones(), 1E-9, SM_SOURCE_FILE_POS);
		sm::eigen::assertNear(m.exp(Vec::Zero(), rand), m.expAtId(rand), 1E-9, SM_SOURCE_FILE_POS);

		UnitQuaternionManifoldConf<>::Manifold u;
		Eigen::Vector4d uqId = Eigen::Vector4d::Zero();
		uqId[3] = 1;
		sm::eigen::assertNear(u.getIdentity(), uqId, 1E-9, SM_SOURCE_FILE_POS);
		sm::eigen::assertNear(u.exp(uqId, Eigen::Vector3d::Zero()), uqId, 1E-9, SM_SOURCE_FILE_POS);
		Eigen::Vector3d randDouble = Eigen::Vector3d::Random();
		sm::eigen::assertNear(u.exp(uqId, randDouble), u.expAtId(randDouble), 1E-9, SM_SOURCE_FILE_POS);
	}

	TEST(DiffManifoldBSplineTestSuite, testManifolds)
	{
		testManifolds<3, Eigen::Dynamic, double>();
		testManifolds<3, 3, double>();
		testManifolds<3, Eigen::Dynamic, float>();
		testManifolds<3, 3, float>();
	}
}

#ifdef TEST_SPLINES

#include "DiffManifoldBSplineTests.hpp"

namespace bsplines {

TEST(DiffManifoldBSplineTestSuite, testInitialization)
{
	TestSpline rbspline;
	TestSplineDS rbsplineDS(splineOrder);

	Eigen::VectorXd x = Eigen::VectorXd::Random(rows);
	rbspline.initConstantUniformSpline(minTime, maxTime, numberOfSegments, x);

	SM_ASSERT_EQ(std::runtime_error, rbspline.getSplineOrder(), splineOrder, "");
	SM_ASSERT_EQ(std::runtime_error, rbsplineDS.getSplineOrder(), splineOrder, "");

	SM_ASSERT_EQ(std::runtime_error, rbspline.getAbsoluteNumberOfSegments() , numberOfSegments + splineOrder * 2 - 1, "");
	SM_ASSERT_EQ(std::runtime_error, maxTime, rbspline.getMaxTime(), "");
	SM_ASSERT_EQ(std::runtime_error, minTime, rbspline.getMinTime(), "");

	SM_ASSERT_EQ(std::runtime_error, minTime, rbspline.getTimeInterval().first, "");
	SM_ASSERT_EQ(std::runtime_error, maxTime, rbspline.getTimeInterval().second, "");

	SM_ASSERT_EQ(std::runtime_error, rbspline.firstRelevantSegment().getKnot(), rbspline.getEvaluatorAt<0>(minTime)._firstRelevantControlVertexIt.getKnot(),"");

	SM_ASSERT_EQ(std::runtime_error, minTime, rbspline.getEvaluatorAt<0>(minTime).getKnot(),"");
	SM_ASSERT_EQ(std::runtime_error, minTime, rbspline.getEvaluatorAt<0>(minTime + 0.01).getKnot(), "");
	SM_ASSERT_EQ(std::runtime_error, rbspline.getEvaluatorAt<0>(maxTime - 0.01).getKnot(), rbspline.getEvaluatorAt<0>(maxTime).getKnot(), ""); //TODO remove this irregularity
	SM_ASSERT_GT(std::runtime_error, maxTime, rbspline.getEvaluatorAt<0>(maxTime - 0.01).getKnot(), "");

	TestSpline::SegmentConstIterator i = rbspline.getAbsoluteBegin(), end = rbspline.getAbsoluteEnd();

	Eigen::MatrixXd zeroM = Eigen::MatrixXd::Zero(splineOrder, splineOrder);

	int c = 0;
	double lastKnot = minTime;

	BSpline bspline = BSpline(splineOrder);
	bspline.initConstantSpline(minTime, maxTime, numberOfSegments, x);
	std::vector<double> knots = bspline.knots();

	for(; i != end; i++){
		assert(i->getControlVertex() == x);

		double tKnot = i.getKnot();

		SM_ASSERT_NEAR(std::runtime_error, knots[c], tKnot, 1e-15, " c = " << c);

		if(!(c < splineOrder - 1|| c >= numberOfSegments + splineOrder - 1)){
			if(tKnot != maxTime){
				SM_ASSERT_EQ(std::runtime_error, rbspline.getEvaluatorAt(tKnot).getKnot(), tKnot, "");
				SM_ASSERT_EQ(std::runtime_error, rbspline.getEvaluatorAt(tKnot + 1e-13).getKnot(), rbspline.getEvaluatorAt(tKnot).getKnot(), "");
			}
			else{
				SM_ASSERT_EQ(std::runtime_error, rbspline.getEvaluatorAt(tKnot).getKnot(), lastKnot, "");
			}
			if(tKnot != minTime)
				SM_ASSERT_EQ(std::runtime_error, rbspline.getEvaluatorAt(tKnot - 1e-13).getKnot(), lastKnot, "");

			lastKnot = rbspline.getEvaluatorAt(tKnot).getKnot();
		}
		c++;
	}
}

TEST(DiffManifoldBSplineTestSuite, testGetBi)
{
	const int numTimeSteps = 100;

	BSpline bspline = BSpline(splineOrder);
	TestSpline rbspline(splineOrder);
	TestSplineDS rbsplineDS(splineOrder);

	rbspline.initConstantUniformSpline(minTime, maxTime, numberOfSegments, zero);
	rbsplineDS.initConstantUniformSpline(minTime, maxTime, numberOfSegments, zero);

	bspline.initConstantSpline(minTime, maxTime, numberOfSegments, zero);
	copyKnots(rbspline, bspline);

	for(int i = 0; i <= numTimeSteps; i ++) {
		double t = minTime + duration * ((double) i / numTimeSteps);
		TestSpline::Evaluator<0> evaluator = rbspline.getEvaluatorAt(t);
		Eigen::VectorXd localBiVector = evaluator.getLocalBi();
		SM_ASSERT_EQ(std::runtime_error, localBiVector.rows(), splineOrder, "");
		SM_ASSERT_NEAR(std::runtime_error, localBiVector.sum(), 1.0, 1e-10, "the bis at a given time should always sum up to 1")

		Eigen::VectorXd localCumulativeBiVector = evaluator.getLocalCumulativeBi();
		SM_ASSERT_EQ(std::runtime_error, localCumulativeBiVector.rows(), splineOrder, "");

		TestSpline::SegmentConstIterator firstIndex = evaluator.getFirstRelevantSegmentIterator();
		TestSpline::SegmentConstIterator lastIndex = evaluator.getLastRelevantSegmentIterator();
		SM_ASSERT_EQ(std::runtime_error, getIteratorDistance(firstIndex, lastIndex, (TestSpline::SegmentConstIterator) rbspline.end()), splineOrder - 1, "the distance between first and last local control point has to be exactly the spline order - 1. (at t = " << t << ")");

		SM_ASSERT_NEAR(std::runtime_error, (double)localCumulativeBiVector[0], 1.0, 1e-10, "");

		sm::eigen::assertNear(bspline.getLocalBiVector(t), rbspline.getEvaluatorAt(t).getLocalBi(), 1E-9, SM_SOURCE_FILE_POS);
		sm::eigen::assertNear(bspline.getLocalBiVector(t), rbsplineDS.getEvaluatorAt(t).getLocalBi(), 1E-9, SM_SOURCE_FILE_POS);
		sm::eigen::assertNear(bspline.getLocalCumulativeBiVector(t), rbspline.getEvaluatorAt(t).getLocalCumulativeBi(), 1E-9, SM_SOURCE_FILE_POS);
		sm::eigen::assertNear(bspline.getLocalCumulativeBiVector(t), rbsplineDS.getEvaluatorAt(t).getLocalCumulativeBi(), 1E-9, SM_SOURCE_FILE_POS);
		sm::eigen::assertNear(bspline.eval(t), rbspline.getEvaluatorAt(t).eval(), 1e-9, SM_SOURCE_FILE_POS);

		for(int j = 0, n = splineOrder; j < n; j++){
			if(j > 0){
				SM_ASSERT_LE(std::runtime_error, localCumulativeBiVector[j], localCumulativeBiVector[j - 1] + 1E-11, " at t = " << t << ", j = " << j);
			}
			SM_ASSERT_NEAR(std::runtime_error, (double)localCumulativeBiVector[j], (double)localBiVector.segment(j, splineOrder - j).sum(), 1e-13, "cumulativeBiVector must be the sum of the localBiVector where it overlaps, but it is not at " << j << " (localBiVector=" << localBiVector << ")");
		}
	}
}

TEST(DiffManifoldBSplineTestSuite, testLongTimeBSplineCompilation)
{
	TestSplineLongTime rbspline;
	TestSplineLongTime::point_t p = zero;
	p[0] = 1.0;
	rbspline.initConstantUniformSpline(minTimeLong, maxTimeLong, numberOfSegments, p);
	sm::eigen::assertEqual(rbspline.getEvaluatorAt(maxTimeLong).eval(), p, SM_SOURCE_FILE_POS);
}

} //namespace bsplines

#include "UnitQuaternionBSplineTests.cpp"
#include "EuclideanBSplineTests.cpp"

#ifdef SPEEDMEASURE
TEST(ZLASTDiffManifoldBSplineTestSuite, printTimings)
{
	sm::timing::Timing::print(std::cout);
}
#endif

#endif
