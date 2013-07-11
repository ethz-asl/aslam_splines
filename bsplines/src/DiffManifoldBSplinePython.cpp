#include <bsplines/DiffManifoldBSpline.hpp>
#include <bsplines/EuclideanBSpline.hpp>
#include <bsplines/UnitQuaternionBSpline.hpp>
#include <bsplines/BSplineFitter.hpp>

#include <numpy_eigen/boost_python_headers.hpp>
using namespace bsplines;
using namespace boost::python;

typedef EuclideanBSpline<>::TYPE PEuclideanBSpline;
typedef UnitQuaternionBSpline<>::TYPE PUnitQuaternionBSpline;


template <typename TSpline>
struct BSplineImporter {
	typedef typename TSpline::template Evaluator<Eigen::Dynamic> Evaluator;
	typedef BSplineFitter<TSpline> Fitter;

	static void addControlVertices(TSpline * bsp, const Eigen::VectorXd knotTimes, const Eigen::MatrixXd controlVertices)
	{
		SM_ASSERT_EQ(std::runtime_error, controlVertices.rows(), knotTimes.size(), "There must be as many knot times as controlVertices!");

		for(int i = 0, n = controlVertices.rows(); i < n; i++){
			bsp->addControlVertex(knotTimes[i], controlVertices.row(i));
		}
	}

	static void setControlVertices(TSpline * bsp, const Eigen::MatrixXd controlVertices)
	{
		SM_ASSERT_EQ(std::runtime_error, controlVertices.rows(), bsp->getNumControlVertices(), "There must be as many knot times as controlVertices!");

		int i = 0;
		for(typename TSpline::SegmentIterator it = bsp->getAbsoluteBegin(), end = bsp->getAbsoluteEnd(); it != end && i < controlVertices.rows(); it ++){
			it->getControlVertex() = controlVertices.row(i);
			i++;
		}
	}

	static typename TSpline::point_t eval(const TSpline * bsp, typename TSpline::time_t t){
		return bsp->template getEvaluatorAt<0>(t).eval();
	}

	static typename TSpline::point_t evalD(const TSpline * bsp, typename TSpline::time_t t, int derivativeOrder){
		return bsp->template getEvaluatorAt<Eigen::Dynamic>(t).evalD(derivativeOrder);
	}

	static void initSpline(TSpline * bsp, const std::vector<typename TSpline::time_t> & times, const std::vector<typename TSpline::point_t> & points, int numberOfSegments, double lambda){
		BSplineFitter<TSpline>().initUniformSpline(*bsp, times, points, numberOfSegments, lambda);
	}

	// Function wrappers turn std::pairs into tuples.
	static boost::python::tuple timeInterval1(const TSpline * bsp)
	{
		std::pair<double,double> ti = bsp->getTimeInterval();
		return boost::python::make_tuple(ti.first,ti.second);
	}

	static boost::python::tuple timeInterval2(const TSpline * bsp, int i)
	{
		std::pair<double,double> ti = bsp->timeInterval(i);
		return boost::python::make_tuple(ti.first,ti.second);
	}


	static class_<Evaluator> importEvaluator(){
		return class_<Evaluator>("Evaluator", init<const TSpline &, const typename TSpline::time_t &>())
			.def("eval", &Evaluator::eval, "Evaluate the spline curve at the evaluators time")
			.def("evalD", &Evaluator::evalD, "Evaluate a spline curve derivative at the evaluators time")
			;
	}

	static class_<TSpline> import(const char * moduleName){
		class_<TSpline> tSplineClass = class_<TSpline>(moduleName, init<int>())
		.def("init", &TSpline::init)
		.def("splineOrder", &TSpline::getSplineOrder, "The order of the spline")
//		.def("polynomialDegree", &TSpline::polynomialDegree, "The degree of the polynomial spline")
		.def("minimumKnotsRequired", &TSpline::minimumKnotsRequired, "The minimum number of knots required based on the spline order")
//		.def("numCoefficientsRequired", &TSpline::numCoefficientsRequired, "The number of coefficients required for a specified number of valid time segments")
//		.def("numKnotsRequired", &TSpline::numKnotsRequired, "The number of knots required for a target number of valid time segments")
//		.def("numValidTimeSegments", numValidTimeSegments1, "The number of valid time segments for a given number of knots")
//		.def("numValidTimeSegments", numValidTimeSegments2, "The number of valid time segments for the current knot sequence")
		.def("addControlVertices", &addControlVertices, "Adds the knots and spline control vertices")
		.def("setControlVertices", &setControlVertices, "Sets the spline control vertices")
//		.def("knots", &TSpline::knotVector, "returns the current knot sequence")
//		.def("coefficients", &TSpline::coefficients, "returns the current coefficient matrix", return_value_policy<copy_const_reference>())
		.def("getMinTime", &TSpline::getMinTime, "The minimum time that the spline is well-defined on")
		.def("getMaxTime", &TSpline::getMaxTime, "The maximum time that the spline is well-defined on")
		.def("eval", eval, "Evaluate the spline curve at a point in time")
		.def("evalD", evalD, "Evaluate a spline curve derivative at a point in time")
		.def("getEvaluatorAt", &TSpline::template getEvaluatorAt<Eigen::Dynamic> , "Get a evaluator at a point in time")
//		.def("Phi", &TSpline::Phi, "Evaluate the local basis matrix at a point in time")
//		.def("localBasisMatrix", &TSpline::localBasisMatrix, "Evaluate the local basis matrix at a point in time")
//		.def("localCoefficientMatrix", &TSpline::localCoefficientMatrix, "Get the matrix of locally-active coefficients for a specified time in matrix form")
//		.def("localCoefficientVector", &TSpline::localCoefficientVector, "Get the stacked vector of locally-active coefficients for a specified time.")
		.def("initSpline", initSpline, "Initialize the spline to interpolate a set of points")
		.def("initSpline", &Fitter::initUniformSpline, "Initialize the spline to interpolate a set of points")
		.def("initUniformSpline", &Fitter::initUniformSpline, "Initialize the spline to interpolate a set of points")
		.def("initUniformSplineSparse", &Fitter::initUniformSplineSparse, "Initialize the spline to interpolate a set of points (Sparse Solution)")
		.def("initSplineSparse", &Fitter::initUniformSplineSparse, "Initialize the spline to interpolate a set of points (Sparse Solution)")
//		.def("basisMatrix", &TSpline::basisMatrix, "Get the basis matrix active on the ith time segment.", return_value_policy<copy_const_reference>())
		.def("timeInterval", &timeInterval1, "Returns a tuple with the time interval that the spline is well-defined on.")
//		.def("timeInterval", &timeInterval2, "Returns a tuple with the time interval of the ith segment.")
		.def("getTimeInterval", &timeInterval1, "Returns a tuple with the time interval that the spline is well-defined on.")
//		.def("getTimeInterval", &timeInterval2, "Returns a tuple with the time interval of the ith segment.")
//		.def("addCurveSegment", &TSpline::addCurveSegment, "Adds a curve segment on the right that interpolates the given point at the given time.")
//		.def("removeCurveSegment", &TSpline::removeCurveSegment, "removes a curve segment on the left")
//		.def("setLocalCoefficientVector", &TSpline::setLocalCoefficientVector, "Sets the local coefficient vector for a specified time")
//		.def("localVvCoefficientVectorIndices", &TSpline::localVvCoefficientVectorIndices, "")
//		.def("localCoefficientVectorIndices", &TSpline::localCoefficientVectorIndices, "For the elements of a local coefficient vector, this gets the indices into the full coefficient vector")
//		.def("setCoefficientVector", &TSpline::setCoefficientVector, "Sets the full stacked coefficient vector of the spline")
//		.def("setCoefficientMatrix", &TSpline::setCoefficientMatrix, "Sets the full coefficient matrix of the spline")
//		.def("addCurveSegment2", &TSpline::addCurveSegment2, "")
//		.def("initSpline2", &TSpline::initSpline2, "")
//		.def("Vi",&TSpline::Vi,"")
//		.def("Mi", &TSpline::Mi, "")
//		.def("Bij", &TSpline::Bij, "")
//		.def("U", &TSpline::U, "U(time, derivativeOrder)")
//		.def("u", &TSpline::u, "")
//		.def("Di", &TSpline::Di, "")
//		.def("Dii", &TSpline::Dii, "")
//		.def("getLocalBi", &TSpline::getLocalBiVector, "getLocalBi(time)")
//		.def("getLocalCumulativeBi", &TSpline::getLocalCumulativeBiVector, "getLocalCumulativeBi(time)")
//		.def("getBiFunction", &getBiFunction, "getBiFunction(time)")
//		.def("getCumulativeBiFunction", &getCumulativeBiFunction, "getBiFunction(time)")
//		.def("segmentIndex", &TSpline::segmentIndex, "")
//		.def("segmentQuadraticIntegral", &TSpline::segmentQuadraticIntegral, "")
//		.def("segmentQuadraticIntegralDiag", &TSpline::segmentQuadraticIntegralDiag, "")
//		.def("curveQuadraticIntegral", &TSpline::curveQuadraticIntegral, "")
//		.def("curveQuadraticIntegralDiag", &TSpline::curveQuadraticIntegralDiag, "")
//		.def("curveQuadraticIntegralSparse", &TSpline::curveQuadraticIntegralSparse, "")
//		.def("curveQuadraticIntegralDiagSparse", &TSpline::curveQuadraticIntegralDiagSparse, "")
//		.def("coefficientVectorLength", &TSpline::coefficientVectorLength, "")
		.def("initConstantUniformSpline", &TSpline::initConstantUniformSpline, "initConstantUniformSpline(double t_min, double t_max, int numSegments, const Eigen::VectorXd & constant")
		.def("initConstantSpline", &TSpline::initConstantUniformSpline, "initConstantSpline(double t_min, double t_max, int numSegments, const Eigen::VectorXd & constant")
		.def("getNumControlVertices", &TSpline::getNumControlVertices, "")
		.def("numVvCoefficients", &TSpline::getNumControlVertices, "")
		;

		return tSplineClass;
	}
};

void import_bspline_diff_manifold_python()
{
	{
		typedef EuclideanBSpline<>::TYPE Spline;
		typedef BSplineImporter<Spline> SplineImporter;
		scope spline = SplineImporter::import("EuclideanBSpline")
			.def("evalI", &Spline::evalIntegral, "")
			.def("evalIntegral", &Spline::evalIntegral, "")
		;
		SplineImporter::importEvaluator();
	}

	{
		typedef UnitQuaternionBSpline<>::TYPE Spline;
		typedef BSplineImporter<Spline> SplineImporter;
		scope spline = SplineImporter::import("UnitQuaternionBSpline");
		SplineImporter::importEvaluator()
			.def("evalAngularVelocity", &SplineImporter::Evaluator::evalAngularVelocity)
			.def("evalAngularAcceleration", &SplineImporter::Evaluator::evalAngularAcceleration)
		;
	}
}
