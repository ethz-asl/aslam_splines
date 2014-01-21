#include <bsplines/python/DiffManifoldBSplineExporter.hpp>
#include <aslam/splines/OPTUnitQuaternionBSpline.hpp>
#include <aslam/splines/OPTEuclideanBSpline.hpp>
#include <aslam/splines/OPTEuclideanBSpline.hpp>

using namespace bspline_exporter;
using namespace bsplines;

void exportOptBSplines() {
  BSplineExporter<aslam::splines::OPTBSpline<EuclideanBSpline<>::CONF>::BSpline>::exportEuclideanSpline("OptEuclideanBSpline");
  BSplineExporter<aslam::splines::OPTBSpline<UnitQuaternionBSpline<>::CONF>::BSpline>::exportUnitQuaternionSpline("OptUnitQuaternionBSpline");
}
