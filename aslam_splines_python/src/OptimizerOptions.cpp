#include <numpy_eigen/boost_python_headers.hpp>
#include <aslam/backend/OptimizerOptions.hpp>
#include <boost/shared_ptr.hpp>

void exportOptimizerOptions()
{
  using namespace boost::python;
  using namespace aslam::backend;
  class_<OptimizerOptions>("OptimizerOptions", init<>())
    .def_readwrite("convergenceDeltaJ",&OptimizerOptions::convergenceDeltaJ)
    .def_readwrite("convergenceDeltaX",&OptimizerOptions::convergenceDeltaX)
    .def_readwrite("levenbergMarquardtLambdaInit",&OptimizerOptions::levenbergMarquardtLambdaInit)   
    .def_readwrite("levenbergMarquardtLambdaBeta", &OptimizerOptions::levenbergMarquardtLambdaBeta)
    .def_readwrite("levenbergMarquardtLambdaP", &OptimizerOptions::levenbergMarquardtLambdaP)
    .def_readwrite("levenbergMarquardtLambdaMuInit",&OptimizerOptions::levenbergMarquardtLambdaMuInit)
    .def_readwrite("levenbergMarquardtEstimateLambdaScale",&OptimizerOptions::levenbergMarquardtEstimateLambdaScale)
    .def_readwrite("doLevenbergMarquardt",&OptimizerOptions::doLevenbergMarquardt) 
    .def_readwrite("doSchurComplement",&OptimizerOptions::doSchurComplement)
    .def_readwrite("maxIterations",&OptimizerOptions::maxIterations)
    .def_readwrite("verbose",&OptimizerOptions::verbose)
    .def_readwrite("linearSolver",&OptimizerOptions::linearSolver)
    .def_readwrite("resetSolverEveryIteration", &OptimizerOptions::resetSolverEveryIteration)
    ;

}
