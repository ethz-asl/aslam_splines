#ifndef ASLAM_PYTHON_EXPORT_FRAME_HPP
#define ASLAM_PYTHON_EXPORT_FRAME_HPP
#include <sstream>
#include <aslam/Frame.hpp>
#include <aslam/backend/CovarianceReprojectionError.hpp>
#include <aslam/backend/HomogeneousExpression.hpp>
#include <aslam/backend/CameraDesignVariable.hpp>
#include <aslam/backend/Scalar.hpp>

namespace aslam {
  namespace python {



    template<typename CAMERA_GEOMETRY_T>
    void exportCovarianceReprojectionError(const std::string & camName)
    {
        std::string name = camName + "ReprojectionErrorAdaptiveCovariance";
      using namespace boost::python;
      using namespace aslam;
      using namespace aslam::backend;
      typedef CAMERA_GEOMETRY_T geometry_t;
      typedef DescriptorBase descriptor_t;
      typedef Frame<geometry_t> frame_t;
      typedef typename frame_t::keypoint_t keypoint_t;

      class_< CovarianceReprojectionError<frame_t>, boost::shared_ptr<CovarianceReprojectionError<frame_t> >, bases< ErrorTerm > >( name.c_str(),
    		  init<const frame_t * , int ,HomogeneousExpression, CameraDesignVariable<geometry_t>, aslam::splines::BSplinePoseDesignVariable*, aslam::backend::Scalar* >( (name + "( frame, keypointIndex, homogeneousPointExpression, CameraDesignVariable, bsplineDesignVariable, lineDelayDv)").c_str()) )
			.def("observationTime", &CovarianceReprojectionError<frame_t>::observationTime)
			.def("covarianceMatrix",  &CovarianceReprojectionError<frame_t>::covarianceMatrix)
					;


    }




  } // namespace python
} // namespace aslam


#endif /* ASLAM_PYTHON_EXPORT_FRAME_HPP */
