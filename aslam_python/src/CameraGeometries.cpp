#include <numpy_eigen/boost_python_headers.hpp>
#include <aslam/cameras.hpp>

#include <aslam/backend/CameraDesignVariable.hpp>
#include <aslam/python/ExportBackendExpressions.hpp>

using namespace aslam::cameras;
using namespace boost::python;

void exportCameraGeometryBase();
void exportCameraProjections();
void exportCameraShutters();



template<typename C>
void exportCameraDesignVariables(std::string name)
{

	using namespace aslam::backend;
	using namespace aslam::python;

	boost::python::class_<CameraDesignVariable<C>, boost::shared_ptr<CameraDesignVariable<C> > >((name + "DesignVariable").c_str() , boost::python::init< boost::shared_ptr<C> >())
			.def("euclideanToKeypoint", &CameraDesignVariable<C>::euclideanToKeypoint)
			.def("homogeneousToKeypoint", &CameraDesignVariable<C>::homogeneousToKeypoint)
			.def("setActive", &CameraDesignVariable<C>::setActive )
			.def("getDesignVariables", &getDesignVariables<CameraDesignVariable<C> >)
			.def("projectionDesignVariable", &CameraDesignVariable<C>::projectionDesignVariable)
			.def("distortionDesignVariable", &CameraDesignVariable<C>::distortionDesignVariable)
			.def("shutterDesignVariable", &CameraDesignVariable<C>::shutterDesignVariable)
			;

}




template<typename C>
void exportCameraGeometry(std::string name)
{
  typedef typename C::shutter_t shutter_t;
  typedef typename C::projection_t projection_t;
  typedef typename C::mask_t mask_t;

  shutter_t & (C::*shutter)() = &C::shutter;
  projection_t & (C::*projection)() = &C::projection;
  mask_t & (C::*mask)() = &C::mask;


  boost::python::class_<C, boost::shared_ptr<C>, boost::python::bases<CameraGeometryBase> >(name.c_str(), boost::python::init<>())
	.def(init<typename C::projection_t>())
	.def(init<typename C::projection_t, typename C::shutter_t > () )
	.def(init<typename C::projection_t, typename C::shutter_t, typename C::mask_t> () )
		 .def("shutter", shutter, return_internal_reference<>())
		 .def("projection", projection, return_internal_reference<>())
		 .def("mask", mask, return_internal_reference<>())
		 ;


  exportCameraDesignVariables<C>(name);

}




void exportCameraGeometries()
{
  exportCameraGeometryBase();
  exportCameraProjections();
  exportCameraShutters();
  //exportPinholeCameraGeometry();
  //exportPinholeRSCameraGeometry();
  //exportOmniCameraGeometry();
  exportCameraGeometry<PinholeCameraGeometry>("PinholeCameraGeometry");
  exportCameraGeometry<DistortedPinholeCameraGeometry>("DistortedPinholeCameraGeometry");
  exportCameraGeometry<PinholeRsCameraGeometry>("PinholeRsCameraGeometry");
  exportCameraGeometry<DistortedPinholeRsCameraGeometry>("DistortedPinholeRsCameraGeometry");

  exportCameraGeometry<OmniCameraGeometry>("OmniCameraGeometry");
  exportCameraGeometry<DistortedOmniCameraGeometry>("DistortedOmniCameraGeometry");

}
