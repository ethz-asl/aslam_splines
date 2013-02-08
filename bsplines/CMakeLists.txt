cmake_minimum_required(VERSION 2.4.6)
include($ENV{ROS_ROOT}/core/rosbuild/rosbuild.cmake)

# Set the build type.  Options are:
#  Coverage       : w/ debug symbols, w/o optimization, w/ code-coverage
#  Debug          : w/ debug symbols, w/o optimization
#  Release        : w/o debug symbols, w/ optimization
#  RelWithDebInfo : w/ debug symbols, w/ optimization
#  MinSizeRel     : w/o debug symbols, w/ optimization, stripped binaries
#set(ROS_BUILD_TYPE Debug)
#set(ROS_BUILD_TYPE RelWithDebInfo)
set(ROS_BUILD_TYPE Release)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/../cmake)
rosbuild_init()
find_package(Eigen REQUIRED)

#set the default path for built executables to the "bin" directory
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/bin)
#set the default path for built libraries to the "lib" directory
set(LIBRARY_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/lib)

FIND_PACKAGE(PythonLibs REQUIRED)
INCLUDE_DIRECTORIES(${PYTHON_INCLUDE_DIRS})

#uncomment if you have defined messages
#rosbuild_genmsg()
#uncomment if you have defined services
#rosbuild_gensrv()

#common commands for building c++ executables and libraries
rosbuild_add_library(${PROJECT_NAME} 
  src/BSpline.cpp
  src/BSplinePose.cpp
)
#target_link_libraries(${PROJECT_NAME} another_library)
rosbuild_add_boost_directories()
#rosbuild_link_boost(${PROJECT_NAME} thread)
#rosbuild_add_executable(example examples/example.cpp)
#target_link_libraries(example ${PROJECT_NAME})




# For building the CHOLMOD / CSPARSE solvers 
SET(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${PROJECT_SOURCE_DIR}/../cmake)
FIND_PACKAGE(SuiteSparse REQUIRED)
include_directories(${CHOLMOD_INCLUDE_DIR})

message(${CHOLMOD_INCLUDE_DIR})

IF(NOT CHOLMOD_FOUND)
   message( FATAL_ERROR "Cholmod was not found")
ENDIF()

target_link_libraries(${PROJECT_NAME} ${CHOLMOD_LIBRARIES})


rosbuild_add_gtest(test/spline_tests 
  test/SplineTests.cpp
  test/BSplinePoseTests.cpp
#  test/DiffManifoldBSplineTests.cpp
  test/test_main.cpp
  )

target_link_libraries(test/spline_tests ${PROJECT_NAME})

#rosbuild_link_boost(test/spline_tests thread )


# Set up the python exports.
SET(PY_PROJECT_NAME bsplines)
SET(PY_PACKAGE_DIR python/bsplines)

rosbuild_find_ros_package(numpy_eigen)
include(${numpy_eigen_PACKAGE_PATH}/cmake/add_python_export_library.cmake)
add_python_export_library(${PY_PROJECT_NAME} ${PY_PACKAGE_DIR}
  src/SplinePython.cpp
  src/BSplinePython.cpp
  src/DiffManifoldBSplinePython.cpp
  src/BSplinePosePython.cpp
  )

rosbuild_link_boost(${PY_PROJECT_NAME} python thread filesystem)
target_link_libraries(${PY_PROJECT_NAME} ${PROJECT_NAME})




#rosbuild_add_pyunit( 
#  test/BSplinePoseTests.py
#)

#rosbuild_add_pyunit( 
#  test/BSplineTests.py TIMEOUT 600
#)

#rosbuild_add_pyunit( 
#  test/SplineTests.py
#)




SET(DOCSOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/doc/)



#--------------------------#
# Documentation
OPTION(BUILD_DOCUMENTATION "Should we build the documentation" ON)

IF(BUILD_DOCUMENTATION)
  FIND_PACKAGE(Doxygen)
  IF(DOXYGEN)
    Set(DOCSOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/doc/)
    SET(DOXYGEN_IMAGE_DIRS "\"${DOCSOURCE_DIR}/figs/\"")
    SET(DOXYGEN_SOURCE_DIRS "\"${CMAKE_CURRENT_SOURCE_DIR}/include/asrl/splines/\" \"${CMAKE_CURRENT_SOURCE_DIR}/src\" \"${CMAKE_CURRENT_SOURCE_DIR}/doc\"")

    Set(HTML_DIR ${CMAKE_CURRENT_SOURCE_DIR}/doc/html)
    SET(DOXYGEN_QUIET YES)
    MAKE_DIRECTORY(${HTML_DIR})
    SET( HAVE_DOT NO )
    
    CONFIGURE_FILE(${DOCSOURCE_DIR}/doxygen.config.in ${DOCSOURCE_DIR}/doxygen.config IMMEDIATE)
    
    ADD_CUSTOM_TARGET(doc ALL ${DOXYGEN} --silent ${DOCSOURCE_DIR}/doxygen.config)
  ENDIF(DOXYGEN)
ENDIF(BUILD_DOCUMENTATION)
#-------------------------#
