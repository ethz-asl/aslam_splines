cmake_minimum_required(VERSION 2.8.3)
project(aslam_bsplines_python)

find_package(catkin_simple REQUIRED)
catkin_simple()

include_directories(${Eigen_INCLUDE_DIRS})

# TODO: 
#FIND_PACKAGE(PythonLibs REQUIRED)
#INCLUDE_DIRECTORIES(${PYTHON_INCLUDE_DIRS})


# Set up the python exports.
SET(PY_PROJECT_NAME bsplines_python)
SET(PY_PACKAGE_DIR python/bsplines)

add_python_export_library(${PY_PROJECT_NAME} ${PY_PACKAGE_DIR}
  src/SplinePython.cpp
  src/BSplinePython.cpp
  src/DiffManifoldBSplinePython.cpp
  src/BSplinePosePython.cpp
)

#rosbuild_link_boost(${PY_PROJECT_NAME} python thread filesystem)
#target_link_libraries(${PY_PROJECT_NAME} ${PROJECT_NAME})


## Add nosetest based cpp test target.
catkin_add_nosetests(test/DiffManifoldBSpline.py)

cs_install()
cs_export()


