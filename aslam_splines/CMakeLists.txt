cmake_minimum_required(VERSION 2.8.3)
project(aslam_splines)

find_package(catkin_simple REQUIRED)
catkin_simple(ALL_DEPS_REQUIRED)

find_package(Boost REQUIRED COMPONENTS system)
include_directories(${Boost_INCLUDE_DIRS})

# enable warnings
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra")
if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
  # deprecated-register warnings are cause by old Eigen versions
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-deprecated-register")
endif()

# enable C++11 support
if(CMAKE_VERSION VERSION_LESS "3.1")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
else()
  set(CMAKE_CXX_STANDARD 11)
endif()

cs_add_library(${PROJECT_NAME}
  src/BSplinePoseDesignVariable.cpp
  src/BSplineRSPoseDesignVariable.cpp
  src/BSplineExpressions.cpp
  src/BSplineRSExpressions.cpp
  src/EuclideanBSplineDesignVariable.cpp
)
target_link_libraries(${PROJECT_NAME} ${Boost_LIBRARIES})

catkin_add_gtest(${PROJECT_NAME}_tests
  test/test_main.cpp
  test/TestBSplineExpressions.cpp
  test/TestOPTBSpline.cpp
  test/TestErrors.cpp
)
target_link_libraries(${PROJECT_NAME}_tests ${PROJECT_NAME})

cs_install()
cs_export()
