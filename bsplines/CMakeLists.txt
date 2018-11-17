cmake_minimum_required(VERSION 2.8.3)
project(bsplines)

find_package(catkin_simple REQUIRED)

catkin_simple(ALL_DEPS_REQUIRED)

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
  src/BSpline.cpp
  src/BSplinePose.cpp
)

if(CATKIN_ENABLE_TESTING)
  catkin_add_gtest(${PROJECT_NAME}_tests
    test/SplineTests.cpp
    test/BSplinePoseTests.cpp
    test/DiffManifoldBSplineTests.cpp
    test/test_main.cpp
  )
  if(TARGET ${PROJECT_NAME}_tests)
    target_link_libraries(${PROJECT_NAME}_tests ${PROJECT_NAME})
  endif()
endif()

cs_install()
cs_export()
