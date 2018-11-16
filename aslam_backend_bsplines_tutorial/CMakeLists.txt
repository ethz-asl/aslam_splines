cmake_minimum_required(VERSION 2.8.3)
project(aslam_backend_bsplines_tutorial)

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
  src/ErrorTermMotionBST.cpp
  src/ErrorTermPriorBST.cpp
  src/exampleBST.cpp
)
target_link_libraries(${PROJECT_NAME} ${Boost_LIBRARIES})


cs_add_executable(${PROJECT_NAME}_exec
  src/exampleBST.cpp
)
target_link_libraries(${PROJECT_NAME}_exec ${PROJECT_NAME})

if(CATKIN_ENABLE_TESTING)
  catkin_add_gtest(${PROJECT_NAME}_tests
    test/test_main.cpp
    test/TestErrorTestMotionBST.cpp
  )
  if(TARGET ${PROJECT_NAME}_tests)
    target_link_libraries(${PROJECT_NAME}_tests ${PROJECT_NAME})
  endif()
endif()

cs_install()
cs_export()
