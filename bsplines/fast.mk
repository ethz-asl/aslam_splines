CONFIG=debug
#CXX=clang

INC=-Iinclude/ -I$(ROS_WORKSPACE)/schweizer_messer/numpy_eigen/include -I$(ROS_WORKSPACE)/ros_comm/tools/rostest/include -I$(ROS_WORKSPACE)/schweizer_messer/sm_timing/include -I$(ROS_WORKSPACE)/schweizer_messer/sm_kinematics/include -I$(ROS_WORKSPACE)/schweizer_messer/sm_eigen/include -I$(ROS_WORKSPACE)/schweizer_messer/sm_common/include -I/usr/include/eigen3 -I/usr/include/opencv -I$(ROS_WORKSPACE)/ros_comm/clients/cpp/roscpp/include -I$(ROS_WORKSPACE)/ros_comm/clients/cpp/roscpp/msg_gen/cpp/include -I$(ROS_WORKSPACE)/ros_comm/clients/cpp/roscpp/srv_gen/cpp/include -I$(ROS_WORKSPACE)/ros_comm/clients/cpp/roscpp_serialization/include -I$(ROS_WORKSPACE)/ros_comm/clients/cpp/roscpp_traits/include -I$(ROS_WORKSPACE)/ros_comm/utilities/xmlrpcpp/src -I$(ROS_WORKSPACE)/ros_comm/tools/rosconsole/include -I$(ROS_WORKSPACE)/ros_comm/utilities/rostime/include -I$(ROS_WORKSPACE)/ros_comm/utilities/cpp_common/include -I$(ROS_WORKSPACE)/ros_comm/messages/rosgraph_msgs/msg_gen/cpp/include -I$(ROS_WORKSPACE)/ros_comm/messages/std_msgs/include -I$(ROS_WORKSPACE)/ros_comm/messages/std_msgs/msg_gen/cpp/include -I$(ROS_WORKSPACE)/ros/core/roslib/msg_gen/cpp/include -I$(ROS_WORKSPACE)/ros/core/roslib/include -I$(ROS_WORKSPACE)/ros/tools/rospack -I$(ROS_WORKSPACE)/ros/tools/rospack/include -I/usr/include/python2.7 -I$(ROS_WORKSPACE)/aslam/sparse_block_matrix/include/ -I/usr/include/suitesparse -I$(ROS_WORKSPACE)/aslam/sparse_block_matrix/include/ 
LIB= -lgtest -L $(ROS_WORKSPACE)/schweizer_messer/sm_timing/lib -L $(ROS_WORKSPACE)/schweizer_messer/sm_kinematics/lib -lsm_timing -lsm_kinematics -L$(ROS_WORKSPACE)/aslam/sparse_block_matrix/lib -lsparse_block_matrix -lcholmod -lcxsparse

TESTS = $(OBJ_DIR)/DiffManifoldBSplineTests.o

CXXFLAGS= -Winvalid-pch -fPIC -W -Wall -fno-strict-aliasing -W -Wall -Wno-unused-parameter -pthread $(INC) 
CXXFLAGS+= -Wfatal-errors

ifeq ($(CONFIG),debug)
	CXXFLAGS += -g
else
	CXXFLAGS += -O3 -DSPEEDMEASURE=1
endif

HEADER=include/bsplines/*.hpp include/bsplines/implementation/*.hpp include/bsplines/manifolds/implementation/*.hpp include/bsplines/manifolds/*.hpp  test/*.hpp $(ROS_WORKSPACE)/schweizer_messer/sm_eigen/include/sm/eigen/*.hpp

OBJ_DIR=build/fast/$(CONFIG)
BIN_DIR=bin/test/$(CONFIG)

.SUFFIXES:

vpath %.o $(OBJ_DIR)
vpath %.cpp test/ src/

all: mkdirs tests

%.dir :
	mkdir -p $(dir $@)
	touch $@
	
mkdirs: $(addsuffix /.dir,$(OBJ_DIR) $(BIN_DIR))

test_part: $(BIN_DIR)/test_part

tests: $(BIN_DIR)/spline_tests_partial

$(OBJ_DIR)/DiffManifoldBSplineTests.o: DiffManifoldBSplineTests.cpp UnitQuaternionBSplineTests.cpp EuclideanBSplineTests.cpp NodeDistributedCacheTests.cpp NumericIntegratorTests.cpp $(HEADER) fast.mk
	@echo compiling $< to $@
	@LC_ALL=C $(CXX) $(CXXFLAGS) -c -o $@ $<
	
$(OBJ_DIR)/%.o : %.cpp
	@echo compiling $< to $@
	@LC_ALL=C $(CXX) $(CXXFLAGS) -c -o $@ $<

$(BIN_DIR)/spline_tests_partial: test_main.cpp $(TESTS) $(OBJ_DIR)/BSpline.o
	@echo compiling $^ to $@
	@LC_ALL=C $(CXX) $(CXXFLAGS) -o $@ $^ $(LIB)
	
run: $(BIN_DIR)/spline_tests_partial
	@echo running config $(CONFIG) 
	LD_LIBRARY_PATH="$(ROS_WORKSPACE)/schweizer_messer/sm_timing/lib:$(ROS_WORKSPACE)/schweizer_messer/sm_kinematics/lib:$(ROS_WORKSPACE)/aslam/sparse_block_matrix/lib" $(BIN_DIR)/spline_tests_partial 
	