# Import the numpy to Eigen type conversion.
import roslib; 
roslib.load_manifest('numpy_eigen'); import numpy_eigen
import os

isCompiled = False
pathToSo = os.path.dirname(os.path.realpath(__file__))
if os.path.isfile(os.path.join(pathToSo,"libaslam_splines_python.so")):    
    roslib.load_manifest('aslam_backend_python'); import aslam_backend
    roslib.load_manifest('bsplines'); import bsplines
    # Import the the C++ exports from your package library.
    from libaslam_splines_python import *
    # Import other files in the directory
    # from mypyfile import *
    isCompiled = True
else:
    print "Warning: the package aslam_splines_python is not compiled. Type 'rosmake aslam_splines_python' if you need this."
    PACKAGE_IS_NOT_COMPILED = True;


