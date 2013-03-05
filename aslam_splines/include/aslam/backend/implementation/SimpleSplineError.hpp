#include <aslam/backend/SimpleSplineError.hpp>
#include <stdio.h>
namespace aslam {
    namespace backend {
        

    
        template<class SPLINE_T>
        SimpleSplineError<SPLINE_T>::SimpleSplineError(spline_t* splineDV, expression_t* splineExpression, Eigen::Matrix<double, spline_t::Dimension,1> y, double t):
        _splineDV(splineDV), _splineExpression(splineExpression), _y(y), _t(t)
        {
            
            // Add the design variables to the error term:        
            std::vector<aslam::backend::DesignVariable*> dvV;
        	for ( unsigned int i = 0; i < _splineDV->numDesignVariables(); i++) {
        		dvV.push_back(_splineDV->designVariable(i));
        	}
        	setDesignVariables(dvV);
            
        }


        template<class SPLINE_T>
        SimpleSplineError<SPLINE_T>::~SimpleSplineError()
        {

            
            
        }


        /// \brief evaluate the error term and return the weighted squared error e^T invR e
        template<class SPLINE_T>
        double SimpleSplineError<SPLINE_T>::evaluateErrorImplementation()
        {
            Eigen::VectorXd error = (_splineDV->spline().eval(_t) - _y);
            parent_t::setError(error);
            return error.transpose()*error;
            
        }


        /// \brief evaluate the jacobians
        template<class SPLINE_T>
        void SimpleSplineError<SPLINE_T>::evaluateJacobiansImplementation()
        {
            
            _splineExpression->evaluateJacobians(parent_t::_jacobians);

            
            
        }
          
          
        
          

    } // namespace backend
} // namespace aslam
