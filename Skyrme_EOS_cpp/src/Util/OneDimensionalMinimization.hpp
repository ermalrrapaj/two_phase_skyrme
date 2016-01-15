/// \file OneDimensionalMinimization.hpp
/// \authorr lroberts
/// \since Aug 18, 2015
///
/// \brief
///
///

#ifndef EOS_ONEDIMENSIONALMINIMIZATION_HPP_
#define EOS_ONEDIMENSIONALMINIMIZATION_HPP_

#include <functional> 
#include <stdexcept> 
#include <gsl/gsl_errno.h> 
#include <gsl/gsl_math.h>
#include <gsl/gsl_min.h> 
#include <iostream> 
#include <string>
#include <sstream> 

class OneDimensionalMinimization { 
public:
  OneDimensionalMinimization(double tol=1.e-8, int maxIter=50) :
      mTol(tol), mMaxIter(maxIter) {} 

  template<class FUNCTION>
  double operator() (FUNCTION func, const double xgin, const double xloin, 
      const double xhiin, const bool max=false) {
    
    // Similar to stack overflow 13289311 
    gsl_function F;
    F.function = [] (double x, void * p)->double { 
      return (*static_cast<FUNCTION*>(p))(x);
    };
    
    if (max) { 
      F.function = [] (double x, void * p)->double { 
        return -(*static_cast<FUNCTION*>(p))(x);
      };
    }
    F.params = &func;
    
    gsl_set_error_handler_off();
    const gsl_min_fminimizer_type *T;
    gsl_min_fminimizer *s;
    
    T = gsl_min_fminimizer_brent; 
    s = gsl_min_fminimizer_alloc(T); 
    gsl_min_fminimizer_set(s, &F, xgin, xloin, xhiin); 
    
    int status, status1, status2; 
    int iter = 0; 
    double xg, xlo, xhi;
    double mAbsTol = 0.0;
    do { 
      iter++; 
      status = gsl_min_fminimizer_iterate(s);
      xg  = gsl_min_fminimizer_x_minimum(s);  
      xlo = gsl_min_fminimizer_x_lower(s);  
      xhi = gsl_min_fminimizer_x_upper(s); 
      status = gsl_min_test_interval(xlo, xhi, mAbsTol, mTol);  
    } while (status == GSL_CONTINUE && iter < mMaxIter);
    
    gsl_set_error_handler(NULL);
    gsl_min_fminimizer_free(s);
    
    return xg; 
  
  }

protected: 
  double mTol; 
  int mMaxIter; 

};
#endif // EOS_ONEDIMENSIONALMINIMIZATION_HPP_

