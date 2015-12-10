/// \file OneDimensionalRoot.hpp
/// \authorr lroberts
/// \since Aug 18, 2015
///
/// \brief
///
///

#ifndef EOS_ONEDIMENSIONALROOT_HPP_
#define EOS_ONEDIMENSIONALROOT_HPP_

#include <functional> 
#include <stdexcept> 
#include <gsl/gsl_errno.h> 
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h> 
#include <iostream> 
#include <string>
#include <sstream> 

class OneDimensionalRoot { 
public:
  OneDimensionalRoot(double tol=1.e-8, int maxIter=50) :
      mTol(tol), mMaxIter(maxIter) {} 

  template<class FUNCTION>
  double operator() (FUNCTION func, double xlo, double xhi) {
    
    // Similar to stack overflow 13289311 
    gsl_function F; 
    F.function = [] (double x, void * p)->double { 
      return (*static_cast<FUNCTION*>(p))(x);
    };
    F.params = &func;
    
    gsl_set_error_handler_off();
    const gsl_root_fsolver_type *T = gsl_root_fsolver_brent; 
    gsl_root_fsolver *s = gsl_root_fsolver_alloc(T); 
    gsl_root_fsolver_set(s, &F, xlo, xhi); 
    if (func(xlo)*func(xhi)>0) {
      std::cerr << "Bad interval " << xlo << " " << xhi << " " << func(xlo) <<
          " " << func(xhi) << std::endl;
      throw std::runtime_error("Bad interval in one-d root find");
    } 
    int status, status2; 
    int iter = 0; 
    double x_lo, x_hi;
    do { 
      iter++; 
      status = gsl_root_fsolver_iterate(s);
      x_lo = gsl_root_fsolver_x_lower(s);  
      x_hi = gsl_root_fsolver_x_upper(s); 
      status = gsl_root_test_interval(x_lo, x_hi, 0, mTol);  
      status2 = gsl_root_test_residual(func(0.5*(x_lo+x_hi)), mTol);  
    } while ((status == GSL_CONTINUE || status2 ==GSL_CONTINUE) && iter < mMaxIter);
    
    gsl_set_error_handler(NULL);
    
    if (status != 0 || iter >= mMaxIter) {
      std::stringstream stout; 
      stout << "Root find did not converge " << x_lo << " " << x_hi << " " 
          << func(x_lo) << " " << func(x_hi) << std::endl;
      throw std::runtime_error(stout.str());
    }
    
    double xcen; 
    if (fabs(func(x_lo))<fabs(func(x_hi))) {
      xcen = x_lo;
    } else {
      xcen = x_hi;
    }
    if (fabs(func(xcen))>fabs(func(0.5*(x_lo+x_hi)))) xcen = 0.5*(x_lo+x_hi);
    
    if (fabs(func(xcen))>mTol*1.e2) { 
        std::stringstream stout; 
        stout << "Root find did not converge 2" << x_lo << " " << x_hi 
            << " " << (x_hi - x_lo)/x_lo
            << " " << func(x_lo) << " " << func(x_hi) << std::endl;
        throw std::runtime_error(stout.str());
    }
    
    return xcen; 
  
  }

protected: 
  double mTol; 
  int mMaxIter; 

};
#endif // EOS_ONEDIMENSIONALROOT_HPP_

