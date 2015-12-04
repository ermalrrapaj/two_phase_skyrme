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
    
    int status; 
    int iter = 0; 
    double x_lo, x_hi;
    do { 
      iter++; 
      status = gsl_root_fsolver_iterate(s);
      x_lo = gsl_root_fsolver_x_lower(s);  
      x_hi = gsl_root_fsolver_x_upper(s); 
      status = gsl_root_test_interval(x_lo, x_hi, 0, mTol);  
    } while (status == GSL_CONTINUE && iter < mMaxIter);
    
    gsl_set_error_handler(NULL);
    
    if (iter >= mMaxIter) throw std::runtime_error("Root find did not converge." + std::to_string(iter));
    double xcen = 0.5*(x_lo + x_hi);
    if (fabs(func(xcen))>mTol*1.e2) 
        throw std::runtime_error("Root find did not converge." + std::to_string(func(xcen)));
    return xcen; 
  }

protected: 
  double mTol; 
  int mMaxIter; 

};
#endif // EOS_ONEDIMENSIONALROOT_HPP_

