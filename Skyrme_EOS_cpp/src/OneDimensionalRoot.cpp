/// \file OneDimensionalRoot.cpp
/// \authorr lroberts
/// \since Aug 18, 2015
///
/// \brief
///
///
#include <stdexcept> 
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h> 

#include "OneDimensionalRoot.hpp" 
template<class FUNCTION> 
double OneDimensionalRoot::FindRoot(FUNCTION func, 
    double xlo, double xhi) {
  
  // Similar to stack overflow 13289311 
  gsl_function F; 
  F.function = [] (double x, void * p)->double { 
    return (*static_cast<FUNCTION*>(p))(x);
  };
  F.params = &func;
  
  const gsl_root_fsolver_type *T; 
  T = gsl_root_fsolver_brent; 
  
  gsl_root_fsolver *s;
  s = gsl_root_fsolver_alloc(T); 
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
  
  if (iter >= mMaxIter) throw std::runtime_error("Root find did not converge.");
  
  return 0.5*(x_lo + x_hi);
}
