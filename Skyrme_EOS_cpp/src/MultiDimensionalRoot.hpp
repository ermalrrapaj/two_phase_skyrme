/// \file MultiDimensionalRoot.hpp
/// \authorr lroberts
/// \since Aug 18, 2015
///
/// \brief
///
///

#ifndef EOS_MULTIDIMENSIONALROOT_HPP_
#define EOS_MULTIDIMENSIONALROOT_HPP_

#include <vector>
#include <functional> 
#include <stdexcept> 
#include <iostream>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multiroots.h> 

class MultiDimensionalRoot { 
public:
  MultiDimensionalRoot(double tol=1.e-8, int maxIter=50) :
      mTol(tol), mMaxIter(maxIter) {} 

  template<class FUNCTION>
  std::vector<double> operator() (FUNCTION func, const std::vector<double>& xg, 
      const int nFunc) {

    // Build the GSL type functions from the passed function
    // Similar to stack overflow 13289311 
    gsl_multiroot_function F;
    F.n = nFunc; 
    F.f = [] (const gsl_vector * x, void * p, gsl_vector * f)->int { 
      std::vector<double> xin; 
      for (unsigned int i=0; i < x->size; ++i) {
        xin.push_back(gsl_vector_get(x, i));
      }
      
      std::vector<double> ff = (*static_cast<FUNCTION*>(p))(xin);
      
      for (unsigned int i=0; i < x->size; ++i) {
        gsl_vector_set(f, i, ff[i]);
      }
      return GSL_SUCCESS;
    };
    F.params = &func;
   
    gsl_vector *x = gsl_vector_alloc(nFunc); 
    for (unsigned int i=0; i<nFunc; ++i) {
      gsl_vector_set(x, i, xg[i]); 
    }

    const gsl_multiroot_fsolver_type *T; 
    gsl_multiroot_fsolver *s;
    T = gsl_multiroot_fsolver_hybrids; 
    s = gsl_multiroot_fsolver_alloc(T, nFunc);
    gsl_multiroot_fsolver_set(s, &F, x);
     
    int status; 
    int iter = 0; 
    do { 
      iter++; 
      status = gsl_multiroot_fsolver_iterate(s);
      if (status) break; // Solver is stuck  
      status = gsl_multiroot_test_residual(s->f, mTol);  
    } while (status == GSL_CONTINUE && iter < mMaxIter);
    
    if (iter >= mMaxIter || status) 
        throw std::runtime_error("Multi root find did not converge.");
    
    std::vector<double> out; 
    for (unsigned int i=0; i<nFunc; ++i) {
      out.push_back(gsl_vector_get(s->x, i)); 
    }
    
    return out;
  }

protected: 
  double mTol; 
  int mMaxIter; 

};
#endif // EOS_MULTIDIMENSIONALROOT_HPP_
