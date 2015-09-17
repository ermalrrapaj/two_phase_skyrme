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
#include <cstdio>
#include <string>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multiroots.h> 

class MultiDRootException : public std::invalid_argument {
public:
  MultiDRootException(std::string message, int iter, int status, 
      std::vector<double> x, std::vector<double> f) : 
      mXfin(x), mFfin(f), mIter(iter), mStatus(status), 
      std::invalid_argument(message) {}
  std::vector<double> GetF() const { return mFfin;}; 
  std::vector<double> GetX() const { return mXfin;}; 
  double GetError() const {
    double err = 0.0;
    for (auto& f : mFfin) err += fabs(f);
    return err;
  }

  int GetIterations() const { return mIter;}; 
  int GetStatus() const { return mStatus;};
   
protected: 
  std::vector<double> mXfin;
  std::vector<double> mFfin;    
  int mIter, mStatus;

}; 

class MultiDimensionalRoot { 
public:
  MultiDimensionalRoot(double tol=1.e-8, int maxIter=50) :
      mTol(tol), mMaxIter(maxIter) {} 

  template<class FUNCTION>
  std::vector<double> operator() (FUNCTION func, const std::vector<double>& xg, 
      const int nFunc) {

    gsl_set_error_handler_off();
    // Build the GSL type functions from the passed function
    // Similar to stack overflow 13289311 
    gsl_multiroot_function F;
    F.n = nFunc; 
    F.f = [] (const gsl_vector * x, void * p, gsl_vector * f)->int { 
      std::vector<double> xin; 
      for (unsigned int i=0; i < x->size; ++i) {
        double xt = gsl_vector_get(x, i);
        if (xt != xt) return GSL_EDOM; 
        xin.push_back(xt);
      }
      std::vector<double> ff = (*static_cast<FUNCTION*>(p))(xin);
      
      for (unsigned int i=0; i < x->size; ++i) {
        if (ff[i] != ff[i]) return GSL_ERANGE;
        gsl_vector_set(f, i, ff[i]);
      }
      return GSL_SUCCESS;
    };
    F.params = &func;
   
    gsl_vector *x = gsl_vector_alloc(nFunc); 
    for (unsigned int i=0; i<nFunc; ++i) {
      gsl_vector_set(x, i, xg[i]); 
    }

    const gsl_multiroot_fsolver_type *T = gsl_multiroot_fsolver_hybrids; 
    //const gsl_multiroot_fsolver_type *T = gsl_multiroot_fsolver_hybrid; 
    //const gsl_multiroot_fsolver_type *T = gsl_multiroot_fsolver_dnewton; 
    gsl_multiroot_fsolver *s = gsl_multiroot_fsolver_alloc(T, nFunc);
    int status = gsl_multiroot_fsolver_set(s, &F, x);
     
    //if (status) printf(" GSL Error: %s\n", gsl_strerror(status)); 
    
    int iter = 0; 
    do { 
      iter++; 
      status = gsl_multiroot_fsolver_iterate(s);
      if (status) break; // Solver is stuck  
      status = gsl_multiroot_test_residual(s->f, mTol);  
    } while (status == GSL_CONTINUE && iter < mMaxIter);
    
    gsl_set_error_handler(NULL);
    
    if (iter >= mMaxIter || status) {
        std::vector<double> ferr, xx;  
        for (int i=0; i<nFunc; i++) {
          ferr.push_back(gsl_vector_get(s->f, i));
          xx.push_back(gsl_vector_get(s->x, i));
        }
        throw MultiDRootException("Multi root find did not converge", iter,
            status, xx, ferr);
    }
    
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
