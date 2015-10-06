/// \file EOSTestSuite.hpp
/// \author lroberts
/// \since Sep 30, 2015
///
/// \brief
///
///


#ifndef EOS_EOSTESTSUITE_HPP_
#define EOS_EOSTESTSUITE_HPP_

#include "EquationsOfState/EOSBase.hpp"

///
/// Run a suite of tests on an EoS to test for thermodynamic consitency, etc. 
/// 
class EOSTestSuite {
public:
  EOSTestSuite(const EOSBase& eos, double tol = 1.e-7, bool verbose = false) : 
      mpEos(eos.MakeUniquePtr()), mVerbose(verbose), mTol(tol) {} 
  
  int CheckThermodynamicConsistency(double T, double nn, double np) const; 
  int CheckAnalyticDerivatives(double T, double nn, double np) const; 
  int CompressionTest(double Ye, double S) const;

protected: 
  std::unique_ptr<EOSBase> mpEos; 
  bool mVerbose;
  double mTol;
};

#endif // EOS_EOSTESTSUITE_HPP_
