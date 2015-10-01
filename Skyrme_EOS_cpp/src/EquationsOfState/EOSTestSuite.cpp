/// \file EOSTestSuite.hpp
/// \author lroberts
/// \since Sep 30, 2015
///
/// \brief
///
///

#include "EOSTestSuite.hpp" 
#include "EquationsOfState/EOSData.hpp" 
#include "EquationsOfState/EOSBase.hpp" 

#define SMALL 1.e-50
#define COMPARE(x,y) fabs((x - y)/(x+SMALL)) 
 
int EOSTestSuite::CheckThermodynamicConsistency(double nn, double np, double T) const { 
  return 1; 
}

int EOSTestSuite::CheckAnalyticDerivatives(double nn, double np, double T) const { 
  double delta = 1.e-6; 
  EOSData state = mpEos->FromNAndT(EOSData::InputFromTNnNp(T, nn, np));
  double errtot = 0.0;  
  // Check temperature derivatives
  {  
    EOSData sp = mpEos->FromNAndT(EOSData::InputFromTNnNp(T*(1.0 + delta), nn, np));
    EOSData sm = mpEos->FromNAndT(EOSData::InputFromTNnNp(T*(1.0 - delta), nn, np));
    double dpdt = (sp.P() - sm.P())/(sp.T() - sm.T());
    double dsdt = (sp.S() - sm.S())/(sp.T() - sm.T());
    double dedt = (sp.E() - sm.E())/(sp.T() - sm.T());
    errtot += COMPARE(dpdt,state.dPdT());
    errtot += COMPARE(dsdt,state.dSdT());
    errtot += COMPARE(dedt,state.dEdT());
    if (mVerbose) {
    std::cout << state.dPdT() << " " << dpdt << " " << COMPARE(dpdt,state.dPdT()) << std::endl; 
    std::cout << state.dSdT() << " " << dsdt << " " << COMPARE(dsdt,state.dSdT()) << std::endl; 
    std::cout << state.dEdT() << " " << dedt << " " << COMPARE(dedt,state.dEdT()) << std::endl; 
    }
  } 

  // Check Nn derivatives  
  {  
    EOSData sp = mpEos->FromNAndT(EOSData::InputFromTNnNp(T, nn*(1.0 + delta), np));
    EOSData sm = mpEos->FromNAndT(EOSData::InputFromTNnNp(T, nn*(1.0 - delta), np));
    double dpdnn = (sp.P() - sm.P())/(sp.Nn() - sm.Nn());
    double dsdnn = (sp.S() - sm.S())/(sp.Nn() - sm.Nn());
    double dednn = (sp.E() - sm.E())/(sp.Nn() - sm.Nn());
    errtot += COMPARE(dpdnn,state.dPdNn());
    errtot += COMPARE(dsdnn,state.dSdNn());
    errtot += COMPARE(dednn,state.dEdNn());
    if (mVerbose) {
    std::cout << state.dPdNn() << " " << dpdnn << " " << COMPARE(dpdnn,state.dPdNn()) << std::endl; 
    std::cout << state.dSdNn() << " " << dsdnn << " " << COMPARE(dsdnn,state.dSdNn()) << std::endl; 
    std::cout << state.dEdNn() << " " << dednn << " " << COMPARE(dednn,state.dEdNn()) << std::endl; 
    }
  } 
  
  // Check Np derivatives  
  {  
    EOSData sp = mpEos->FromNAndT(EOSData::InputFromTNnNp(T, nn, np*(1.0 + delta)));
    EOSData sm = mpEos->FromNAndT(EOSData::InputFromTNnNp(T, nn, np*(1.0 - delta)));
    double dpdnp = (sp.P() - sm.P())/(sp.Np() - sm.Np());
    double dsdnp = (sp.S() - sm.S())/(sp.Np() - sm.Np());
    double dednp = (sp.E() - sm.E())/(sp.Np() - sm.Np());
    errtot += COMPARE(dpdnp,state.dPdNp());
    errtot += COMPARE(dsdnp,state.dSdNp());
    errtot += COMPARE(dednp,state.dEdNp());
    if (mVerbose) {
    std::cout << state.dPdNp() << " " << dpdnp << " " << COMPARE(dpdnp,state.dPdNp()) << std::endl; 
    std::cout << state.dSdNp() << " " << dsdnp << " " << COMPARE(dsdnp,state.dSdNp()) << std::endl; 
    std::cout << state.dEdNp() << " " << dednp << " " << COMPARE(dednp,state.dEdNp()) << std::endl; 
    } 
  } 
  if (mVerbose) std::cout << errtot/9.0 << std::endl;
  if (mVerbose < mTol) return 0;
  return 1; 
}

int EOSTestSuite::CompressionTest(double Ye, double S) const {
  return 1; 
}

