/// \file EOSSkyrme.cpp
/// \author lroberts
/// \since Aug 18, 2015
///
/// \brief
///
///

#include <math.h>

#include "EOSSkyrme.hpp" 
#include "OneDimensionalRoot.hpp"

static double HBC = 197.3; 
static double MNUC = 949.565/HBC;
static double PI = 3.14159; 

extern "C" {
  double ifermi12_(double* scale_density);
  double zfermi32_(double* eta);
}

EOSSkyrme::EOSSkyrme() :
    mA(-286.1/HBC), 
    mB(-107.1/HBC), 
    mC(968.0/HBC),
    mD(0.0), 
    mF(0.0), 
    mG(0.0),
    mDelta(2.002) {}

EOSData EOSSkyrme::FromNpMunAndT(const EOSData& eosIn) const {
  
  auto root_func = [&eosIn, this](double nn)->double {  
      EOSData out = BaseEOSCall(eosIn.T(), nn, eosIn.Np()); 
      return (out.Mun() - eosIn.Mun()) / (eosIn.Mun() + 1.e-30);
  }; 
  
  OneDimensionalRoot rootFinder(1.e-10);
  double nn_lo = 1.e-30;
  double nn_hi = 1.e3;
  double nn = rootFinder.FindRoot(root_func, nn_lo, nn_hi);
  
  return BaseEOSCall(eosIn.T(), nn, eosIn.Np());  
}

EOSData EOSSkyrme::FromNAndT(const EOSData& eosIn) const {
  return BaseEOSCall(eosIn.T(), eosIn.Nn(), eosIn.Np()); 
} 

EOSData EOSSkyrme::BaseEOSCall(const double T, const double nn, 
    const double np) const {
  
  const double nt = nn + np; 
  const double xp = np/(nt+1.e-40);
  
  double momsp = 1.0 + mF*(nn + np) + mG*(nn - np);
  double momsn = 1.0 + mF*(nn + np) - mG*(nn - np);
  
  double invetan = 2.0*PI*PI * nn * pow(2.0*MNUC/momsn*T, 1.5); 
  double invetap = 2.0*PI*PI * np * pow(2.0*MNUC/momsp*T, 1.5); 
  double etan = ifermi12_(&invetan); 
  double etap = ifermi12_(&invetap); 
  double taup = pow(2.0*MNUC/momsp*T, 2.5)/(2.0*PI*PI) * zfermi32_(&etap);
  double taun = pow(2.0*MNUC/momsn*T, 2.5)/(2.0*PI*PI) * zfermi32_(&etan);
  
  double Up = (taup*(mF - mG) + taun*(mF + mG)) / (2.0*MNUC) 
      + 2.0*mA*nt + 4.0*mB*nn + mC*(1.0 + mDelta)*pow(nt, mDelta) 
      + 4.0*mD*pow(nt, mDelta-2.0)*(nn*nt + (mDelta-1.0)*nn*np);
  double Un = (taup*(mF + mG) + taun*(mF - mG)) / (2.0*MNUC) 
      + 2.0*mA*nt + 4.0*mB*np + mC*(1.0 + mDelta)*pow(nt, mDelta) 
      + 4.0*mD*pow(nt, mDelta-2.0)*(np*nt + (mDelta-1.0)*nn*np);
  
  double mup = etap*T + Up; 
  double mun = etan*T + Un; 
  
  double ee = taup*momsp/(2.0*MNUC) + taun*momsn/(2.0*MNUC) 
      + (mA + 4.0*mB*xp*(1.0-xp))*nt*nt 
      + (mC + 4.0*mD*xp*(1.0-xp))*pow(nt, mDelta + 1.0);
  double pp = (5.0/6.0*momsp - 0.5)/MNUC*taup + (5.0/6.0*momsn - 0.5)/MNUC*taun
      + (mA + 4.0*mB*xp*(1.0-xp))*nt*nt // Check that there shouldn't be a factor fo two here 
      + (mC + 4.0*mD*xp*(1.0-xp))*mDelta*pow(nt, mDelta + 1.0);
  double ss = 5.0/(6.0*MNUC*T)*(momsp*taup + momsn*taun) - nt*etap - nn*etan; 
   
  return EOSData::Output(T, nn, np, mun, mup, pp, ss/nt, ee/nt);  
}



