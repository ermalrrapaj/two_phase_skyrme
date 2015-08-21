/// \file GibbsPhaseConstruct.cpp
/// \authorr lroberts
/// \since Aug 18, 2015
///
/// \brief
///
///

#include <math.h>
#include <iostream> 

#include "GibbsPhaseConstruct.hpp" 
#include "MultiDimensionalRoot.hpp"
#include "OneDimensionalRoot.hpp"

GibbsPhaseConstruct::GibbsPhaseConstruct(const EOSBase& eos) {
  mpEos = eos.MakeUniquePtr();   
}
 
std::vector<EOSData> GibbsPhaseConstruct::FindPhasePoint(double T, double mu, 
    double NLoG, double NHiG, bool doMun) {
  
  if (NLoG>NHiG) {
    double tmp = NHiG;
    NHiG = NLoG;
    NLoG = tmp;  
  }

  // Declare pointers so we can switch between neutron and proton chemical
  // potentials 
  EOSData (EOSBase::*eosCall)(const EOSData&) const;
  EOSData (*eosDat)(double, double, double);
  if (doMun) {
    eosDat = EOSData::InputFromTNpMun;
    eosCall = &EOSBase::FromNpMunAndT;
  } else {
    eosDat = EOSData::InputFromTNnMup;
    eosCall = &EOSBase::FromNnMupAndT;
  }  
  
  // Using these guesses, search for the phase boundary points
  bool initial;
  auto root_f = [this, &mu, &T, initial, &eosCall, &eosDat, &doMun]
      (std::vector<double> xx) -> std::vector<double> {
    EOSData eLo = (mpEos.get()->*eosCall)(eosDat(T, exp(xx[0]), mu)); 
    EOSData eHi = (mpEos.get()->*eosCall)(eosDat(T, exp(xx[0])+exp(xx[1]), mu));
    if (doMun) {
      if (initial) {
        return {eHi.P() - eLo.P(), eHi.Mup() - eLo.Mup()}; 
      } else {  
        return {(eHi.P() - eLo.P())/(eLo.P() + 1.e-8), 
            (eHi.Mup() - eLo.Mup())/(eHi.Mup() + 1.e-8)}; 
      }
    } else {
      if (initial) {
        return {eHi.P() - eLo.P(), eHi.Mun() - eLo.Mun()}; 
      } else {  
        return {(eHi.P() - eLo.P())/(eLo.P() + 1.e-8), 
            (eHi.Mun() - eLo.Mun())/(eHi.Mun() + 1.e-8)}; 
      }
    }
  };

  // First get close with linear equations  
  MultiDimensionalRoot rootFinder(1.e-10, 100);
  initial = true;
  std::vector<double> logN = rootFinder(root_f, 
      {log(NLoG), log(NHiG-NLoG)}, 2); 

  // Now get to high precision with non-linear scaled version of equations
  initial = false;
  logN = rootFinder(root_f, logN, 2);
  
  EOSData eosLo = (mpEos.get()->*eosCall)(eosDat(T, exp(logN[0]), mu)); 
  EOSData eosHi = (mpEos.get()->*eosCall)(eosDat(T, exp(logN[0])+exp(logN[1]), mu));
  
  if (fabs(eosHi.Np()/eosLo.Np()-1.0) < 1.e-4 &&
      fabs(eosHi.Nn()/eosLo.Nn()-1.0) < 1.e-4) {
    throw std::runtime_error("GibbsPhaseConstruct converged to the same point.");
  }
  
  return {eosLo, eosHi}; 
}  
