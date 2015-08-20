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

std::vector<EOSData> GibbsPhaseConstruct::FindPhasePoint(double T, double mun, 
    double NpLoG, double NpHiG) {
  
  if (NpLoG>NpHiG) {
    double tmp = NpHiG;
    NpHiG = NpLoG;
    NpLoG = tmp;  
  }
   
  // Using these guesses, search for the phase boundary points
  bool initial;
  auto root_f = [this, &mun, &T, initial](std::vector<double> xx) 
      -> std::vector<double> {
    EOSData eLo = mpEos->FromNpMunAndT(
        EOSData::InputFromTNpMun(T, exp(xx[0]), mun)); 
    EOSData eHi = mpEos->FromNpMunAndT(
        EOSData::InputFromTNpMun(T, exp(xx[0]) + exp(xx[1]), mun));
    if (initial) {
      return {eHi.P() - eLo.P(), eHi.Mup() - eLo.Mup()}; 
    } else {  
      return {(eHi.P() - eLo.P())/(eLo.P() + 1.e-8), 
          (eHi.Mup() - eLo.Mup())/(eHi.Mup() + 1.e-8)}; 
    }
  };
   
  // First get close with linear equations  
  MultiDimensionalRoot rootFinder(1.e-10, 200);
  initial = true;
  std::vector<double> logNp = rootFinder(root_f, 
      {log(NpLoG), log(NpHiG-NpLoG)}, 2); 

  // Now get to high precision with non-linear scaled version of equations
  initial = false;
  logNp = rootFinder(root_f, logNp, 2);
  
  EOSData eosLo = mpEos->FromNpMunAndT(
      EOSData::InputFromTNpMun(T, exp(logNp[0]), mun));
  EOSData eosHi = mpEos->FromNpMunAndT(
      EOSData::InputFromTNpMun(T, exp(logNp[0]) + exp(logNp[1]), mun));
  
  if (fabs(eosHi.Np()/eosLo.Np()-1.0) < 1.e-4 &&
      fabs(eosHi.Nn()/eosLo.Nn()-1.0) < 1.e-4) {
    throw std::runtime_error("GibbsPhaseConstruct converged to the same point.");
  }
  
  return {eosLo, eosHi}; 
}  
