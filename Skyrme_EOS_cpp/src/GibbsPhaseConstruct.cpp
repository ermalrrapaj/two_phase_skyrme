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

std::vector<EOSData> GibbsPhaseConstruct::FindPhasePoint(double T, double mun) {
  
  // Using these guesses, search for the phase boundary points
  bool initial;
  auto root_f = [this, &mun, &T, initial](std::vector<double> xx) 
      -> std::vector<double> {
    EOSData eLo = mpEos->FromNpMunAndT(
        EOSData::InputFromTMunNp(T, mun, exp(xx[0]))); 
    EOSData eHi = mpEos->FromNpMunAndT(
        EOSData::InputFromTMunNp(T, mun, exp(xx[0]) + exp(xx[1])));
    if (initial) {
      return {eHi.P() - eLo.P(), eHi.Mup() - eLo.Mup()}; 
    } else {  
      return {(eHi.P() - eLo.P())/(eLo.P() + 1.e-8), 
          (eHi.Mup() - eLo.Mup())/(eHi.Mup() + 1.e-8)}; 
    }
  };
   
  // First get close with linear equations  
  MultiDimensionalRoot rootFinder(1.e-10, 100);
  initial = true;
  std::vector<double> logNp = rootFinder(root_f, {log(1.e-80), log(0.08)}, 2); 
  
  // Now get to high precision with non-linear scaled version of equations
  initial = false;
  logNp = rootFinder(root_f, logNp, 2); 
  
  return {mpEos->FromNpMunAndT(EOSData::InputFromTMunNp(T, mun, exp(logNp[0]))),
          mpEos->FromNpMunAndT(EOSData::InputFromTMunNp(T, mun, 
          exp(logNp[0]) + exp(logNp[1])))};
}  
