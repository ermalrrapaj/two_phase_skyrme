/// \file GibbsPhaseConstruct.cpp
/// \authorr lroberts
/// \since Aug 18, 2015
///
/// \brief
///
///

#include <math.h>
#include <iostream> 
#include <algorithm> 

#include "GibbsPhaseConstruct.hpp" 
#include "MultiDimensionalRoot.hpp"
#include "OneDimensionalRoot.hpp"

GibbsPhaseConstruct::GibbsPhaseConstruct(const EOSBase& eos) {
  mpEos = eos.MakeUniquePtr();   
}


std::vector<std::pair<EOSData, EOSData>>
GibbsPhaseConstruct::FindPhaseRange(double T, bool doMun, 
    double muStart, double muEnd, double deltaMu, double NLoG, double NHiG) {
  
  std::vector<std::pair<EOSData, EOSData>> phasePoints; 
   
  for (double mu = muStart; mu < muEnd; mu += deltaMu) {
    
    std::vector<EOSData> outDat;
    
    try { 
      outDat = FindPhasePoint(T, mu, NLoG, NHiG, doMun); 
    } catch (...) {
      continue;
    }
         
    if (((outDat[0].Nn()>outDat[0].Np()) && doMun) 
        || ((outDat[0].Np()>outDat[0].Nn()) && !doMun)) 
      phasePoints.push_back(std::pair<EOSData, EOSData>(outDat[0], outDat[1]));
    
    if (doMun) { 
      NLoG = outDat[0].Np()*0.95;
      NHiG = outDat[1].Np()*1.05;
    } else { 
      NLoG = outDat[0].Nn()*0.95;
      NHiG = outDat[1].Nn()*1.05;
    }
  }
  
  return phasePoints; 

} 

std::vector<std::pair<EOSData, EOSData>>
GibbsPhaseConstruct::FindFixedTPhaseBoundary(double T) {
  
  std::vector<std::pair<EOSData, EOSData>> phaseBound; 
  
  // Scan over neutron chemical potentials 
  double NpLoG = 1.e-20;
  double NpHiG = 0.08;
  auto phaseRange = FindPhaseRange(T, true, 0.0, 5.0 * T, 0.05, NpLoG, NpHiG);
  phaseBound.insert(phaseBound.end(), phaseRange.begin(), phaseRange.end()); 
  
  phaseRange = FindPhaseRange(T, true, 0.0, -15.0 * T, -0.05, NpLoG, NpHiG);
  phaseBound.insert(phaseBound.end(), phaseRange.begin(), phaseRange.end()); 
  
  phaseRange = FindPhaseRange(T, false, 0.0, 5.0 * T, 0.05, NpLoG, NpHiG);
  phaseBound.insert(phaseBound.end(), phaseRange.begin(), phaseRange.end()); 
  
  phaseRange = FindPhaseRange(T, false, 0.0, -15.0 * T, -0.05, NpLoG, NpHiG);
  phaseBound.insert(phaseBound.end(), phaseRange.begin(), phaseRange.end()); 
  
  // Sort the phase boundary data by the neutron chemical potential
  std::sort (phaseBound.begin(), phaseBound.end(), 
      [](std::pair<EOSData, EOSData> a, std::pair<EOSData, EOSData> b) { 
      return ((a.first).Mun() < (b.first).Mun());});
    
  return phaseBound; 
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
  double (EOSData::*getPotential)(void) const;
  EOSData (*eosDat)(double, double, double);
  if (doMun) {
    eosDat = EOSData::InputFromTNpMun;
    eosCall = &EOSBase::FromNpMunAndT;
    getPotential = &EOSData::Mup;
  } else {
    eosDat = EOSData::InputFromTNnMup;
    eosCall = &EOSBase::FromNnMupAndT;
    getPotential = &EOSData::Mun;
  }  
  
  // Using these guesses, search for the phase boundary points
  bool initial;
  auto root_f = [this, &mu, &T, initial, &eosCall, &eosDat, &doMun, &getPotential]
      (std::vector<double> xx) -> std::vector<double> {
    EOSData eLo = (mpEos.get()->*eosCall)(eosDat(T, exp(xx[0]), mu)); 
    EOSData eHi = (mpEos.get()->*eosCall)(eosDat(T, exp(xx[0])+exp(xx[1]), mu));
    if (initial) {
      return {eHi.P() - eLo.P(), (eHi.*getPotential)() - (eLo.*getPotential)()}; 
    } else {  
      return {(eHi.P() - eLo.P())/(eLo.P() + 1.e-8), 
          ((eHi.*getPotential)() - (eLo.*getPotential)()) 
          / ((eHi.*getPotential)() + 1.e-8)}; 
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
