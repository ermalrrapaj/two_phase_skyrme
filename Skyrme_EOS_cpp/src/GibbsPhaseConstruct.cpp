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
#include <iterator> 

#include "GibbsPhaseConstruct.hpp" 
#include "MultiDimensionalRoot.hpp"
#include "OneDimensionalRoot.hpp"

GibbsPhaseConstruct::GibbsPhaseConstruct(const EOSBase& eos) {
  mpEos = eos.MakeUniquePtr();   
}

EOSData GibbsPhaseConstruct::FromNAndT(const EOSData& eosIn) const {
  double T = eosIn.T(); 
  double nb = eosIn.Nb();
  double np = eosIn.Np(); 

  auto phaseBound = FindFixedTPhaseBoundary(T); 
  
  std::vector<double> u, npArr;
  u.reserve(phaseBound.size()); 
  npArr.reserve(phaseBound.size()); 
  
  for (auto &p : phaseBound) {
    double uC = (nb - p.first.Nb())/(p.second.Nb() - p.first.Nb());
    u.push_back(uC);
    npArr.push_back((1.0 - uC)*p.first.Np() + uC*p.second.Np()); 
  }

  // Find the bracketing points
  std::vector<double>::iterator lo = npArr.begin();
  std::vector<double>::iterator up = npArr.end();
  while (lo+1 != up) {
    auto mid = std::distance(lo, up)/2 + lo; 
    if (*mid > np) {
      up = mid;
    } else { 
      lo = mid;
    } 
  } 

}

std::vector<std::pair<EOSData, EOSData>>
GibbsPhaseConstruct::FindPhaseRange(double T, bool doMun, 
    double muStart, double muEnd, double deltaMu, double NLoG, double NHiG) const {
  
  std::vector<std::pair<EOSData, EOSData>> phasePoints;
   
  auto comp_func = [&](double a, double b) -> bool {
      if (muStart < muEnd) {return (a<b);} else {return (a>b);}};
  
  auto iter_func = [&](double dmu) -> void { 
    for (double mu = muStart; comp_func(mu, muEnd); mu += dmu) {
      
      std::pair<EOSData, EOSData> outDat;
      try { 
        outDat = FindPhasePoint(T, mu, NLoG, NHiG, doMun); 
      } catch (...) {
        continue;
      }

      phasePoints.push_back(outDat);
      
      if (doMun) { 
        NLoG = outDat.first.Np()*0.95;
        NHiG = outDat.second.Np()*1.05;
      } else { 
        NLoG = outDat.first.Nn()*0.95;
        NHiG = outDat.second.Nn()*1.05;
      }
    }
    return;
  };
  
  // Do the basic iteration 
  iter_func(deltaMu); 

  // ok, that was a total failure 
  if (phasePoints.size() == 0) return phasePoints; 
  
  // Work forward with smaller steps from the last point 
  double muStartOld = muStart; 
  double muEndOld = muEnd; 
  if (doMun) {
    muStart = phasePoints.back().first.Mun();
  } else {
    muStart = phasePoints.back().first.Mup();
  }
  muEnd = muStart + 5.0 * deltaMu;
  iter_func(deltaMu*0.1); 
  muStart = muStartOld;
  muEnd = muEndOld;
   
  // Now work back towards muStart from the first succesful point 
  muEnd = muStart; 
  if (doMun) {
    muStart = phasePoints[0].first.Mun(); 
    NLoG = phasePoints[0].first.Np(); 
    NHiG = phasePoints[0].second.Np(); 
  } else { 
    muStart = phasePoints[0].first.Mup();
    NLoG = phasePoints[0].first.Nn(); 
    NHiG = phasePoints[0].second.Nn(); 
  } 
  
  if ((muStart-muEnd)/(muEnd + 1.e-4) < -0.02
      || (muStart-muEnd)/(muEnd + 1.e-4) > 0.02) 
      iter_func(-deltaMu);
  
  return phasePoints; 

} 

std::vector<std::pair<EOSData, EOSData>>
GibbsPhaseConstruct::FindFixedTPhaseBoundary(double T, double NLoG, 
    double NHiG, double deltaMu) const {
  
  std::vector<std::pair<EOSData, EOSData>> phaseBound; 
  
  // Scan over chemical potentials 
  deltaMu = deltaMu * T;
  double murange = 25.0 * T;  
  auto phaseRange = FindPhaseRange(T, true, 0.0, murange, deltaMu, NLoG, NHiG);
  phaseBound.insert(phaseBound.end(), phaseRange.begin(), phaseRange.end()); 
  
  phaseRange = FindPhaseRange(T, false, 0.0, murange, deltaMu, NLoG, NHiG);
  phaseBound.insert(phaseBound.end(), phaseRange.begin(), phaseRange.end()); 
  
  phaseRange = FindPhaseRange(T, true,  0.0, -murange, -deltaMu, NLoG, NHiG);
  phaseBound.insert(phaseBound.end(), phaseRange.begin(), phaseRange.end()); 
  
  phaseRange = FindPhaseRange(T, false, 0.0, -murange, -deltaMu, NLoG, NHiG);
  phaseBound.insert(phaseBound.end(), phaseRange.begin(), phaseRange.end()); 
  
  // Sort the phase boundary data by the neutron chemical potential
  std::sort (phaseBound.begin(), phaseBound.end(), 
      [](std::pair<EOSData, EOSData> a, std::pair<EOSData, EOSData> b) { 
      return ((a.first).Mun() < (b.first).Mun());});
  
  // Now check for and remove (likely) bad points
  int i = 1; 
  while (i < phaseBound.size()) {
    double dnb = phaseBound[i-1].second.Nb()/phaseBound[i].second.Nb();
    if (dnb<0.80 || dnb>1.2) {
      phaseBound.erase(phaseBound.begin()+i);
    } else {
      i++;
    }
  }  

  return phaseBound; 
}

EOSData GibbsPhaseConstruct::GetState(const EOSData& eosIn, 
    const EOSData& lo, const EOSData& hi, double ug) const { 
  double T = eosIn.T();
   
  auto root_f = [this, T, &eosIn]
      (std::vector<double> xx) -> std::vector<double> {
    
    EOSData eLo = mpEos->FromNAndT(
        EOSData::InputFromTNnNp(T, exp(xx[0]), exp(xx[1]))); 
    EOSData eHi = mpEos->FromNAndT(
        EOSData::InputFromTNnNp(T, exp(xx[2]), exp(xx[3])));
    
    return { (eHi.P() - eLo.P()) / (eLo.P() + 1.e-10),
         (eHi.Mun() - eLo.Mun()) / (eLo.Mun() + 1.e-10),
         (eHi.Mup() - eLo.Mup()) / (eLo.Mup() + 1.e-10),
         ((1.0 - xx[4])*eLo.Np() + xx[4]*eHi.Np())/eosIn.Np() - 1.0,
         ((1.0 - xx[4])*eLo.Nn() + xx[4]*eHi.Nn())/eosIn.Nn() - 1.0};
  };
  
  MultiDimensionalRoot rootFinder = MultiDimensionalRoot(1.e-11, 100);
  
  auto pars = rootFinder(root_f, {log(lo.Nn()), log(lo.Np()), log(hi.Nn()), 
      log(hi.Np()), ug}, 2);

  return EOSData();
} 
 
std::pair<EOSData, EOSData> GibbsPhaseConstruct::FindPhasePoint(double T, 
    double mu, double NLoG, double NHiG, bool doMun) const {
  
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
      return {(eHi.P() - eLo.P())/(eLo.P() + 1.e-12), 
          ((eHi.*getPotential)() - (eLo.*getPotential)()) 
          / ((eHi.*getPotential)() + 1.e-12)}; 
    }
  };

  // First get close with linear equations  
  MultiDimensionalRoot rootFinder(1.e-8, 200);
  initial = true;
  std::vector<double> logN = rootFinder(root_f, 
      {log(NLoG), log(NHiG-NLoG)}, 2); 

  // Now get to high precision with non-linear scaled version of equations
  initial = false;
  rootFinder = MultiDimensionalRoot(1.e-11, 100);
  logN = rootFinder(root_f, logN, 2);
  
  EOSData eosLo = (mpEos.get()->*eosCall)(eosDat(T, exp(logN[0]), mu)); 
  EOSData eosHi = (mpEos.get()->*eosCall)(eosDat(T, exp(logN[0])+exp(logN[1]), mu));
  
  if (fabs(eosHi.Np()/eosLo.Np()-1.0) < 1.e-8 &&
      fabs(eosHi.Nn()/eosLo.Nn()-1.0) < 1.e-8) {
    throw std::runtime_error("GibbsPhaseConstruct converged to the same point.");
  }
  
  if (eosHi.Nb() > 1.e2) {
    throw std::runtime_error("GibbsPhaseConstruct converged to a huge density.");
  } 
  
  if (eosHi.Nb() < 1.e-3) {
    throw std::runtime_error("GibbsPhaseConstruct converged to a small density.");
  } 

  return std::pair<EOSData, EOSData>(eosLo, eosHi); 
}  

