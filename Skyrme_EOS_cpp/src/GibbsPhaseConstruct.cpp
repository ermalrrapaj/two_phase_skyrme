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
  auto outDat = FindPhasePoint(2.0/197.3, 0.0); 
  std::cout << std::endl;
  std::cout << outDat[0].Nn() << " " << outDat[1].Nn() << std::endl;
  std::cout << outDat[0].Np() << " " << outDat[1].Np() << std::endl;
  std::cout << outDat[0].P() << " " << outDat[1].P() << std::endl;
  std::cout << outDat[0].Mup() << " " << outDat[1].Mup() << std::endl;
}

std::vector<EOSData> GibbsPhaseConstruct::FindPhasePoint(double T, double mun) {
  
  MultiDimensionalRoot rootFinder(1.e-10, 1000);

  // Find the zero pressure point for this chemical potential 
  auto zero_P = [this](double logNp) -> double { 
    EOSData eOut= mpEos->FromNpMunAndT(
        EOSData::InputFromTMunNp(T, mun, exp(logNp))); 
    return eOut.P()/(eOut.Nb()*eOut.T());  
  };
  double lNp_pzero = rootFinder([this](double logNp) -, {log(0.2)}, 1);


  
  // Using these guesses, search for the phase boundary points 
  auto root_f = [this, &mun, &T](std::vector<double> xx) 
      -> std::vector<double> {
    std::cout << exp(xx[0]) << " " << exp(xx[1]) << std::endl;
    EOSData eLo = mpEos->FromNpMunAndT(
        EOSData::InputFromTMunNp(T, mun, exp(xx[0]))); 
    EOSData eHi = mpEos->FromNpMunAndT(
        EOSData::InputFromTMunNp(T, mun, exp(xx[0]) + exp(xx[1])));
    return {eHi.P() - eLo.P(), eHi.Mup() - eLo.Mup()}; 
  };  
  
  std::vector<double> logNp = rootFinder(root_f, {log(1.e-80), log(0.08)}, 2); 
  
  return {mpEos->FromNpMunAndT(EOSData::InputFromTMunNp(T, mun, exp(logNp[0]))),
          mpEos->FromNpMunAndT(EOSData::InputFromTMunNp(T, mun, 
          exp(logNp[0]) + exp(logNp[1])))};
}  
