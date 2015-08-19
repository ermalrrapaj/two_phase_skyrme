/// \file GibbsPhaseConstruct.cpp
/// \authorr lroberts
/// \since Aug 18, 2015
///
/// \brief
///
///

#include "GibbsPhaseConstruct.hpp" 
#include "MultiDimensionalRoot.hpp"

GibbsPhaseConstruct::GibbsPhaseConstruct(const EOSBase& eos) {
  mpEos = eos.MakeUniquePtr();   
}

//void GibbsPhaseConstruct::FindPhasePoint(double T, double mun) {
//  
//  MultiDimensionalRoot rootFinder(1.e-10); 
//  auto root_f = [&mpEos](std::vector<double> xx) -> std::vector<double> ff {
//    EOSData eLo = mpEos->FromMunNpandT(
//        EOSData::InputFromTMunN(T, mun, exp(xx[0]))); 
//    EOSData eHi = mpEos->FromMunNpandT(
//        EOSData::InputFromTMunN(T, mun, exp(xx[1])));
//    return {eHi.P() - eLo.P(), eHi.Mup() - eLo.Mup()}; 
//  };  
//  
//  
//}  
