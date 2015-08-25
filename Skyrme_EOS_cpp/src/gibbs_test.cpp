#include <iostream> 
#include <math.h> 
#include <vector> 
#include <utility>
#include <algorithm> 

#include "EOSData.hpp" 
#include "EOSSkyrme.hpp" 
#include "Constants.hpp"
#include "GibbsPhaseConstruct.hpp"

int main() {
  EOSSkyrme eos;
  GibbsPhaseConstruct gibbs(eos);
  double TMeV = 2.0/Constants::HBCFmMeV; 
  double delta = 0.05;
  double HBC = 1.0; //Constants::HBCFmMeV;
  delta = 0.025 * TMeV;
  
  auto phaseBound2 = gibbs.FindFixedTPhaseBoundary(TMeV/HBC); 
  
  std::vector<std::pair<EOSData,EOSData>> phaseBound;
  
  // Scan over neutron chemical potentials 
  double NpLoG = 1.e-20;
  double NpHiG = 0.08;
  for (double mun = 0.0; mun < 4.4 * TMeV; mun += delta) {
    
    std::vector<EOSData> outDat;
    double T =  TMeV/HBC; 
    double mu = mun/HBC;
    try { 
      outDat = gibbs.FindPhasePoint(T, mu, NpLoG, NpHiG); 
    } catch (...) {
      //std::cerr << mun <<  " failed.\n";
      continue;
    }
         
    if (outDat[0].Nn()>outDat[0].Np()) 
        phaseBound.push_back(std::pair<EOSData, EOSData>(outDat[0], outDat[1]));
     
    NpLoG = outDat[0].Np()*0.95;
    NpHiG = outDat[1].Np()*1.05;
  } 
  
  NpLoG = 1.e-20;
  NpHiG = 0.08;
  for (double mun = 0.0; mun > -7.5 * TMeV; mun -= delta) {
    
    std::vector<EOSData> outDat;
    double T =  TMeV/HBC; 
    double mu = mun/HBC;
    try { 
      outDat = gibbs.FindPhasePoint(T, mu, NpLoG, NpHiG); 
    } catch (...) {
      //std::cerr << mun <<  " failed.\n";
      continue;
    }
    
    if (outDat[0].Nn()>outDat[0].Np()) 
        phaseBound.push_back(std::pair<EOSData, EOSData>(outDat[0], outDat[1]));
     
    NpLoG = outDat[0].Np()*0.95;
    NpHiG = outDat[1].Np()*1.05;
  } 
  std::cerr << phaseBound.size() << std::endl; 
  
  // Scan over proton chemical potentials
  double NnLoG = 1.e-20;
  double NnHiG = 0.08;
  for (double mup = 0.0; mup < 4.4 * TMeV; mup += delta) {
    
    std::vector<EOSData> outDat;
    
    double T =  TMeV/HBC; 
    double mu = mup/HBC;
    try { 
      outDat = gibbs.FindPhasePoint(T, mu, NnLoG, NnHiG, false); 
    } catch (...) {
      //std::cerr << mup <<  " failed.\n";
      continue;
    }

    if (outDat[0].Nn()<outDat[0].Np()) 
        phaseBound.push_back(std::pair<EOSData, EOSData>(outDat[0], outDat[1]));
    NnLoG = outDat[0].Nn()*0.95;
    NnHiG = outDat[1].Nn()*1.05;
  } 
  
  NnLoG = 1.e-20;
  NnHiG = 0.08;
  for (double mup = 0.0; mup > -7.5 * TMeV; mup -= delta) {
    
    std::vector<EOSData> outDat;
    
    double T =  TMeV/HBC; 
    double mu = mup/HBC;
    try { 
      outDat = gibbs.FindPhasePoint(T, mu, NnLoG, NnHiG, false); 
    } catch (...) {
      //std::cerr << mup <<  " failed.\n";
      continue;
    }

    if (outDat[0].Nn()<outDat[0].Np())
        phaseBound.push_back(std::pair<EOSData, EOSData>(outDat[0], outDat[1]));
     
    NnLoG = outDat[0].Nn()*0.95;
    NnHiG = outDat[1].Nn()*1.05;
  }
  std::cerr << phaseBound.size() << std::endl; 
  std::cerr << phaseBound2.size() << std::endl; 
  
  std::sort (phaseBound.begin(), phaseBound.end(), 
      [](std::pair<EOSData, EOSData> a, std::pair<EOSData, EOSData> b) { 
        return ((a.first).Mun() < (b.first).Mun());});
  
  for (auto &a : phaseBound2) {
    std::cout << (a.first).Nn() << " "; 
    std::cout << (a.second).Nn() << " "; 
    std::cout << (a.first).Np() << " "; 
    std::cout << (a.second).Np() << " ";
    std::cout << (a.first).P()*HBC << " ";
    std::cout << (a.first).Mun()*HBC << " ";
    std::cout << (a.first).Mup()*HBC << std::endl;
  }  
   
    
  return 0;
}

