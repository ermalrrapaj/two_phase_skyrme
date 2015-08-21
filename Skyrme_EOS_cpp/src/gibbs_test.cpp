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
  double TMeV = 2.0; 
  double delta = 0.05;

  std::vector<std::pair<EOSData,EOSData>> phaseBound;
  
  // Scan over neutron chemical potentials 
  double NpLoG = 1.e-20;
  double NpHiG = 0.08;
  for (double mun = 0.0; mun < 8.8; mun += delta) {
    
    std::vector<EOSData> outDat;
    
    try { 
      outDat = gibbs.FindPhasePoint(TMeV/Constants::HBCFmMeV, 
          mun/Constants::HBCFmMeV, NpLoG, NpHiG); 
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
  for (double mun = 0.0; mun > -15.0; mun -= delta) {
    
    std::vector<EOSData> outDat;
    
    try { 
      outDat = gibbs.FindPhasePoint(TMeV/Constants::HBCFmMeV, 
          mun/Constants::HBCFmMeV, NpLoG, NpHiG); 
    } catch (...) {
      //std::cerr << mun <<  " failed.\n";
      continue;
    }
    
    if (outDat[0].Nn()>outDat[0].Np()) 
        phaseBound.push_back(std::pair<EOSData, EOSData>(outDat[0], outDat[1]));
     
    NpLoG = outDat[0].Np()*0.95;
    NpHiG = outDat[1].Np()*1.05;
  } 
  
  // Scan over proton chemical potentials
  double NnLoG = 1.e-20;
  double NnHiG = 0.08;
  for (double mup = 0.0; mup < 8.8; mup += delta) {
    
    std::vector<EOSData> outDat;
    
    try { 
      outDat = gibbs.FindPhasePoint(TMeV/Constants::HBCFmMeV, 
          mup/Constants::HBCFmMeV, NnLoG, NnHiG, false); 
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
  for (double mup = 0.0; mup > -15.0; mup -= delta) {
    
    std::vector<EOSData> outDat;
    
    try { 
      outDat = gibbs.FindPhasePoint(TMeV/Constants::HBCFmMeV, 
          mup/Constants::HBCFmMeV, NnLoG, NnHiG, false); 
    } catch (...) {
      //std::cerr << mup <<  " failed.\n";
      continue;
    }

    if (outDat[0].Nn()<outDat[0].Np())
        phaseBound.push_back(std::pair<EOSData, EOSData>(outDat[0], outDat[1]));
     
    NnLoG = outDat[0].Nn()*0.95;
    NnHiG = outDat[1].Nn()*1.05;
  }
  
  std::sort (phaseBound.begin(), phaseBound.end(), 
      [](std::pair<EOSData, EOSData> a, std::pair<EOSData, EOSData> b) { 
        return ((a.first).Mun() < (b.first).Mun());});
  
  for (auto &a : phaseBound) {
    std::cout << (a.first).Nn() << " "; 
    std::cout << (a.second).Nn() << " "; 
    std::cout << (a.first).Np() << " "; 
    std::cout << (a.second).Np() << " ";
    std::cout << (a.first).P()*Constants::HBCFmMeV << " ";
    std::cout << (a.first).Mun()*Constants::HBCFmMeV << " ";
    std::cout << (a.first).Mup()*Constants::HBCFmMeV << std::endl;
  }  
   
  return 0;
}

