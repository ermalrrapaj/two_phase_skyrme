#include <iostream> 
#include <math.h> 

#include "EOSData.hpp" 
#include "EOSSkyrme.hpp" 
#include "Constants.hpp"
#include "GibbsPhaseConstruct.hpp"

int main() {
  EOSSkyrme eos;
  GibbsPhaseConstruct gibbs(eos);
  double TMeV = 2.0; 
  double delta = 0.01;
   
  // Scan over neutron chemical potentials 
  double NpLoG = 1.e-20;
  double NpHiG = 0.08;
  for (double mun = 0.0; mun < 8.8; mun += delta) {
    
    std::vector<EOSData> outDat;
    
    try { 
      outDat = gibbs.FindPhasePoint(TMeV/Constants::HBCFmMeV, 
          mun/Constants::HBCFmMeV, NpLoG, NpHiG); 
    } catch (...) {
      std::cerr << mun <<  " failed.\n";
      continue;
    }
    
    if (outDat[0].Nn()>outDat[0].Np()) {
      std::cout << outDat[0].Nn() << " "; 
      std::cout << outDat[1].Nn() << " "; 
      std::cout << outDat[0].Np() << " "; 
      std::cout << outDat[1].Np() << " ";
      std::cout << outDat[1].P()*Constants::HBCFmMeV << " ";
      std::cout << outDat[1].Mun()*Constants::HBCFmMeV << " ";
      std::cout << outDat[1].Mup()*Constants::HBCFmMeV << std::endl;
    } 
     
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
      std::cerr << mun <<  " failed.\n";
      continue;
    }
    
    if (outDat[0].Nn()>outDat[0].Np()) {
      std::cout << outDat[0].Nn() << " "; 
      std::cout << outDat[1].Nn() << " "; 
      std::cout << outDat[0].Np() << " "; 
      std::cout << outDat[1].Np() << " ";
      std::cout << outDat[1].P()*Constants::HBCFmMeV << " ";
      std::cout << outDat[1].Mun()*Constants::HBCFmMeV << " ";
      std::cout << outDat[1].Mup()*Constants::HBCFmMeV << std::endl;
    } 
     
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
      std::cerr << mup <<  " failed.\n";
      continue;
    }

    if (outDat[0].Nn()<outDat[0].Np()) {
      std::cout << outDat[0].Nn() << " "; 
      std::cout << outDat[1].Nn() << " "; 
      std::cout << outDat[0].Np() << " "; 
      std::cout << outDat[1].Np() << " ";
      std::cout << outDat[1].P()*Constants::HBCFmMeV << " ";
      std::cout << outDat[1].Mun()*Constants::HBCFmMeV << " ";
      std::cout << outDat[1].Mup()*Constants::HBCFmMeV << std::endl;
    }
     
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
      std::cerr << mup <<  " failed.\n";
      continue;
    }

    if (outDat[0].Nn()<outDat[0].Np()) {
      std::cout << outDat[0].Nn() << " "; 
      std::cout << outDat[1].Nn() << " "; 
      std::cout << outDat[0].Np() << " "; 
      std::cout << outDat[1].Np() << " ";
      std::cout << outDat[1].P()*Constants::HBCFmMeV << " ";
      std::cout << outDat[1].Mun()*Constants::HBCFmMeV << " ";
      std::cout << outDat[1].Mup()*Constants::HBCFmMeV << std::endl;
    }
     
    NnLoG = outDat[0].Nn()*0.95;
    NnHiG = outDat[1].Nn()*1.05;
  }
  
  return 0;
}

