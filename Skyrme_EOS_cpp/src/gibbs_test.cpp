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
  double NpLoG = 1.e-20;
  double NpHiG = 0.08;
   
  for (double mun = 0.0; mun < 8.8; mun += 0.005) {
    
    std::vector<EOSData> outDat;
    
    try { 
      outDat = gibbs.FindPhasePoint(TMeV/Constants::HBCFmMeV, 
          mun/Constants::HBCFmMeV, NpLoG, NpHiG); 
    } catch (...) {
      std::cerr << mun <<  " failed.\n";
      continue;
    }

    std::cout << outDat[0].Nn() << " "; 
    std::cout << outDat[1].Nn() << " "; 
    std::cout << outDat[0].Np() << " "; 
    std::cout << outDat[1].Np() << " ";
    std::cout << outDat[1].P()*Constants::HBCFmMeV << " ";
    std::cout << outDat[1].Mun()*Constants::HBCFmMeV << " ";
    std::cout << outDat[1].Mup()*Constants::HBCFmMeV << std::endl;
     
    NpLoG = outDat[0].Np()*0.95;
    NpHiG = outDat[1].Np()*1.05;
  } 
  return 0;
}

