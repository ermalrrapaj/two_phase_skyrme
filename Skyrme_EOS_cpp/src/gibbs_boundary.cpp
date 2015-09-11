#include <iostream>
#include <fstream> 
#include <math.h> 
#include <vector> 
#include <utility>
#include <algorithm> 

#include "EOSData.hpp" 
#include "EOSSkyrme.hpp" 
#include "Constants.hpp"
#include "GibbsPhaseConstruct.hpp"
#include "SkyrmeParameters.hpp"    

int main() {
  double HBC = Constants::HBCFmMeV;
  
  //EOSSkyrme eos;
	//std::vector<double> Ska35S2009{-172.485, 172.087, -1767.71, 0.282732, 
  //    12899.2, 0.413266, 0.35};
  EOSSkyrme eos = EOSSkyrme::FromErmalSkyrme(ErmalSkyrmeParameters::Ska35S2009);
  GibbsPhaseConstruct gibbs(eos, false);
  
  std::ofstream outFile; 
  outFile.open("gibbs_boundary.out");
   
  for (double T = 3.0/HBC; T <= 22.0/HBC; T += 1.0/HBC) {
    auto phaseBound = gibbs.FindFixedTPhaseBoundary(T); 
    
    std::cout << T * HBC << " " << phaseBound.size() << std::endl; 
    if (phaseBound.size() > 0) {
      // Write out the phase boundary
      outFile << "#[1] Nn_1 (fm^-3)" << std::endl; 
      outFile << "#[2] Nn_2 (fm^-3)" << std::endl; 
      outFile << "#[3] Np_1 (fm^-3)" << std::endl; 
      outFile << "#[4] Np_2 (fm^-3)" << std::endl; 
      outFile << "#[5] P (MeV fm^-3)" << std::endl; 
      outFile << "#[6] Mun (MeV fm^-3)" << std::endl; 
      outFile << "#[7] Mup (MeV fm^-3)" << std::endl; 
      outFile << "#[8] T (MeV)" << std::endl; 
      
      for (auto &a : phaseBound) {
        outFile << (a.first).Nn() << " "; 
        outFile << (a.second).Nn() << " "; 
        outFile << (a.first).Np() << " "; 
        outFile << (a.second).Np() << " ";
        outFile << (a.first).P()*HBC << " ";
        outFile << (a.first).Mun()*HBC << " ";
        outFile << (a.first).Mup()*HBC << " ";
        outFile << (a.first).T()*HBC << std::endl;
      }
      
      outFile << std::endl; 
    } 
  }  
   
  outFile.close(); 
  return 0;
}

