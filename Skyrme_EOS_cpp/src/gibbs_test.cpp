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
  double HBC = Constants::HBCFmMeV;
  
  for (double T = 1.0/HBC; T < 3.0/HBC; T += 1.0/HBC) {
    auto phaseBound = gibbs.FindFixedTPhaseBoundary(T); 
    
    std::cerr << T * HBC << " " << phaseBound.size() << std::endl; 
    // Write out the phase boundary
    std::cout << "#[1] Nn_1 (fm^-3)" << std::endl; 
    std::cout << "#[2] Nn_2 (fm^-3)" << std::endl; 
    std::cout << "#[3] Np_1 (fm^-3)" << std::endl; 
    std::cout << "#[4] Np_2 (fm^-3)" << std::endl; 
    std::cout << "#[5] P (MeV fm^-3)" << std::endl; 
    std::cout << "#[6] Mun (MeV fm^-3)" << std::endl; 
    std::cout << "#[7] Mup (MeV fm^-3)" << std::endl; 
    std::cout << "#[8] T (MeV)" << std::endl; 
    
    for (auto &a : phaseBound) {
      std::cout << (a.first).Nn() << " "; 
      std::cout << (a.second).Nn() << " "; 
      std::cout << (a.first).Np() << " "; 
      std::cout << (a.second).Np() << " ";
      std::cout << (a.first).P()*HBC << " ";
      std::cout << (a.first).Mun()*HBC << " ";
      std::cout << (a.first).Mup()*HBC << " ";
      std::cout << (a.first).T()*HBC << std::endl;
    }
    
     
    std::reverse(phaseBound.begin(), phaseBound.end());  
    
    for (auto &a : phaseBound) {
      std::cout << (a.second).Nn() << " "; 
      std::cout << (a.first).Nn() << " "; 
      std::cout << (a.second).Np() << " ";
      std::cout << (a.first).Np() << " "; 
      std::cout << (a.first).P()*HBC << " ";
      std::cout << (a.first).Mun()*HBC << " ";
      std::cout << (a.first).Mup()*HBC << " ";
      std::cout << (a.first).T()*HBC << std::endl;
    }
    std::cout << (phaseBound.back().first).Nn() << " "; 
    std::cout << (phaseBound.back().second).Nn() << " "; 
    std::cout << (phaseBound.back().first).Np() << " ";
    std::cout << (phaseBound.back().second).Np() << " "; 
    std::cout << (phaseBound.back().first).P()*HBC << " ";
    std::cout << (phaseBound.back().first).Mun()*HBC << " ";
    std::cout << (phaseBound.back().first).Mup()*HBC << " ";
    std::cout << (phaseBound.back().first).T()*HBC << std::endl;
    
    std::cout << std::endl; 
  
  }   
    
  return 0;
}

