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
  
  for (double T = 0.5/HBC; T < 25.0/HBC; T += 0.5/HBC) {
    auto phaseBound = gibbs.FindFixedTPhaseBoundary(T); 
    std::cerr << T * HBC << " " << phaseBound.size() << std::endl; 
    
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

    std::cout << std::endl; 
  
  }   
    
  return 0;
}

