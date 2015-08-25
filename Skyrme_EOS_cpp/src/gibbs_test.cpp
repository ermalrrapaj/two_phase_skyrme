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
  double T = 2.0/HBC; 
  auto phaseBound2 = gibbs.FindFixedTPhaseBoundary(T); 
  
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

