#include <iostream>
#include <math.h> 
#include <vector> 
#include <utility>

#include "Util/Constants.hpp"
#include "EquationsOfState/EOSData.hpp" 
#include "EquationsOfState/EOSSkyrme.hpp" 
#include "EquationsOfState/GibbsPhaseConstruct.hpp"
#include "EquationsOfState/SkyrmeParameters.hpp"    

int main() {
  double HBC = Constants::HBCFmMeV;
  
  EOSSkyrme eos = EOSSkyrme::FromErmalSkyrme(ErmalSkyrmeParameters::Ska35S2009);
  GibbsPhaseConstruct gibbs(eos, false);
   
  double T = 3.0/HBC;
  auto phaseBound = gibbs.FindFixedTPhaseBoundary(T); 
  std::cout << T * HBC << " " << phaseBound.size() << std::endl; 
    
   
  if (phaseBound.size() < 100) return 1;
  return 0;
}

