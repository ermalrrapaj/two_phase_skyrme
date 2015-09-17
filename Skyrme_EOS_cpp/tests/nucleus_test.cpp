#include <iostream> 
#include <math.h> 
#include <vector>

#include "Util/Constants.hpp"
#include "EquationsOfState/EOSData.hpp" 
#include "EquationsOfState/EOSSkyrme.hpp" 
#include "EquationsOfState/GibbsPhaseConstruct.hpp"
#include "EquationsOfState/SkyrmeParameters.hpp"
#include "EquationsOfState/NucleusBase.hpp"

int main() {
  const double HBC = Constants::HBCFmMeV;
  
  EOSSkyrme eos;

  LDNucleus nuc(28, 56, eos); 
   
  for (double ln = -3.0; ln < log10(0.2); ln += 0.01) { 
    EOSData extState = 
        eos.FromNAndT(EOSData::InputFromTNbYe(0.1/HBC, pow(10.0, ln), 0.01));
    std::cout << pow(10.0, ln) << " " << HBC/56.0 * nuc.GetBindingEnergy(extState, 0.1) << std::endl; 
  }
  return 0;
}

