#include <iostream> 
#include <math.h> 
#include <vector>
#include <memory>

#include "Util/Constants.hpp"
#include "EquationsOfState/EOSData.hpp" 
#include "EquationsOfState/EOSSkyrme.hpp" 
#include "EquationsOfState/GibbsPhaseConstruct.hpp"
#include "EquationsOfState/SkyrmeParameters.hpp"
#include "EquationsOfState/NucleusBase.hpp"
#include "EquationsOfState/EOSNSE.hpp"

int main() {
  const double HBC = Constants::HBCFmMeV;
  
  EOSSkyrme eos; // = EOSSkyrme::FreeGas();
  
  std::vector<std::unique_ptr<NucleusBase>> nuclei; 
  nuclei.push_back(StaticNucleus(28, 56, 56.0*8.0/197.3, {}, {}, 9.2*56.0 ).MakeUniquePtr());
  EOSNSE nseEos(nuclei, eos);
  
  double nb = 1.e-2;
  double T = 0.5/HBC;
  double ye = 0.4;
  for (nb = 0.01; nb < 0.2; nb *=1.05) { 
  try {
  auto densities = 
      nseEos.GetExteriorDensities(EOSData::InputFromTNbYe(T, nb, ye)); 
  
  std::cout << nb << " " << T*HBC << " " << ye << " ";
  std::cout 
    << 56.0*densities[2]/nb << " " 
    << densities[0] << " " << densities[1] << std::endl;
  } catch(...) {
  }
  }
  return 0;
}

