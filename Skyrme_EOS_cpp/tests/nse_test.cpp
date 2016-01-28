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
#include "EquationsOfState/LDNucleus.hpp"
#include "EquationsOfState/EOSNSE.hpp"

int main() {
  const double HBC = Constants::HBCFmMeV;
  
  EOSSkyrme eos;//  = EOSSkyrme::FreeGas();
  EOSSkyrme eosInside; 
   
  std::vector<std::unique_ptr<NucleusBase>> nuclei; 
  //nuclei.push_back(StaticNucleus(28, 56, 56.0*8.0/HBC, {}, {}, 9.2*56.0 ).MakeUniquePtr());
  for (int i=50; i<62; i++) {
  nuclei.push_back(LDNucleus(28, i, eosInside).MakeUniquePtr());
  //nuclei.push_back(StaticNucleus(28, i, i*8.0/HBC, {}, {}, 9.2*i ).MakeUniquePtr());
  }
  EOSNSE nseEos(nuclei, eos);
  
  double nb = 0.02;
  double T = 2.0/HBC;
  double ye = 0.4;
  for (nb = 0.001; nb < 0.16; nb *=1.1) { 
  try {
  auto densities = 
      nseEos.GetExteriorDensities(EOSData::InputFromTNbYe(T, nb, ye)); 
  
  std::cout << nb << " " << T*HBC << " " << ye << " ";
  std::cout 
    << 56.0*densities[2]/nb << " " 
    << densities[0]/nb << " " << densities[1]/nb << " "
    << densities[3] << " " << densities[4] << " " << densities[5] << std::endl;
  } catch(std::exception& e ) {
  std::cout << nb << " " << T*HBC << " " << ye << " fail " << e.what() << std::endl;
  }
  }
  return 0;
}

