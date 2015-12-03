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
  
  EOSSkyrme eos  = EOSSkyrme::FreeGas();
  //EOSSkyrme eosInside; 
   
  std::vector<std::unique_ptr<NucleusBase>> nuclei; 
  //nuclei.push_back(StaticNucleus(28, 56, 56.0*8.0/HBC, {}, {}, 9.2*56.0 ).MakeUniquePtr());
  for (int i=50; i<62; i++) {
  //nuclei.push_back(LDNucleus(28, i, eosInside).MakeUniquePtr());
  nuclei.push_back(StaticNucleus(28, i, i*8.0/HBC, {}, {}, 9.2*i ).MakeUniquePtr());
  }
  EOSNSE nseEos(nuclei, eos);
  
  
  double np0, nn0;
  double T = 2.0/HBC;
  double ye = 0.4;
  for (np0 = 0.001; np0 < 0.064; np0 *=1.1) { 
  try {
	  nn0 = (1.0-ye)/ye*np0;
	EOSData  eostot = nseEos.GetTotalEOS(EOSData::InputFromTNnNp(T, nn0, np0)); 
  std::cout << nn0+np0 << " " << T*HBC << " " << ye << " ";
  std::cout  << eostot.Np() << " " << eostot.Nn() << " "<< std::endl;
  } catch(std::exception& e ) {
  std::cout << nn0+np0  << " " << T*HBC << " " << ye << " fail " << e.what() << std::endl;
  }
  }
  return 0;
}

