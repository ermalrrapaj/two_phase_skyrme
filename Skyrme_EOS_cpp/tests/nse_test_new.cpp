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
  
  EOSSkyrme eos; //= EOSSkyrme::FreeGas();
  //EOSSkyrme eosInside; 
   
  std::vector<std::unique_ptr<NucleusBase>> nuclei; 
  //nuclei.push_back(StaticNucleus(28, 56, 56.0*8.0/HBC, {}, {}, 9.2*56.0 ).MakeUniquePtr());
  for (int i=50; i<62; i++) {
    //nuclei.push_back(LDNucleus(28, i, eosInside).MakeUniquePtr());
    nuclei.push_back(StaticNucleus(28, i, i*8.0/HBC, {}, {}, i/0.16 ).MakeUniquePtr());
  }
  EOSNSE nseEos(nuclei, eos);
  
  
  double np0, nn0;
  double T = 1.3/HBC;
  double ye = 1.e-6;
  for (nn0 = 0.001; nn0 < 0.05; nn0 *=1.01) { 
    //try {
      np0 = nn0*ye/(1.0-ye);
      std::vector<double> ns = nseEos.GetTotalDensities(EOSData::InputFromTNnNp(T, nn0, np0)); 
      std::cout << ns[0] + ns[1] << " " << nn0 + np0 << " " 
          << ns[1]/(ns[0] + ns[1]) << " " << ye << " " 
          << ns[0] << " " << ns[1] << " " << ns[2] << std::endl; 
      if (ns[2] > 1.0) break; 
      //std::cout << nn0+np0 << " " << T*HBC << " " << ye << " ";
      //std::cout  << eostot.Np() << " " << eostot.Nn() << " "<< std::endl;
    //} catch(std::exception& e ) {
    //  std::cout << nn0+np0  << " " << T*HBC << " " << ye << " fail " << e.what() << std::endl;
    //}
  }
  return 0;
}

