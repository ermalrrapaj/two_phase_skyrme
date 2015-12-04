#include <iostream> 
#include <math.h> 
#include <vector>
#include <memory>
#include <algorithm> 

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
  for (int i=53; i<65; i++) {
    //nuclei.push_back(LDNucleus(28, i, eosInside).MakeUniquePtr());
    //nuclei.push_back(StaticNucleus(27, i, i*7.5/HBC, {}, {}, 1.0*i/0.16 ).MakeUniquePtr());
    nuclei.push_back(StaticNucleus(28, i, i*8.0/HBC, {}, {}, 1.0*i/0.16 ).MakeUniquePtr());
    //nuclei.push_back(StaticNucleus(29, i, i*7.5/HBC, {}, {}, 1.0*i/0.16 ).MakeUniquePtr());
  }
  EOSNSE nseEos(nuclei, eos);
  
  
  double T = 1.0/HBC;
  for (double np0 = 1.e-11; np0<1.e-1; np0 *= 2.0) {
    double nn0 = 1.e-11;
    double delta = 0.01;
    double npo = np0;
    double nno = nn0;
    while (nn0<1.e-1) { 
      try {
        std::vector<double> ns = 
            nseEos.GetTotalDensities(EOSData::InputFromTNnNp(T, nn0, np0)); 
        auto eosExt = eos.FromNAndT(EOSData::InputFromTNnNp(T, nn0, np0)); 
        if (ns[2] < 1.0) { 
            std::cout << ns[0] << " " << ns[1] << " " 
            << nn0 << " " << np0 << " " << ns[2] << std::endl;
          if (fabs(ns[0]/nno - 1.0) > 0.02 || fabs(ns[1]/npo - 1.0)>0.02) delta *= 0.8;
          if (fabs((ns[0] + ns[1])/(npo + nno)-1.0) < 1.e-3) delta /= 0.8;
          delta = std::max(1.e-4, delta);
          delta = std::min(1.e-1, delta);
          npo = ns[1];
          nno = ns[0];
        }
      } catch(std::exception& e ) {}
      nn0 *= 1.0 + delta; 
    }
    std::cout << std::endl;
    std::cerr << np0 << std::endl;
  }
  return 0;
}

