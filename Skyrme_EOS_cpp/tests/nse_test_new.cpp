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
  
  EOSSkyrme eos; // = EOSSkyrme::FreeGas();
  //EOSSkyrme eosInside; 
   
  std::vector<std::unique_ptr<NucleusBase>> nuclei; 
  //nuclei.push_back(StaticNucleus(28, 56, 56.0*8.0/HBC, {}, {}, 9.2*56.0 ).MakeUniquePtr());
  for (int i=53; i<78; i++) {
    //nuclei.push_back(LDNucleus(28, i, eosInside).MakeUniquePtr());
    //nuclei.push_back(StaticNucleus(27, i, i*7.5/HBC, {}, {}, 1.0*i/0.16 ).MakeUniquePtr());
    nuclei.push_back(StaticNucleus(27, i, i*7.3/HBC, {}, {}, 0.7*i/0.16 ).MakeUniquePtr());
    nuclei.push_back(StaticNucleus(28, i, i*8.0/HBC, {}, {}, 0.7*i/0.16 ).MakeUniquePtr());
    nuclei.push_back(StaticNucleus(29, i, i*7.5/HBC, {}, {}, 0.7*i/0.16 ).MakeUniquePtr());
    nuclei.push_back(StaticNucleus(50, i, i*8.0/HBC, {}, {}, 0.7*i/0.16 ).MakeUniquePtr());
  }
  EOSNSE nseEos(nuclei, eos);
  
  std::cout << "[1] n_{n,t} " << std::endl; 
  std::cout << "[2] n_{p,t} " << std::endl; 
  std::cout << "[3] n_{n,o} " << std::endl; 
  std::cout << "[4] n_{p,o} " << std::endl; 
  std::cout << "[5] sum v_i n_i " << std::endl; 
  std::cout << "[6] E_c (first nucleus) " << std::endl; 
  std::cout << "[7] mu_{n,o} " << std::endl; 
  std::cout << "[8] mu_{p,o} " << std::endl; 
  std::cout << "[9] P_o " << std::endl; 
  double T = 0.2/HBC;
  for (double np0 = 1.e-2; np0<1.e-1; np0 *= 1.2) {
    double nn0 = 1.e-20;
    double delta = 0.1;
    double npo = np0;
    double nno = nn0;
    while (nn0<5.e-1) { 
      try {
        //NSEProperties nse = 
        //    nseEos.GetTotalDensities(EOSData::InputFromTNnNp(T, nn0, np0)); 
        NSEProperties nse = 
            nseEos.GetExteriorProtonDensity(np0, nn0, T);
        double npt = nse.npTot;
        double nnt = nse.nnTot;
        double npe = nse.eosExterior.Np();
        double nne = nse.eosExterior.Nn();
         
        std::cout << nnt << " " << npt << " " 
          << nne << " " << npe << " " 
          << nse.unuc << " " << 
          nuclei[0]->CoulombEnergy(
          nuclei[0]->GetVolume(nse.eosExterior, npt), npe, npt)
          << " " << nse.eosExterior.Mun() << " " << nse.eosExterior.Mup() 
          << " " << nse.eosExterior.P() << std::endl;
        if (fabs(nnt/nno - 1.0) > 0.02 || fabs(npt/npo - 1.0)>0.02) 
            delta *= 0.8;
        if (fabs((nnt+npt)/(npo + nno)-1.0) < 1.e-3) delta /= 0.8;
        delta = std::max(1.e-4, delta);
        delta = std::min(1.e-1, delta);
        npo = npt;
        nno = nnt;
      } catch(std::exception& e ) {
        std::cerr << e.what();
      }
      nn0 *= 1.0 + delta; 
    }
    std::cout << std::endl;
    std::cerr << np0 << std::endl;
  }
  return 0;
}

