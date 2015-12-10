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
  EOSSkyrme eosInside; 
   
  std::vector<std::unique_ptr<NucleusBase>> nuclei; 
  std::vector<std::unique_ptr<NucleusBase>> nucleiStatic; 
  for (int Z=28; Z<=28; Z++) {
    for (int N=28; N<=28; N++) {
      nuclei.push_back(LDNucleus(Z, Z+N, eosInside).MakeUniquePtr());
      nucleiStatic.push_back(
          LDNucleus(Z, Z+N, eosInside).GetStaticNucleus().MakeUniquePtr());
    }
  }
  EOSNSE nseEos(nuclei, eos);
  EOSNSE nseEosStatic(nucleiStatic, eos);
  
  std::cout << "[1] n_{n,t} " << std::endl; 
  std::cout << "[2] n_{p,t} " << std::endl; 
  std::cout << "[3] n_{n,o} " << std::endl; 
  std::cout << "[4] n_{p,o} " << std::endl; 
  std::cout << "[5] sum v_i n_i " << std::endl; 
  std::cout << "[6] n_{n,t} (static) " << std::endl; 
  std::cout << "[7] n_{p,t} (static) " << std::endl; 
  std::cout << "[8] n_{n,o} (static) " << std::endl; 
  std::cout << "[9] n_{p,o} (static) " << std::endl; 
  std::cout << "[10] sum v_i n_i (static) " << std::endl; 
  std::cout << "[11] E_c (first nucleus) " << std::endl; 
  std::cout << "[12] mu_{n,o} " << std::endl; 
  std::cout << "[13] mu_{p,o} " << std::endl; 
  std::cout << "[14] P_o " << std::endl; 
  std::cout << "[15] BE (first nucleus) " << std::endl; 
  std::cout << "[16] BE (first nucleus, static) " << std::endl; 
  std::cout << "[17] rho (first nucleus) " << std::endl; 
  
  double T = 1.5/HBC;
  
  for (double np0 = 1.e-2; np0<1.1e-2; np0 *= 1.2) {
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
        NSEProperties nseStat = 
            nseEosStatic.GetExteriorProtonDensity(np0, nn0, T);
        double npt = nse.npTot;
        double nnt = nse.nnTot;
        double npe = nse.eosExterior.Np();
        double nne = nse.eosExterior.Nn();
         
        std::cout << nnt << " " << npt << " " 
          << nne << " " << npe << " " 
          << nse.unuc << " " 
          << nseStat.nnTot << " " << nseStat.npTot << " " 
          << nseStat.eosExterior.Nn() << " " << nseStat.eosExterior.Np() << " " 
          << nseStat.unuc << " " << 
          nuclei[0]->CoulombEnergy(
          nuclei[0]->GetVolume(nse.eosExterior, npt), npe, npt)
          << " " << nse.eosExterior.Mun() << " " << nse.eosExterior.Mup() 
          << " " << nse.eosExterior.P()  
          << " " << nuclei[0]->GetBindingEnergy(nse.eosExterior, npt)
          * 197.3 / (double) nuclei[0]->GetA()
          << " " << nucleiStatic[0]->GetBindingEnergy(nse.eosExterior, npt)
          * 197.3 / (double) nucleiStatic[0]->GetA()
          << " " << nuclei[0]->GetA()/nuclei[0]->GetVolume(nse.eosExterior, npt)
          << std::endl;
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

