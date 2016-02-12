#include <iostream> 
#include <fstream> 
#include <ostream> 
#include <istream> 
#include <math.h> 
#include <vector>
#include <memory>
#include <algorithm>
#include <string> 
#include <sys/stat.h> 

#include "boost/format.hpp"

#include "Util/Constants.hpp"
#include "EquationsOfState/EOSData.hpp" 
#include "EquationsOfState/EOSSkyrme.hpp" 
#include "EquationsOfState/GibbsPhaseConstruct.hpp"
#include "EquationsOfState/SkyrmeParameters.hpp"
#include "EquationsOfState/NucleusBase.hpp"
#include "EquationsOfState/LDNucleus.hpp"
#include "EquationsOfState/EOSNSE.hpp"

int main() {
  
  struct stat buffer; // Not sure what this is supposed to do
  
  const double HBC = Constants::HBCFmMeV;
  EOSSkyrme eos; // = EOSSkyrme::FreeGas();
  EOSSkyrme eosInside;

  std::vector<std::unique_ptr<NucleusBase>> nuclei; 
  std::vector<std::unique_ptr<NucleusBase>> nucleiStatic; 
  
  // Create some random nuclei 
  for (int Z=20; Z<=100; Z+=10) {
    std::cout << Z << std::endl;
    for (int N=Z-8; N<=Z+100; N+=10) {
      nuclei.push_back(LDNucleus(Z, Z+N, eosInside).MakeUniquePtr());
      nucleiStatic.push_back(
          StaticLDNucleus(Z, Z+N, eosInside).MakeUniquePtr());
    }
  }
  
  EOSNSE nseEos(nucleiStatic, eos);
  
  std::ofstream ofile("nse_test_new.out", std::ofstream::out);  
  ofile << "#[1] n_{n,t} " << std::endl; 
  ofile << "#[2] n_{p,t} " << std::endl; 
  ofile << "#[3] n_{n,o} " << std::endl; 
  ofile << "#[4] n_{p,o} " << std::endl; 
  ofile << "#[5] F " << std::endl;  
  ofile << "#[6] E " << std::endl;
  ofile << "#[7] P " << std::endl; 
  ofile << "#[8] S " << std::endl; 
  ofile << "#[9] mun_tot " << std::endl; 
  ofile << "#[10] mup_tot " << std::endl; 
  ofile << "#[11] P_{exterior} " << std::endl; 
  ofile << "#[12] unuc " << std::endl; 
  ofile << "#[13] <Ec> [fm^-1]" << std::endl;
  ofile << "#[14] <Be> [fm^-1]" << std::endl;
  ofile << "#[15] Ffree" << std::endl;
  
  std::string frmt;
  for (int ii=0; ii<15; ++ii) frmt.append("% 10.3e "); 
  auto writeProperties = [&ofile, &eos, frmt](
    const NSEProperties& a, double T)->void{
    auto bulk = eos.FromNAndT(EOSData::InputFromTNnNp(T, a.nnTot, a.npTot)); 
    double fbulk = bulk.E() - T*bulk.S();
    ofile << boost::format(frmt) 
      % a.nnTot
      % a.npTot 
      % a.eosExterior.Nn() 
      % a.eosExterior.Np() 
      % a.F
      % a.E
      % a.P
      % a.S
      % a.mun
      % a.mup
      % a.eosExterior.P() 
      % a.unuc
      % a.avgEc
      % a.avgBe
      % fbulk 
      << std::endl;
  };
  
  double T = 1.0/HBC;
  for (double np0 = 5.e-3; np0<2.e-2; np0 *= 20.0) {
    std::cout << np0 << std::endl;

    // Get the array of guesses for given temperature and electron fraction
    auto guessArr = nseEos.GetValidPoints(np0, T, 1.e-14, 0.5); 
 	  
    // Now calculate on uniform nno grid for a range of temperatures
    for (double TT = 0.5/HBC; TT < 10.0/HBC; TT *=2.0) {
      std::cout << "Calculating T = " << TT*HBC << std::endl;
      std::vector<NSEProperties> uniGrid; 
      for (double nnt = 1.e-12; nnt<0.4; nnt *=1.5) {
        std::cout << nnt << " ";
        auto lb = std::lower_bound(guessArr.begin(), guessArr.end(), nnt,
            [](NSEProperties a, double b) {
              return a.nnTot < b;
        });
        if (lb == guessArr.end()) {
          std::cout << std::endl;
          continue;
        }
        lb = lb - 1; 
        std::cout << lb->eosExterior.Nn() << " " << lb->nnTot << std::endl;
        try {
        auto gp = nseEos.GetExteriorDensities(
            EOSData::InputFromTNnNp(TT, nnt, np0), lb->eosExterior);
        uniGrid.push_back(gp); 
        } catch (...) { 
          std::cout << "Failed at nnt = " << nnt << " " << TT << std::endl;
        }
      }

      for (auto& pt : uniGrid) 
          writeProperties(pt, TT); 
      ofile << std::endl;
    }
  }  // End of loop over nptot

  return 0;
}

