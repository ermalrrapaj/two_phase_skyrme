#include <iostream> 
#include <fstream> 
#include <math.h> 
#include <vector>
#include <memory>
#include <algorithm>
#include <string> 

#include "boost/format.hpp"

#include "Util/Constants.hpp"
#include "EquationsOfState/EOSData.hpp" 
#include "EquationsOfState/EOSSkyrme.hpp" 
#include "EquationsOfState/GibbsPhaseConstruct.hpp"
#include "EquationsOfState/SkyrmeParameters.hpp"
#include "EquationsOfState/NucleusBase.hpp"
#include "EquationsOfState/EOSNSE.hpp"

int main() {
  const double HBC = Constants::HBCFmMeV;
  
  EOSSkyrme eos;// = EOSSkyrme::FreeGas();
  EOSSkyrme eosInside; 
   
  std::vector<std::unique_ptr<NucleusBase>> nuclei; 
  std::vector<std::unique_ptr<NucleusBase>> nucleiStatic; 
  for (int Z=28; Z<=28; Z++) {
    for (int N=27; N<=29; N++) {
      nuclei.push_back(LDNucleus(Z, Z+N, eosInside).MakeUniquePtr());
      nucleiStatic.push_back(
          LDNucleus(Z, Z+N, eosInside).GetStaticNucleus().MakeUniquePtr());
    }
  }
  EOSNSE nseEos(nuclei, eos);
  EOSNSE nseEosStatic(nucleiStatic, eos);
  
  std::ofstream ofile("nse_test_new.out", std::ofstream::out);  
  ofile << "[1] n_{n,t} " << std::endl; 
  ofile << "[2] n_{p,t} " << std::endl; 
  ofile << "[3] n_{n,o} " << std::endl; 
  ofile << "[4] n_{p,o} " << std::endl; 
  ofile << "[5] sum v_i n_i " << std::endl; 
  ofile << "[6] n_{n,t} (static) " << std::endl; 
  ofile << "[7] n_{p,t} (static) " << std::endl; 
  ofile << "[8] n_{n,o} (static) " << std::endl; 
  ofile << "[9] n_{p,o} (static) " << std::endl; 
  ofile << "[10] sum v_i n_i (static) " << std::endl; 
  ofile << "[11] E_c (first nucleus) " << std::endl; 
  ofile << "[12] mu_{n,o} " << std::endl; 
  ofile << "[13] mu_{p,o} " << std::endl; 
  ofile << "[14] P_o " << std::endl; 
  ofile << "[15] BE (first nucleus) " << std::endl; 
  ofile << "[16] BE (first nucleus, static) " << std::endl; 
  ofile << "[17] rho (first nucleus) " << std::endl; 
  
  double T = 1.0/HBC;
  
  for (double np0 = 4.e-2; np0<8.e-2; np0 *= 1.5) {
    std::cout << np0 << std::endl;
    double nn0 = 1.e-14;
    double delta = 0.01;
    double deltaMin = 1.e-8;
    double deltaMax = 1.0;
    NSEProperties nseT = 
      nseEosStatic.GetExteriorProtonDensity(np0, nn0, T)[0];
    double npo = nseT.npTot;
    double nno = nseT.nnTot;
    std::vector<NSEProperties> allPts;
    
    // Find all of the points
    while (nn0<0.08) { 
      try {
        //NSEProperties nse = 
        //    nseEos.GetTotalDensities(EOSData::InputFromTNnNp(T, nn0, np0)); 
        std::vector<NSEProperties> nsev = 
            nseEosStatic.GetExteriorProtonDensity(np0, nn0, T);
        
        allPts.insert(allPts.end(), nsev.begin(), nsev.end());
         
        NSEProperties nse = nsev.back();
        double npt = nse.npTot;
        double nnt = nse.nnTot;
        double npe = nse.eosExterior.Np();
        double nne = nse.eosExterior.Nn();

        if (fabs(nnt/nno - 1.0) > 0.1 && delta>deltaMin) {
          nn0 /= 1.0 + delta; 
          delta *= 0.5;
          std::cerr << " Back step \n"; 
        } else { 
          if (fabs(nnt/nno - 1.0)<0.01) delta /= 0.5;
          nno = nnt;
          npo = npt;
        }
        delta = std::max(deltaMin, delta);
        delta = std::min(deltaMax, delta);
      } catch(std::exception& e ) {
        std::cerr << "Failed at exterior neutron density " << nn0 << std::endl;
        std::cerr << e.what();
      }
      nn0 *= 1.0 + delta; 
    }
    
    // The combined sorting procedure makes plots continous.  The exterior 
    // neutron density is triple valued for fixed total neutron density, so the
    // first sort results in a messy line.  When calculating the EoS, only one 
    // of the three points will have the minimum free energy.
     
    // Sort the points by total neutron density
    std::sort(allPts.begin(), allPts.end(), [](NSEProperties a, NSEProperties b) {
        return b.nnTot > a.nnTot;
    });
    
    // Find the minimum exterior proton density
    auto minEl = std::min_element(allPts.begin(), allPts.end(), [](NSEProperties a, NSEProperties b) {
        return a.eosExterior.Np() < b.eosExterior.Np();
    }); 
    
    // Sort all element above minimum exterior density by exterior density  
    std::sort(minEl, allPts.end(), [](NSEProperties a, NSEProperties b) {
        return b.eosExterior.Nb() > a.eosExterior.Nb();
    });
      
    // Write out the points  
    for (auto& nse : allPts) { 
      double npt = nse.npTot;
      double nnt = nse.nnTot;
      double npe = nse.eosExterior.Np();
      double nne = nse.eosExterior.Nn();
      
      double BePerNuc = 0.0;
      double nucRho = 0.0;
      double Ec = 0.0;
      try { 
        BePerNuc = nuclei[0]->GetBindingEnergy(nse.eosExterior, npt)
          * 197.3 / (double) nuclei[0]->GetA();
        nucRho = nuclei[0]->GetA()/nuclei[0]->GetVolume(nse.eosExterior, npt);
        Ec = nuclei[0]->CoulombEnergy(nuclei[0]->GetVolume(nse.eosExterior, npt), npe, npt);
      } catch(...) {}
      
      double BePerNucStat = nucleiStatic[0]->GetBindingEnergy(nse.eosExterior, npt)
          * 197.3 / (double) nuclei[0]->GetA();
      
      
      std::string frmt;
      for (int ii=0; ii<12; ++ii) frmt.append("% 10.3e "); 
      ofile << boost::format(frmt) 
        % nnt 
        % npt 
        % nne 
        % npe 
        % nse.unuc
        % Ec
        % nse.eosExterior.Mun() 
        % nse.eosExterior.Mup() 
        % nse.eosExterior.P()  
        % BePerNuc 
        % BePerNucStat 
        % nucRho
        << std::endl;
    }
    ofile << std::endl;
  }
  return 0;
}

