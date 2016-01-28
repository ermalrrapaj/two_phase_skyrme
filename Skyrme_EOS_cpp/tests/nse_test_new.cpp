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
  struct stat buffer;
  const double HBC = Constants::HBCFmMeV;
  EOSSkyrme eos; // = EOSSkyrme::FreeGas();
  EOSSkyrme eosInside;
  std::vector<std::unique_ptr<NucleusBase>> nuclei; 
  std::vector<std::unique_ptr<NucleusBase>> nucleiStatic; 
  for (int Z=28; Z<=28; Z++) {
    for (int N=20; N<=40; N++) {
      nuclei.push_back(LDNucleus(Z, Z+N, eosInside).MakeUniquePtr());
      nucleiStatic.push_back(
          LDNucleus(Z, Z+N, eosInside).GetStaticNucleus().MakeUniquePtr());
      //nucleiStatic.push_back(
      //    StaticLDNucleus(Z, Z+N, eosInside).MakeUniquePtr());
    }
  }
  
  EOSNSE nseEos(nuclei, eos);
  EOSNSE nseEosStatic(nucleiStatic, eos);
  
  std::ofstream ofile("nse_test_new.out", std::ofstream::out);  
  ofile << "[1] n_{n,t} " << std::endl; 
  ofile << "[2] n_{p,t} " << std::endl; 
  ofile << "[3] n_{n,o} " << std::endl; 
  ofile << "[4] n_{p,o} " << std::endl; 
  ofile << "[5] F " << std::endl;  
  ofile << "[6] E " << std::endl;
  ofile << "[7] P " << std::endl; 
  ofile << "[8] S " << std::endl; 
  ofile << "[9] mun_tot " << std::endl; 
  ofile << "[10] mup_tot " << std::endl; 
  ofile << "[11] P_{exterior} " << std::endl; 
  ofile << "[12] unuc " << std::endl; 
  ofile << "[13] <Ec> [fm^-1]" << std::endl;
  ofile << "[14] <Be> [fm^-1]" << std::endl;
  
  double T = 1.0/HBC;
  
  for (double np0 = 2.e-2; np0<8.e-2; np0 *= 20.0) {
    std::cout << np0 << std::endl;
    double nn0 = 1.e-14;
    double delta = 0.01;
    double deltaMin = 1.e-1;
    double deltaMax = 1.0;
    double nno = nn0;
    std::vector<NSEProperties> allPts;
    std::vector<NSEProperties> allPtsNSE;
    std::vector<EOSData> alldata; 
    
    // Find all of the points
    while (nn0<0.1) { 
      std::cout << nn0;
      try {
        //NSEProperties nse = 
        //    nseEos.GetTotalDensities(EOSData::InputFromTNnNp(T, nn0, np0)); 
        std::vector<NSEProperties> nsev = 
            nseEosStatic.GetExteriorProtonDensity(np0, nn0, T);
        
        allPts.insert(allPts.end(), nsev.begin(), nsev.end());
        
        if (nsev.size() > 0) { 
          std::cout << " " << nsev.size() << std::endl;
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
          }
        } else {
          std::cout << std::endl;
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
    std::sort(allPts.begin(), allPts.end(), 
        [](NSEProperties a, NSEProperties b) {
        return b.nnTot > a.nnTot;
    });
    
    // Find the minimum exterior proton density
    auto minEl = std::min_element(allPts.begin(), allPts.end(), 
        [](NSEProperties a, NSEProperties b) {
        return a.eosExterior.Np() < b.eosExterior.Np();
    }); 
    
    // Sort all element above minimum exterior density by exterior density  
    std::sort(minEl, allPts.end(), [](NSEProperties a, NSEProperties b) {
        return b.eosExterior.Nb() > a.eosExterior.Nb();
    });
    
	
	  for (auto& nse : allPts) 
        allPtsNSE.push_back(nseEosStatic.GetStateNSEprop(nse));
	 
 	  
    std::string frmt;
    for (int ii=0; ii<14; ++ii) frmt.append("% 10.3e "); 
    for (int i =0; i< allPtsNSE.size(); i++){
      ofile << boost::format(frmt) 
        % allPtsNSE[i].nnTot
        % allPtsNSE[i].npTot 
        % allPtsNSE[i].eosExterior.Nn() 
        % allPtsNSE[i].eosExterior.Np() 
        % allPtsNSE[i].F
        % allPtsNSE[i].E
        % allPtsNSE[i].P
        % allPtsNSE[i].S
        % allPtsNSE[i].mun
        % allPtsNSE[i].mup
        % allPtsNSE[i].eosExterior.P() 
        % allPtsNSE[i].unuc
        % allPtsNSE[i].avgEc
        % allPtsNSE[i].avgBe 
        << std::endl;
    }
    ofile << std::endl;

    // Test serialization repeatedly
    nseEosStatic.SeTNSEdata(allPtsNSE);
    
    std::ofstream ofs("nse.dat"); 
    boost::archive::binary_oarchive oa(ofs); 
    oa << nseEosStatic; 
    ofs.close();

	  EOSNSE nsedat(nucleiStatic, eos); 
	  std::ifstream ifs("nse.dat");
	  boost::archive::binary_iarchive ia(ifs); 
	  ia >> nsedat;
    ifs.close();
  
  }  // End of loop over nptot

  return 0;
}

