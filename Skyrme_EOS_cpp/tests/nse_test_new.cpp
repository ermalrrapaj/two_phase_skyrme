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
  for (int Z=20; Z<=100; Z+=10) {
    std::cout << Z << std::endl;
    for (int N=Z-8; N<=Z+100; N+=10) {
      nuclei.push_back(LDNucleus(Z, Z+N, eosInside).MakeUniquePtr());
      nucleiStatic.push_back(
          StaticLDNucleus(Z, Z+N, eosInside).MakeUniquePtr());
      //nucleiStatic.push_back(
      //    StaticLDNucleus(Z, Z+N, eosInside).MakeUniquePtr());
    }
  }
  
  EOSNSE nseEos(nuclei, eos);
  EOSNSE nseEosStatic(nucleiStatic, eos);
  
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
  
  double T = 1.0/HBC;
  
  for (double np0 = 5.e-3; np0<2.e-2; np0 *= 20.0) {
    std::cout << np0 << std::endl;
    double nn0 = 1.e-14;
    double delta = 0.01;
    double deltaMin = 1.e-2;
    double deltaMax = 1.e-1;
    double nno = nn0;
    
    std::vector<std::vector<NSEProperties>> allPts;
    std::vector<NSEProperties> allPtsNSE;
    std::vector<EOSData> alldata; 
    
    // Find all of the points
    while (nn0<5.e-1) { 
      std::cout << nn0;
      try {
        //NSEProperties nse = 
        //    nseEos.GetTotalDensities(EOSData::InputFromTNnNp(T, nn0, np0)); 
        std::vector<NSEProperties> nsev = 
            nseEosStatic.GetExteriorProtonDensity(np0, nn0, T);
        allPts.push_back(nsev);  
        //allPts.insert(allPts.end(), nsev.begin(), nsev.end());
        
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

    std::vector<std::vector<NSEProperties>> branches(100);
    
    std::sort(allPts.begin(), allPts.end(), 
        [](std::vector<NSEProperties> a, std::vector<NSEProperties> b) {
        return b[0].eosExterior.Nn() > a[0].eosExterior.Nn();
    });

    int flip = 0;
    int multiple = 0;
    if (allPts[1].size()>1) multiple = 1;
     
    for (int i=1; i<allPts.size()-1; i++) {
      if (allPts[i].size() != allPts[i-1].size() &&
          allPts[i].size() != allPts[i+1].size()) {
        std::cout << "Skipping point" << std::endl;
        continue;
      }
       
      if (allPts[i].size() != allPts[i-1].size()) flip+= 1; 
      
      // Remove all but the highest neutron density branch  
      //double nntmax = 0.0;
      //int midx = 0; 
      //for (int j=0; j<allPts[i].size(); --j) {
      //  if (allPts[i][j].nnTot>nntmax) {
      //    nntmax = allPts[i][j].nnTot;
      //    midx = j;
      //  }
      //} 
      //branches[flip].push_back(allPts[i][midx]);
      
      // Keep all branches 
      if (allPts[i].size() != allPts[i-1].size()) {
        flip+= allPts[i-1].size();
        if (allPts[i].size()>1) {
          multiple += 1;
        } 
      }
      
      if (multiple%2==1) { 
        for (int j=allPts[i].size()-1; j>=0; --j)  
          branches[flip+j].push_back(allPts[i][j]);
      } else {
        for (int j=0; j<allPts[i].size(); ++j)  
          branches[flip+j].push_back(allPts[i][j]);
      }
    } 

    // The combined sorting procedure makes plots continous.  The exterior 
    // neutron density is triple valued for fixed total neutron density, so the
    // first sort results in a messy line.  When calculating the EoS, only one 
    // of the three points will have the minimum free energy.
     
    // Sort the points by total neutron density
    //std::sort(allPts.begin(), allPts.end(), 
    //    [](NSEProperties a, NSEProperties b) {
    //    return b.nnTot > a.nnTot;
    //});
    // 
    //// Find the minimum exterior proton density
    //auto minEl = std::min_element(allPts.begin(), allPts.end(), 
    //    [](NSEProperties a, NSEProperties b) {
    //    return a.eosExterior.Np() < b.eosExterior.Np();
    //}); 
    //
    //// Sort all element above minimum exterior density by exterior density  
    //std::sort(minEl, allPts.end(), [](NSEProperties a, NSEProperties b) {
    //    return b.eosExterior.Nb() > a.eosExterior.Nb();
    //});
    
	  std::vector<std::vector<NSEProperties>> branchState(branches.size()); 
    for (int i=0; i<branches.size(); ++i) {
      for (auto& pt: branches[i]) 
          branchState[i].push_back(nseEosStatic.GetStateNSEprop(pt));
    } 
	  
    std::string frmt;
    for (int ii=0; ii<15; ++ii) frmt.append("% 10.3e "); 
    auto writeProperties = [&ofile, &eos, frmt, T](const NSEProperties& a)->void{
      auto bulk = eos.FromNAndT(EOSData::InputFromTNnNp(T, a.nnTot, a.npTot)); 
      double fbulk = bulk.E() - T*bulk.S();
      //if (a.F < fbulk) {
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
      //}
    };
 	  
    std::vector<NSEProperties> mixed;
    for (auto& branch : branchState) {
      for (auto& pt : branch) {
        auto bulk = eos.FromNAndT(EOSData::InputFromTNnNp(T, pt.nnTot, pt.npTot)); 
        double fbulk = bulk.E() - T*bulk.S();
        if (fbulk < 0.0) fbulk *= 1.0 + 1.e-4;
        else fbulk *= 1.0 - 1.e-4;
        auto nProp = nseEosStatic.GetTotalDensities(bulk, pt.npTot);
        if (fbulk > pt.F 
            || (fabs((nProp.nnTot + nProp.npTot)/bulk.Nb()-1.0)>1.e-6 
                && bulk.Nb() < 0.04) || bulk.Nn() < 1.e-5)
            mixed.push_back(pt);
        //writeProperties(pt);  
      }
      //ofile << std::endl;
    } 
    
    std::sort(mixed.begin(), mixed.end(), 
        [](NSEProperties a, NSEProperties b) {
        return b.nnTot > a.nnTot;
    });

    // Tack on bulk solution at high density
    for (double nnb = mixed.back().nnTot*1.01; nnb <= 1.0; nnb *=1.01) {
      mixed.push_back(NSEProperties(eos.FromNAndT(
          EOSData::InputFromTNnNp(T, nnb, np0))));
    }
    // Tack on bulk solution at low density  
    for (double nnb = mixed[0].nnTot*0.99; nnb >= 1.e-14; nnb *= 0.99) {
      mixed.push_back(NSEProperties(eos.FromNAndT(
          EOSData::InputFromTNnNp(T, nnb, np0))));
    }

    //Resort 
    std::sort(mixed.begin(), mixed.end(), 
        [](NSEProperties a, NSEProperties b) {
        return b.nnTot > a.nnTot;
    });
    
    // Now calculate on uniform grid
    for (double TT = 0.5; TT < 10.0; TT *=1.5) {
      std::vector<NSEProperties> uniGrid; 
      for (double nnt = 1.e-12; nnt<0.4; nnt *=1.5) {
        std::cout << nnt << " ";
        auto lb = std::lower_bound(mixed.begin(), mixed.end(), nnt,
            [](NSEProperties a, double b) {
              return a.nnTot < b;
        });
        if (lb == mixed.end()) {
          std::cout << std::endl;
          continue;
        }
        
        std::cout << lb->eosExterior.Nn() << " " << lb->nnTot << std::endl;
        try {
        auto gp = nseEosStatic.GetExteriorDensities(
            EOSData::InputFromTNnNp(TT, nnt, np0), lb->eosExterior);
        uniGrid.push_back(gp); 
        } catch (...) { 
          std::cout << "Failed at nnt = " << nnt << " " << TT << std::endl;
        }
      }

      for (auto& pt : uniGrid) 
          writeProperties(pt); 
      ofile << std::endl;
    }

    // Test serialization repeatedly
    //nseEosStatic.SetNSEdata(allPtsNSE);
    
    //std::ofstream ofs("nse.dat"); 
    //boost::archive::binary_oarchive oa(ofs); 
    //oa << nseEosStatic; 
    //ofs.close();

	  //EOSNSE nsedat(nucleiStatic, eos); 
	  //std::ifstream ifs("nse.dat");
	  //boost::archive::binary_iarchive ia(ifs); 
	  //ia >> nsedat;
    //ifs.close();
  
  }  // End of loop over nptot

  return 0;
}

