#include <iostream> 
#include <fstream> 
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
#include "EquationsOfState/EOSNSE.hpp"

int main() {
  struct stat buffer;
  const double HBC = Constants::HBCFmMeV;
  EOSSkyrme eos; //= EOSSkyrme::FreeGas();
  EOSSkyrme eosInside;
  std::vector<std::unique_ptr<NucleusBase>> nuclei; 
  std::vector<std::unique_ptr<NucleusBase>> nucleiStatic; 
  for (int Z=28; Z<=28; Z++) {
    for (int N=20; N<=40; N++) {
      nuclei.push_back(LDNucleus(Z, Z+N, eosInside).MakeUniquePtr());
      nucleiStatic.push_back(
          LDNucleus(Z, Z+N, eosInside).GetStaticNucleus().MakeUniquePtr());
    }
  }
  //if((stat ("nse.dat", &buffer) != 0)){
  EOSNSE nseEos(nuclei, eos);
  EOSNSE nseEosStatic(nuclei, eos);
  /*
  std::ofstream ofile("nse_test_new.out", std::ofstream::out);  
  ofile << "[1] n_{n,t}_1 " << std::endl; 
  ofile << "[2] n_{p,t}_1 " << std::endl; 
  ofile << "[3] n_{n,o} " << std::endl; 
  ofile << "[4] n_{p,o} " << std::endl; 
  ofile << "[5] sum v_i n_i " << std::endl;  
  ofile << "[6] E_c (first nucleus) " << std::endl; 
  ofile << "[7] P_o " << std::endl;
  ofile << "[8] BE (first nucleus) " << std::endl; 
  ofile << "[9] BE (first nucleus, static) " << std::endl; 
  ofile << "[10] rho (first nucleus) " << std::endl;
  ofile << "[11] n_{n,t}_2 " << std::endl; 
  ofile << "[12] n_{p,t}_2 " << std::endl; 
  ofile << "[13] mun_tot " << std::endl; 
  ofile << "[14] mup_tot " << std::endl; 
  ofile << "[15] P_tot " << std::endl; 
  ofile << "[16] E_tot " << std::endl; 
  ofile << "[17] S_tot " << std::endl;
   */
  double T = 1.0/HBC;
  std::vector<NSEProperties> allPtsNSE;
  
 // for (double np0 = 4.e-2; np0<8.e-2; np0 *= 1.5) {
   double np0=2.e-2;
    std::cout << np0 << std::endl;
    double nn0 = 1.e-12;
    double delta = 0.01;
    double deltaMin = 1.e-1;
    double deltaMax = 1.0;
    NSEProperties nseT = 
      nseEosStatic.GetExteriorProtonDensity(np0, nn0, T)[0];
    double npo = nseT.npTot;
    double nno = nseT.nnTot;
    std::vector<NSEProperties> allPts;
    std::vector<EOSData> alldata; 
    
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
     double npt[2],nnt[2],npe,nne,mun,mup,BePerNuc,BePerNucStat,nucRho,Ec,Unuc,P[2],Etot,Stot;  
    // Write out the points  
    for (auto& nse : allPts) { 
	  alldata.push_back(nseEosStatic.GetState(nse));}
	
	for (auto& nse : allPts) { 
	  allPtsNSE.push_back(nseEosStatic.GetStateNSEprop(nse));}  
	 
 	  
 /* 
  for (int i =0; i< allPts.size(); i++){
	  npt[0] = allPts[i].npTot;
      nnt[0] = allPts[i].nnTot;
      P[0] = allPts[i].eosExterior.P();
      npe = allPts[i].eosExterior.Np();
      nne = allPts[i].eosExterior.Nn();
      Unuc = allPts[i].unuc;
      
      BePerNuc = 0.0;
      nucRho = 0.0;
      Ec = 0.0;
      try { 
        BePerNuc = nuclei[0]->GetBindingEnergy(allPts[i].eosExterior, npt[0])
          * 197.3 / (double) nuclei[0]->GetA();
        nucRho = nuclei[0]->GetA()/nuclei[0]->GetVolume(allPts[i].eosExterior, npt[0]);
        Ec = nuclei[0]->CoulombEnergy(nuclei[0]->GetVolume(allPts[i].eosExterior, npt[0]), npe, npt[0]);
      } catch(...) {}
      
      BePerNucStat = nucleiStatic[0]->GetBindingEnergy(allPts[i].eosExterior, npt[0])
          * 197.3 / (double) nuclei[0]->GetA();
	  
	  
      npt[1] = alldata[i].Np();
      nnt[1] = alldata[i].Nn();
      P[1] = alldata[i].P();
      Etot = alldata[i].E();
      Stot = alldata[i].S();
      mun = alldata[i].Mun();
      mup = alldata[i].Mup();
      
      std::string frmt;
      for (int ii=0; ii<17; ++ii) frmt.append("% 10.3e "); 
      ofile << boost::format(frmt) 
        % nnt[0] 
        % npt[0] 
        % nne 
        % npe 
        % Unuc
        % Ec
        % P[0]  
        % BePerNuc 
        % BePerNucStat 
        % nucRho
        % nnt[1]
        % npt[1]
        % mun
        % mup
        % P[1]
        % Etot
        % Stot
        << std::endl;
    }
    ofile << std::endl;
  */ //}
  nseEosStatic.SeTNSEdata(allPtsNSE);
 // std::ofstream ofile("nse_test_new.out", std::ofstream::out);
  //EOSSingleNucleus gibbs(eos);
  std::ofstream ofs("nse.dat"); 
  boost::archive::text_oarchive oa(ofs); 
  oa << nseEosStatic; 
  ofs.close();
  std::ofstream ofile("nse_test_new_nb_F.out", std::ofstream::out);
  for (auto& nse : allPtsNSE) { 
	  ofile << nse.eosExterior.Nn() << "\t" << nse.nnTot <<"\t" 
	  		<<nse.F << "\t" << nse.P <<std::endl;// allPtsNSE.push_back(nseEosStatic.GetStateNSEprop(nse));
	  }
	  ofile.close();
 /* }
  else {
	  EOSNSE nsedat(nucleiStatic, eos); 
	  std::ifstream ifs("nse.dat");
	  boost::archive::text_iarchive ia(ifs); 
	  ia >> nsedat;
      ifs.close();
  
  }*/
  
  return 0;
}

