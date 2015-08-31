#include <iostream>
#include <fstream> 
#include <math.h> 
#include <vector> 
#include <utility>
#include <algorithm> 

#include "EOSData.hpp" 
#include "EOSSkyrme.hpp" 
#include "Constants.hpp"
#include "GibbsPhaseConstruct.hpp"
#include "EOSSingleNucleus.hpp"

int main() {
  double HBC = Constants::HBCFmMeV;
  
  EOSSkyrme eos;
  GibbsPhaseConstruct gibbs(eos, false);
  
  std::ifstream ifs("eos_data.xml");
  boost::archive::text_iarchive ia(ifs); 
  ia >> gibbs;
  ifs.close();
  
  for (double lye = -2.0; lye <= log10(0.5); lye += 0.1) { 
    double ye = pow(10.0, lye); 
    std::ofstream outFile; 
    outFile.open("State" + std::to_string(ye) + ".out"); 
    
    for (double lT = log10(0.1); lT < log10(30.0); lT += 0.05) {
      double T  = pow(10.0, lT)/HBC;
      for (double lnb = -8.0; lnb<0.5; lnb += 0.005) {  
        double nb = pow(10.0, lnb); 
        EOSData stateGibbs = 
            gibbs.FromNAndT(EOSData::InputFromTNbYe(T, nb, ye));
        
        EOSData stateBulk = 
            eos.FromNAndT(EOSData::InputFromTNbYe(T, nb, ye));
        
        outFile << nb << " " << T << " " << stateGibbs.P() << " " << stateBulk.P(); 
        outFile << " " << stateGibbs.S() << " " << stateBulk.S() << std::endl;
      }
      outFile << std::endl;
      outFile << T*HBC << " done \n";
    }
    outFile.close(); 
  }
  return 0;
   
  for (double T = 0.01/HBC; T <= 0.1/HBC; T += 0.01/HBC) {
    auto phaseBound = gibbs.FindFixedTPhaseBoundary(T); 
    
    std::cerr << T * HBC << " " << phaseBound.size() << std::endl; 
    if (phaseBound.size() > 0) {
      // Write out the phase boundary
      std::cout << "#[1] Nn_1 (fm^-3)" << std::endl; 
      std::cout << "#[2] Nn_2 (fm^-3)" << std::endl; 
      std::cout << "#[3] Np_1 (fm^-3)" << std::endl; 
      std::cout << "#[4] Np_2 (fm^-3)" << std::endl; 
      std::cout << "#[5] P (MeV fm^-3)" << std::endl; 
      std::cout << "#[6] Mun (MeV fm^-3)" << std::endl; 
      std::cout << "#[7] Mup (MeV fm^-3)" << std::endl; 
      std::cout << "#[8] T (MeV)" << std::endl; 
      
      for (auto &a : phaseBound) {
        std::cout << (a.first).Nn() << " "; 
        std::cout << (a.second).Nn() << " "; 
        std::cout << (a.first).Np() << " "; 
        std::cout << (a.second).Np() << " ";
        std::cout << (a.first).P()*HBC << " ";
        std::cout << (a.first).Mun()*HBC << " ";
        std::cout << (a.first).Mup()*HBC << " ";
        std::cout << (a.first).T()*HBC << std::endl;
      }
      
       
      //std::reverse(phaseBound.begin(), phaseBound.end());  
      //
      //for (auto &a : phaseBound) {
      //  std::cout << (a.second).Nn() << " "; 
      //  std::cout << (a.first).Nn() << " "; 
      //  std::cout << (a.second).Np() << " ";
      //  std::cout << (a.first).Np() << " "; 
      //  std::cout << (a.first).P()*HBC << " ";
      //  std::cout << (a.first).Mun()*HBC << " ";
      //  std::cout << (a.first).Mup()*HBC << " ";
      //  std::cout << (a.first).T()*HBC << std::endl;
      //}
      //std::cout << (phaseBound.back().first).Nn() << " "; 
      //std::cout << (phaseBound.back().second).Nn() << " "; 
      //std::cout << (phaseBound.back().first).Np() << " ";
      //std::cout << (phaseBound.back().second).Np() << " "; 
      //std::cout << (phaseBound.back().first).P()*HBC << " ";
      //std::cout << (phaseBound.back().first).Mun()*HBC << " ";
      //std::cout << (phaseBound.back().first).Mup()*HBC << " ";
      //std::cout << (phaseBound.back().first).T()*HBC << std::endl;
      
      std::cout << std::endl; 
    } 
  }   
    
  return 0;
}

