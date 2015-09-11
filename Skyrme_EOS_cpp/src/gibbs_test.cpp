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
#include "SkyrmeParameters.hpp"    

int main() {
  double HBC = Constants::HBCFmMeV;
  
	//std::vector<double> Ska35S2009{-172.485, 172.087, -1767.71, 0.282732, 
  //    12899.2, 0.413266, 0.35};
  EOSSkyrme eos = EOSSkyrme::FromErmalSkyrme(ErmalSkyrmeParameters::Ska35S2009);
  EOSSingleNucleus eosSN(eos, false);
  GibbsPhaseConstruct gibbs(eos, false);
  
  // Read in pre-calculated EoS data 
  std::ifstream ifs("eos_data.xml");
  boost::archive::text_iarchive ia(ifs); 
  ia >> gibbs;
  ifs.close();
  
  std::ifstream ifs2("eos_data.xml");
  boost::archive::text_iarchive ib(ifs2); 
  ib >> eosSN;
  ifs2.close(); 
   
  for (double lye = log10(0.3); lye <= log10(0.5); lye += 5.0) { 
    double ye = pow(10.0, lye); 
    std::ofstream outFile; 
    outFile.open("State" + std::to_string(ye) + ".out"); 
    
    for (double lT = log10(0.1); lT < log10(30); lT += 0.25) {
      double T  = pow(10.0, lT)/HBC;
      std::cout << "Point : " << ye << " " << T << std::endl;
      for (double lnb = -6.0; lnb<log10(0.15); lnb += 0.01) {  
        double nb = pow(10.0, lnb); 
        EOSData stateGibbs = 
            gibbs.FromNAndT(EOSData::InputFromTNbYe(T, nb, ye));
        
        EOSData stateBulk = 
            eos.FromNAndT(EOSData::InputFromTNbYe(T, nb, ye));
        
        EOSData stateSN = 
            eosSN.FromNAndT(EOSData::InputFromTNbYe(T, nb, ye));
        
        double u = 1.0; 
        if (stateSN.Phases().size()>1) {
          u = (nb - stateSN.Phases()[0].Nb()) 
              / (stateSN.Phases()[1].Nb() - stateSN.Phases()[0].Nb());
        }
        
        outFile << nb << " " << T << " " << stateGibbs.P() << " " << stateBulk.P(); 
        outFile << " " << stateSN.P() ;
        outFile << " " << stateGibbs.S() << " " << stateBulk.S() ;
        outFile << " " << stateSN.S() << " " << u << std::endl;
      }
      outFile << std::endl;
      std::cout << T*HBC << " done \n";
    }
    outFile.close(); 
  }
  return 0;
   
}

