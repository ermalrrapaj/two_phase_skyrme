#include <iostream>
#include <fstream> 
#include <math.h> 
#include <vector> 
#include <utility>
#include <algorithm> 

#include "EquationsOfState/EOSData.hpp" 
#include "EquationsOfState/EOSSkyrme.hpp" 
#include "Util/Constants.hpp"
#include "EquationsOfState/GibbsPhaseConstruct.hpp"
#include "EquationsOfState/EOSSingleNucleus.hpp"
#include "EquationsOfState/SkyrmeParameters.hpp"    
#include "EquationsOfState/EOSElectron.hpp"    

#include <boost/archive/text_iarchive.hpp> 
#include <boost/archive/text_oarchive.hpp> 
#include <boost/archive/binary_iarchive.hpp> 
#include <boost/archive/binary_oarchive.hpp> 

int main() {
  const double HBC = Constants::HBCFmMeV;
  const double PI = Constants::Pi;
  
	//std::vector<double> Ska35S2009{-172.485, 172.087, -1767.71, 0.282732, 
  //    12899.2, 0.413266, 0.35};
  double Tmin = 1.e-2;
  double Tmax = 2.e2;
  EOSSkyrme eos = EOSSkyrme::FromErmalSkyrme(ErmalSkyrmeParameters::Ska35S2009);
  EOSSingleNucleus eosSN(eos, false, Tmin, Tmax);
  GibbsPhaseConstruct gibbs(eos, false, Tmin, Tmax);
  EOSElectron eosEl; 
   
  // Read in pre-calculated EoS data 
  std::ifstream ifs("eos_data.arx");
  if (ifs) {
    boost::archive::binary_iarchive ia(ifs); 
    ia >> gibbs;
    ifs.close();
  } else {
    ifs.close();
    gibbs = GibbsPhaseConstruct(eos, true, Tmin, Tmax);
    
    std::ofstream ofs("eos_data.arx"); 
    boost::archive::binary_oarchive oa(ofs); 
    oa << gibbs; 
    ofs.close();  
  }
  
  {
    std::cout << "Trying to read in \n"; 
    std::ifstream ifs2("eos_data.arx");
    boost::archive::binary_iarchive ib(ifs2); 
    ib >> eosSN;
    ifs2.close(); 
  }

  std::ofstream outFile; 
  outFile.open("State.out"); 
  
  for (double lye = log10(0.01); lye <= log10(0.51); lye += (log(0.5)-log10(0.01))/3) { 
    double ye = pow(10.0, lye); 
    outFile << "# [ 1] nb \n"; 
    outFile << "# [ 2] T  \n"; 
    outFile << "# [ 3] Ye  \n"; 
    outFile << "# [ 4] P (Gibbs) \n"; 
    outFile << "# [ 5] P (Bulk)  \n"; 
    outFile << "# [ 6] P (SN)  \n";
    outFile << "# [ 7] P (Electron)  \n";
    outFile << "# [ 8] P (Gamma)  \n";
    outFile << "# [ 9] S (Gibbs)  \n";
    outFile << "# [10] S (Bulk)  \n";
    outFile << "# [11] S (SN)  \n";
    outFile << "# [12] S (Electron)  \n";
    outFile << "# [13] u \n";
     
    for (double lT = log10(0.01); lT < log10(10.0); lT += 10.0) {
      double T  = pow(10.0, lT)/HBC;
      std::cout << "Point : " << ye << " " << T*HBC << std::endl;
      for (double lnb = log10(0.2); lnb>log10(1.e-12); lnb -= 0.1) {  
        double nb = pow(10.0, lnb);
        std::cout << nb << " " << T*HBC << " " << ye << std::endl;
        EOSData stateGibbs = 
            gibbs.FromNAndT(EOSData::InputFromTNbYe(T, nb, ye));
        
        EOSData stateBulk = 
            eos.FromNAndT(EOSData::InputFromTNbYe(T, nb, ye));
        
        EOSData stateSN = 
            eosSN.FromNAndT(EOSData::InputFromTNbYe(T, nb, ye));
        
        EOSData stateEl = 
            eosEl.FromNAndT(EOSData::InputFromTNbYe(T, nb, ye));
        
        double u = 1.0; 
        if (stateSN.Phases().size()>1) {
          u = (nb - stateSN.Phases()[0].Nb()) 
              / (stateSN.Phases()[1].Nb() - stateSN.Phases()[0].Nb());
        }
        
        outFile << nb << " " << T*HBC << " " << ye << " " << stateGibbs.P() << " " << stateBulk.P(); 
        outFile << " " << stateSN.P() << " ";
        outFile << stateEl.P() << " " << PI*PI/45.0*pow(T,4);
        outFile << " " << stateGibbs.S() << " " << stateBulk.S() ;
        outFile << " " << stateSN.S() << " " << stateEl.S() << " " << u << std::endl;
      }
      outFile << std::endl;
      std::cout << T*HBC << " done \n";
    }
    outFile << std::endl;
  }
  outFile.close(); 
  return 0;
   
}

