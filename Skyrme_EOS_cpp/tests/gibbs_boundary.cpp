#include <iostream>
#include <fstream> 
#include <math.h> 
#include <vector> 
#include <utility>
#include <algorithm> 

#include "Util/Constants.hpp"
#include "EquationsOfState/EOSData.hpp" 
#include "EquationsOfState/EOSSkyrme.hpp" 
#include "EquationsOfState/GibbsPhaseConstruct.hpp"
#include "EquationsOfState/SkyrmeParameters.hpp"    

int main() {
  double HBC = Constants::HBCFmMeV;
  
  //EOSSkyrme eos;
	//std::vector<double> Ska35S2009{-172.485, 172.087, -1767.71, 0.282732, 
  //    12899.2, 0.413266, 0.35};
  EOSSkyrme eos = EOSSkyrme::FromErmalSkyrme(ErmalSkyrmeParameters::Ska35S2009);
  GibbsPhaseConstruct gibbs(eos, false);
  
  // Calculat self bound liquid densities
  std::ofstream selfBoundFile; 
  selfBoundFile.open("self_bound.out");
  selfBoundFile << "# [1] Nb " << std::endl;
  selfBoundFile << "# [2] T (MeV) " << std::endl;
  selfBoundFile << "# [3] Ye " << std::endl;
  selfBoundFile << "# [4] P " << std::endl;
  selfBoundFile << "# [5] Mun " << std::endl;
  selfBoundFile << "# [6] Mup " << std::endl;
  selfBoundFile << "# [7] Nn (low density) " << std::endl;
  selfBoundFile << "# [8] Np (low density) " << std::endl;
  selfBoundFile << "# [9] P (low density) " << std::endl;
  
  std::ofstream outFile; 
  outFile.open("gibbs_boundary.out");
  outFile << "#[1] Nn_1 (fm^-3)" << std::endl; 
  outFile << "#[2] Nn_2 (fm^-3)" << std::endl; 
  outFile << "#[3] Np_1 (fm^-3)" << std::endl; 
  outFile << "#[4] Np_2 (fm^-3)" << std::endl; 
  outFile << "#[5] P (MeV fm^-3)" << std::endl; 
  outFile << "#[6] Mun (MeV fm^-3)" << std::endl; 
  outFile << "#[7] Mup (MeV fm^-3)" << std::endl; 
  outFile << "#[8] T (MeV)" << std::endl; 
  
  double TStart = 0.1/HBC; 
  double TEnd = 0.01/HBC;
  auto phaseBoundFull = gibbs.FindFixedTPhaseBoundary(TStart);
  std::cout << phaseBoundFull.size() << std::endl;
   
  for (double lT = log10(TStart); lT>=log10(TEnd); lT -= 1.0) {
    double T = pow(10.0, lT);
    std::cout << T * HBC << " "; 
    
    //auto selfBound = gibbs.SelfBoundPoints(T);  
    //for (auto &ptpair : selfBound) {
    //  auto pt = ptpair.second;
    //  selfBoundFile << pt.Nb() << " " << pt.T()*HBC << " " << pt.Ye() << " "  
    //      << pt.P() << " " << pt.Mun() << " " << pt.Mup() << " "; 
    //  auto eosPts = 
    //      eos.FromMuAndT(EOSData::InputFromTMunMup(pt.T(), pt.Mun(), pt.Mup()));
    //  if (eosPts.size()>0) {
    //    selfBoundFile << eosPts[0].Np() << " " << eosPts[0].Nn() << " " 
    //        << eosPts[0].P(); 
    //  }
    //  selfBoundFile << std::endl;
    //} 
    //selfBoundFile << std::endl; 
   
    // Calculate phase boundaries
    auto phaseBound = phaseBoundFull; 
    if (T < TStart)  
      phaseBound = gibbs.FindFixedTPhaseBoundary(T, phaseBoundFull); 
    
    std::cout << phaseBound.size() << std::endl; 
    if (phaseBound.size() > 0) {
      // Write out the phase boundary
      for (auto &a : phaseBound) {
        outFile << (a.first).Nn() << " "; 
        outFile << (a.second).Nn() << " "; 
        outFile << (a.first).Np() << " "; 
        outFile << (a.second).Np() << " ";
        outFile << (a.first).P()*HBC << " ";
        outFile << (a.first).Mun()*HBC << " ";
        outFile << (a.first).Mup()*HBC << " ";
        outFile << (a.first).T()*HBC << std::endl;
      }
      outFile << std::endl; 
    } 
  }  
   
  selfBoundFile.close();
  outFile.close(); 
  return 0;
}

