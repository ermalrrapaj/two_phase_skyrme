#include <iostream> 
#include <math.h> 
#include <vector>

#include "Util/Constants.hpp"
#include "EquationsOfState/EOSData.hpp" 
#include "EquationsOfState/EOSSkyrme.hpp" 
#include "EquationsOfState/GibbsPhaseConstruct.hpp"
#include "EquationsOfState/SkyrmeParameters.hpp"
#include "EquationsOfState/EOSTestSuite.hpp" 

int main() {
  const double HBC = Constants::HBCFmMeV;
  
/* Test for whether the conversion is right. -> I got the right numbers!	
	double HBC=HBC;
  EOSSkyrme eos2 =EOSSkyrme:: FromErmalSkyrme (0.0, 0.0, -1048.53, -0.0914293, 15488.0, -0.5, 1.002), eos;
  std::cout<<-286.1/HBC<<", "<<-107.1/HBC<<", "<<968.0/HBC<<", "<<0<<", "<<0<<", "<<0<<", "<<2.002<<"\n";
 */
	
  EOSSkyrme eos2 = EOSSkyrme::FromErmalSkyrme(ErmalSkyrmeParameters::Ska35S2009);
  EOSSkyrme eos4 = EOSSkyrme::FromSaturation(SaturationSkyrmeParameters::Ska35S2009);

  EOSData eosIn = EOSData::InputFromTNnNp(5.0/HBC, 0.05, 0.05); 
	EOSData state2 = eos2.FromNAndT(eosIn);
  EOSData state4 = eos4.FromNAndT(eosIn);

  // Check that the two ways of getting the parameters work  
  std::cout << "Ska35S2009 from parameters:\n"; 
  std::cout << "Np = "<<state2.Np() << ", Nn =" << state2.Nn() << ", Ye = " << state2.Ye() << ", T =  "
  << state2.T()*HBC << ", E =  " << state2.E()*HBC << ", S =  " << state2.S()  << ", Mun = "
  << state2.Mun()*HBC << ", Mup =  " << state2.Mup()*HBC  << ", P = " 
  << state2.P()*HBC << " " << std::endl;
  std :: cout << "dMundNn = "<<state2.dMundNn() <<", dMundNp = "<<state2.dMundNp() <<", dMundT = "<<state2.dMundT()
  << ", dMupdNn = "<<state2.dMupdNn() <<", dMupdNp = "<<state2.dMupdNp() <<", dMupdT = "<<state2.dMupdT() 
  << std::endl;
  std :: cout << "dPdNn = "<<state2.dPdNn() <<", dPdNp = "<<state2.dPdNp() <<", dPdT = "<<state2.dPdT()
  << ", dSdNn = "<<state2.dSdNn() <<", dSdNp = "<<state2.dSdNp() <<", dSdT = "<<state2.dMupdT() 
  << std::endl;
  std :: cout << "dEdNn = "<<state2.dEdNn() <<", dEdNp = "<<state2.dEdNp() <<", dEdT = "<<state2.dEdT() 
  << std::endl << std::endl;
 std::cout << "Ska35S2009 from saturation data:\n";
  std::cout << "Np = "<< state4.Np() << ", Nn = " << state4.Nn() << ", Ye =  " << state4.Ye() << ", T = "
  << state4.T()*HBC << ", E = " << state4.E()*HBC << ", S = " << state4.S() << ", Mun = " 
  << state4.Mun()*HBC << ", Mup = " << state4.Mup()*HBC << ", P = " 
  << state4.P()*HBC << " " << std::endl;
  std :: cout << "dMundNn = "<<state4.dMundNn() <<", dMundNp = "<<state4.dMundNp() <<", dMundT = "<<state4.dMundT()
  << ", dMupdNn = "<<state4.dMupdNn() <<", dMupdNp = "<<state4.dMupdNp() <<", dMupdT = "<<state4.dMupdT() 
  << std::endl;
   std :: cout << "dPdNn = "<<state4.dPdNn() <<", dPdNp = "<<state4.dPdNp() <<", dPdT = "<<state4.dPdT()
  << ", dSdNn = "<<state4.dSdNn() <<", dSdNp = "<<state4.dSdNp() <<", dSdT = "<<state4.dMupdT() 
  << std::endl;  
   std :: cout << "dEdNn = "<<state4.dEdNn() <<", dEdNp = "<<state4.dEdNp() <<", dEdT = "<<state4.dEdT() 
  << std::endl << std::endl;

  
  double tol = 1.e-2; // We seem to only get about 1% agreement
  //if (fabs(state4.E()/state2.E() - 1.0)>tol) return 1;
  //if (fabs(state4.S()/state2.S() - 1.0)>tol) return 1;
  //if (fabs(state4.P()/state2.P() - 1.0)>tol) return 1;
  //if (fabs(state4.Mun()/state2.Mun() - 1.0)>tol) return 1;
  //if (fabs(state4.Mup()/state2.Mup() - 1.0)>tol) return 1;
  
  // Check that we can succesfully find a chemical potential 
  EOSData forward = eos2.FromNAndT(eosIn); 
  EOSData reverse = eos2.FromNnMupAndT(
      EOSData::InputFromTNnMup(forward.T(), forward.Nn(), forward.Mup()));
  std::cout << "Mu test : " << reverse.Np() << " " << forward.Np() << std::endl;
  tol = 1.e-6;
  if (fabs(reverse.Np()/forward.Np() - 1.0) > tol) return 1;
  
  // Check that we can succesfully find an entropy 
  forward = eos2.FromNAndT(eosIn); 
  reverse = eos2.FromNAndS(
      EOSData::InputFromSNnNp(forward.S(), forward.Nn(), forward.Np()));
  if (fabs(reverse.T()/forward.T() - 1.0) > tol) return 1;
  
  // Check for consistency of the EoS
  EOSSkyrme free(0, 0, 0, 0, 0, 0, 0); 
  EOSTestSuite test(eos2, 1.e-3, true);  
  std::cout << test.CheckAnalyticDerivatives(2.0/HBC, 0.08, 0.08) << std::endl;
  
  
  return 0;
}

