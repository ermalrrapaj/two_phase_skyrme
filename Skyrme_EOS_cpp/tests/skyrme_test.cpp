#include <iostream> 
#include <math.h> 
#include <vector>

#include "Util/Constants.hpp"
#include "EquationsOfState/EOSData.hpp" 
#include "EquationsOfState/EOSSkyrme.hpp" 
#include "EquationsOfState/GibbsPhaseConstruct.hpp"
#include "EquationsOfState/SkyrmeParameters.hpp"

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
  std::cout << state2.Np() << " " << state2.Nn() << " " << state2.Ye() << " "
  << state2.T()*HBC << " " << state2.E()*HBC << " " << state2.S()  << " "
  << state4.Mun()*HBC << " " << state4.Mup()*HBC  << " " 
  << state4.P()*HBC << " " << std::endl;
  
  std::cout << state4.Np() << " " << state4.Nn() << " " << state4.Ye() << " "
  << state4.T()*HBC << " " << state4.E()*HBC << " " << state4.S() << " " 
  << state4.Mun()*HBC << " " << state4.Mup()*HBC << " " 
  << state4.P()*HBC << " " << std::endl;
  
  const double tol = 1.e-2; // We seem to only get about 1% agreement
  if (fabs(state4.E()/state2.E() - 1.0)>tol) return 1;
  if (fabs(state4.S()/state2.S() - 1.0)>tol) return 1;
  if (fabs(state4.P()/state2.P() - 1.0)>tol) return 1;
  if (fabs(state4.Mun()/state2.Mun() - 1.0)>tol) return 1;
  if (fabs(state4.Mup()/state2.Mup() - 1.0)>tol) return 1;
  
  // Check that we can succesfully find a chemical potential 
  EOSData forward = eos2.FromNAndT(eosIn); 
  EOSData reverse = eos2.FromNnMupAndT(
      EOSData::InputFromTNnMup(forward.T(), forward.Nn(), forward.Mup()));
  std::cout << "Mu test : " << reverse.Np() << " " << forward.Np() << std::endl;
  if (fabs(reverse.Np()/forward.Np() - 1.0) > tol) return 1;
  
  // Check that we can succesfully find an entropy 
   
  return 0;
}

