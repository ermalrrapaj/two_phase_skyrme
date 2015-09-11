#include <iostream> 
#include <math.h> 
#include <vector>

#include "EOSData.hpp" 
#include "EOSSkyrme.hpp" 
#include "GibbsPhaseConstruct.hpp"
#include "SkyrmeParameters.hpp"

int main() {

/* Test for whether the conversion is right. -> I got the right numbers!	
	double HBC=197.3;
  EOSSkyrme eos2 =EOSSkyrme:: FromErmalSkyrme (0.0, 0.0, -1048.53, -0.0914293, 15488.0, -0.5, 1.002), eos;
  std::cout<<-286.1/HBC<<", "<<-107.1/HBC<<", "<<968.0/HBC<<", "<<0<<", "<<0<<", "<<0<<", "<<2.002<<"\n";
 */
	
  static constexpr std::array<const double, 7> Ska35S2009 = {{-172.485, 172.087, -1767.71, 0.282732, 12899.2, 0.413266, 0.35}};
 
  EOSSkyrme eos2 = EOSSkyrme::FromErmalSkyrme(Ska35S2009);
  EOSSkyrme eos4 = EOSSkyrme::FromSaturation(SaturationSkyrmeParameters::Ska35S2009);

  EOSData eosIn = EOSData::InputFromTNnNp(5.0/197.3, 0.05, 0.05); 
	EOSData state2 = eos2.FromNAndT(eosIn);
  EOSData state4 = eos4.FromNAndT(eosIn);
    
  std::cout << state2.Np() << " " << state2.Nn() << " " << state2.Ye() << " "
  << state2.T()*197.3 << " " << state2.E()*197.3 << " " << state2.S()  << " "
  << state4.Mun()*197.3 << " " << state4.Mup()*197.3  << " " 
  << state4.P()*197.3 << " " << std::endl;
  
  std::cout << state4.Np() << " " << state4.Nn() << " " << state4.Ye() << " "
  << state4.T()*197.3 << " " << state4.E()*197.3 << " " << state4.S() << " " 
  << state4.Mun()*197.3 << " " << state4.Mup()*197.3 << " " 
  << state4.P()*197.3 << " " << std::endl;

  return 0;
}

