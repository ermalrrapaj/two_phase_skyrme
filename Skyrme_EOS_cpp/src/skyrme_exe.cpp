#include <iostream> 
#include <math.h> 

#include "EOSData.hpp" 
#include "EOSSkyrme.hpp" 
#include "GibbsPhaseConstruct.hpp"

int main() {
	
/* Test for whether the conversion is right. -> I got the right numbers!	
	double HBC=197.3;
  EOSSkyrme eos2 =EOSSkyrme:: FromErmalSkyrme (0.0, 0.0, -1048.53, -0.0914293, 15488.0, -0.5, 1.002), eos;
  std::cout<<-286.1/HBC<<", "<<-107.1/HBC<<", "<<968.0/HBC<<", "<<0<<", "<<0<<", "<<0<<", "<<2.002<<"\n";
 */
  EOSSkyrme eos2=EOSSkyrme:: FromErmalSkyrme (-172.485, 172.087, -1767.71, 0.282732, 12899.2, 0.413266, 0.35), eos;
 // GibbsPhaseConstruct gibbs(eos);

  //for (double lrho = -5.0; lrho < 0.0; lrho += 0.01) {
   // EOSData eosDat = EOSData::InputFromTMunNp(5.0/197.3, 0.0, pow(10.0, lrho)); 
     EOSData eosDat = EOSData::InputFromTNnNp(5.0/197.3, 0.05, 0.05); 
	 EOSData eos3 =eos2.FromNAndT(eosDat);
    //EOSData eosOut = eos.FromNpMunAndT(eosDat);
    
   // std::cout << eosOut.Np() << " " << eosOut.Ye() << " " << eosOut.P() << std::endl;
   std::cout << eos3.Np() << " " << eos3.Nn() << " " << eos3.Ye() <<" "
   << eos3.T()*197.3 << " "<< eos3.E()*197.3 << " "<< eos3.S() 
   <<" "<< eos3.Mun()*197.3<<" "<< eos3.Mup()*197.3<<" "<< eos3.P()*197.3<< " "<< std::endl;
   
  //}


  return 0;
}

