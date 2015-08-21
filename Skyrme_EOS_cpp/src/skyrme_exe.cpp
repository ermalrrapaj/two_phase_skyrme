#include <iostream> 
#include <math.h> 
#include <vector>

#include "EOSData.hpp" 
#include "EOSSkyrme.hpp" 
#include "GibbsPhaseConstruct.hpp"

int main() {
	std::vector<double> Ska35S2009{-172.485, 172.087, -1767.71, 0.282732, 12899.2, 0.413266, 0.35};
	std::vector<double> Ska25S2009{-137.191, 153.932, -2181.28, 0.340552, 14605.8, 0.46161, 0.25};
	std::vector<double> SKT109{-112.324, 142.467, -1810.72, 0.28223, 12863., 0.392585, 1.0/3.0};
	std::vector<double> SKT209{-113.857, 143.999, -1807.87, 0.267778, 12802.4, 0.366144, 1.0/3.0};
	std::vector<double> SKT309{-124.432, 148.492, -1812.16, 0.288584, 12906.6, 0.416129, 1.0/3.0};
	
	std::vector<double> Ska35S201{-2.41114, -0.507978, -1767.92, 0.247025, 12910.2, 0.220377, 0.35};
	std::vector<double> Ska25S201{7.88383, -0.431231, -2182.96, 0.293491, 14651.7, 0.304538, 0.25};
	std::vector<double> SKT11{25.175, -12.2601, -1815.65, 0.24924, 12984.8, 0.217553, 1.0/3.0};
	std::vector<double> SKT21{-28.3961, 16.7196, -1822.27, 0.258133, 13155.8, 0.24668, 1.0/3.0};
	std::vector<double> SKT31{4.51143, 0.0220502, -1816.13, 0.256993, 13025.9, 0.247513, 1.0/3.0};
	
	std::vector<double> data_Ska35S2009{0.157915, -16.0445, 1.00038, 239.851, 32.1508, 53.5764, -137.329};
	std::vector<double> data_Ska25S2009{0.160492, -16.0374, 0.984062, 220.085, 32.123, 50.8728, -134.829};
	std::vector<double> data_SKT109{0.158366, -16.0629, 0.972029, 237.932, 32.6151, 55.3149, -135.069};
	std::vector<double> data_SKT209{0.159957, -16.2195, 0.971756, 240.127, 33.1405, 57.9679, -130.311};
	std::vector<double> data_SKT309{0.158289, -16.0429, 0.977556, 237.422, 32.0233, 52.6175, -138.493};
	
	std::vector<double> data_Ska35S201{0.160944, -16.0461, 0.992821, 219.76, 32.0304, 46.3568, -164.742};
	std::vector<double> data_Ska25S201{0.158052, -16.0526, 1.00279, 239.875, 32.2029, 50.6062, -161.435};
	std::vector<double> data_SKT11{0.158519, -16.065, 0.987809, 237.264, 32.5327, 51.5294, -160.419};
	std::vector<double> data_SKT21{0.158706, -16.0573, 1.0113, 236.173, 32.3098, 50.4438, -160.611};
	std::vector<double> data_SKT31{0.158687, -16.049, 0.995682, 236.752, 31.9537, 48.8985, -163.314};
/* Test for whether the conversion is right. -> I got the right numbers!	
	double HBC=197.3;
  EOSSkyrme eos2 =EOSSkyrme:: FromErmalSkyrme (0.0, 0.0, -1048.53, -0.0914293, 15488.0, -0.5, 1.002), eos;
  std::cout<<-286.1/HBC<<", "<<-107.1/HBC<<", "<<968.0/HBC<<", "<<0<<", "<<0<<", "<<0<<", "<<2.002<<"\n";
 */
 
  EOSSkyrme eos2 = EOSSkyrme:: FromErmalSkyrme(Ska35S2009),
  eos4 = EOSSkyrme:: FromSaturation(data_Ska35S2009);
 // GibbsPhaseConstruct gibbs(eos);

  //for (double lrho = -5.0; lrho < 0.0; lrho += 0.01) {
   // EOSData eosDat = EOSData::InputFromTMunNp(5.0/197.3, 0.0, pow(10.0, lrho)); 
     EOSData eosDat = EOSData::InputFromTNnNp(5.0/197.3, 0.05, 0.05); 
	 EOSData eos3 =eos2.FromNAndT(eosDat),eos5 =eos4.FromNAndT(eosDat);
   // EOSData eosOut = eos.FromNpMunAndT(eosDat);
    
   // std::cout << eosOut.Np() << " " << eosOut.Ye() << " " << eosOut.P() << std::endl;
   std::cout << eos3.Np() << " " << eos3.Nn() << " " << eos3.Ye() <<" "
   << eos3.T()*197.3 << " "<< eos3.E()*197.3 << " "<< eos3.S() 
   <<" "<< eos5.Mun()*197.3<<" "<< eos5.Mup()*197.3<<" "<< eos5.P()*197.3<< " "<< std::endl;
   std::cout << eos5.Np() << " " << eos5.Nn() << " " << eos5.Ye() <<" "
   << eos5.T()*197.3 << " "<< eos5.E()*197.3 << " "<< eos5.S() 
   <<" "<< eos5.Mun()*197.3<<" "<< eos5.Mup()*197.3<<" "<< eos5.P()*197.3<< " "<< std::endl;

  //}


  return 0;
}

