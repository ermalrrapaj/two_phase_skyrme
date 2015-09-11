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
  
	std::vector<double> Ska35S2009{-172.485, 172.087, -1767.71, 0.282732, 
      12899.2, 0.413266, 0.35};
  EOSSkyrme eos = EOSSkyrme::FromErmalSkyrme(Ska35S2009);
  //EOSSkyrme eos;
  EOSSingleNucleus gibbs(eos);
  std::ofstream ofs("eos_data.xml"); 
  boost::archive::text_oarchive oa(ofs); 
  oa << gibbs; 
  ofs.close(); 

  //GibbsPhaseConstruct reloadGibbs(eos, false); 
  //
  //std::ifstream ifs("eos_data.xml");
  //boost::archive::text_iarchive ia(ifs); 
  //ia >> reloadGibbs;
  //ifs.close();
  
   
   
  return 0;
}

