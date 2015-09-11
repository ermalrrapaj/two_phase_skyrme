#include <iostream>
#include <fstream> 
#include <math.h> 
#include <vector> 
#include <utility>
#include <algorithm> 

#include "Util/Constants.hpp"
#include "EquationsOfState/EOSData.hpp" 
#include "EquationsOfState/EOSSkyrme.hpp" 
#include "EquationsOfState/SkyrmeParameters.hpp"    
#include "EquationsOfState/GibbsPhaseConstruct.hpp"
#include "EquationsOfState/EOSSingleNucleus.hpp"

int main() {
  double HBC = Constants::HBCFmMeV;
  
  EOSSkyrme eos = EOSSkyrme::FromErmalSkyrme(ErmalSkyrmeParameters::Ska35S2009);
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

