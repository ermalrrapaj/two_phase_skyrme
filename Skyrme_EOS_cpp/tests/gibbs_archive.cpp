#include <iostream>
#include <string>
#include <fstream> 
#include <math.h> 
#include <vector> 
#include <utility>
#include <algorithm>
#include <sys/stat.h> 

#include "Util/Constants.hpp"
#include "EquationsOfState/EOSData.hpp" 
#include "EquationsOfState/EOSSkyrme.hpp" 
#include "EquationsOfState/SkyrmeParameters.hpp"    
#include "EquationsOfState/GibbsPhaseConstruct.hpp"
#include "EquationsOfState/EOSSingleNucleus.hpp"

int main() 
{
  struct stat buffer;   
  double HBC = Constants::HBCFmMeV;
  std:: array<std::string,5> Skyrmenames09 = {"Ska35s2009.xml","Ska25s2009.xml","SKT109.xml","SKT209.xml","SKT309.xml"};
  std:: array<std::string,5> Skyrmenames1 = {"Ska35s201.xml","Ska25s201.xml","SKT11.xml","SKT21.xml","SKT31.xml"};
  std:: string eos_dataName_ermal("eos_data_"),eos_dataName_saturation("eos_data_sat_"), eos_dataName; 
  for (int counter = 0; counter < 5; counter++)
  {
	eos_dataName=eos_dataName_ermal+Skyrmenames09[counter];
	if((stat (eos_dataName.c_str(), &buffer) != 0))
	{
			EOSSkyrme eos = EOSSkyrme::FromErmalSkyrme(ErmalSkyrmeParameters::ErmalAllSkyrme09[counter]);
			EOSSingleNucleus gibbs(eos);
			std::ofstream ofs(eos_dataName.c_str()); 
			boost::archive::text_oarchive oa(ofs); 
			oa << gibbs; 
			ofs.close();
	}
	eos_dataName=eos_dataName_ermal+Skyrmenames1[counter];
	if((stat (eos_dataName.c_str(), &buffer) != 0))
	{
			EOSSkyrme eos = EOSSkyrme::FromErmalSkyrme(ErmalSkyrmeParameters::ErmalAllSkyrme1[counter]);
			EOSSingleNucleus gibbs(eos);
			std::ofstream ofs(eos_dataName.c_str()); 
			boost::archive::text_oarchive oa(ofs); 
			oa << gibbs; 
			ofs.close();
	}
	eos_dataName=eos_dataName_saturation+Skyrmenames09[counter];
	if((stat (eos_dataName.c_str(), &buffer) != 0))
	{
			EOSSkyrme eos = EOSSkyrme::FromSaturation(SaturationSkyrmeParameters::SaturationSkyrme09[counter]);
			EOSSingleNucleus gibbs(eos);
			std::ofstream ofs(eos_dataName.c_str()); 
			boost::archive::text_oarchive oa(ofs); 
			oa << gibbs; 
			ofs.close();
	}
	eos_dataName=eos_dataName_saturation+Skyrmenames1[counter];
	if((stat (eos_dataName.c_str(), &buffer) != 0))
	{
			EOSSkyrme eos = EOSSkyrme::FromSaturation(SaturationSkyrmeParameters::SaturationSkyrme1[counter]);
			EOSSingleNucleus gibbs(eos);
			std::ofstream ofs(eos_dataName.c_str()); 
			boost::archive::text_oarchive oa(ofs); 
			oa << gibbs; 
			ofs.close();
	}
  }
  //GibbsPhaseConstruct reloadGibbs(eos, false); 
  //
  //std::ifstream ifs("eos_data.xml");
  //boost::archive::text_iarchive ia(ifs); 
  //ia >> reloadGibbs;
  //ifs.close();
  
   
   
  return 0;
}

