#include <iostream> 
#include <math.h> 
#include <vector>
#include <stdio.h> 

#include "EquationsOfState/EOSData.hpp" 
#include "EquationsOfState/EOSSkyrme.hpp" 
#include "Output/EOSTable.hpp"
#include "Util/Constants.hpp"

#include <H5Cpp.h> 
 
int main() {
   
  EOSSkyrme eos;
  EOSTable tab(eos, 0.1, 100.0, 1.e-8, 1.0, 1.e-4, 0.7, 20, 20, 10); 
  
  auto eosDat = eos.FromNAndS(EOSData::InputFromSNbYe(1.0, 1.e-3, 0.4)); 
  std::cout << eosDat.S() << " " << eosDat.T() * Constants::HBCFmMeV << std::endl;
  
  H5::H5File h5File("temporary_file.h5", H5F_ACC_TRUNC);
     
  tab.WriteToH5(h5File); 
  h5File.close(); 

  if (remove("temporary_file.h5") != 0) return 1;
  std::cout << "Everything worked" << std::endl;
  return 0;
}

