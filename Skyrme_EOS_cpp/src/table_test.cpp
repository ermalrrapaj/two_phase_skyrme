#include <iostream> 
#include <math.h> 
#include <vector>

#include "EOSData.hpp" 
#include "EOSSkyrme.hpp" 
#include "EOSTable.hpp"

#include <H5Cpp.h> 
 
int main() {
   
  EOSSkyrme eos;
  EOSTable tab(eos, 0.1, 100.0, 1.e-8, 1.0, 1.e-4, 0.7, 20, 20, 10); 
  
  H5::H5File h5File("file.h5", H5F_ACC_TRUNC);
     
  tab.WriteToH5(h5File); 

  return 0;
}

