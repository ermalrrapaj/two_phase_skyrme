#include <iostream> 
#include <math.h> 
#include <vector>

#include "EOSData.hpp" 
#include "EOSSkyrme.hpp" 
#include "EOSTable.hpp"

int main() {
   
  EOSSkyrme eos;
  EOSTable tab(eos, 0.1, 100.0, 1.e-8, 1.0, 1.e-4, 0.7, 20, 20, 10); 

  return 0;
}

