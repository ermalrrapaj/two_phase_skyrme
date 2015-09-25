#include <iostream> 
#include <math.h> 
#include <vector>

#include "EquationsOfState/EOSData.hpp" 

int main() {
  EOSData eos = EOSData::InputFromTNbYe(1.0, 2.0, 3.0);
  
  std::cout << eos.T() << " "; 
  
  eos.Set("S", 4.0);
  std::cout << eos.S() << " "; 
  std::cout << eos.Get("Np") << std::endl; 
   
  try {
    double P = eos.P(); 
  } catch(...) {
    return 0;
  }
  return 1;
}
