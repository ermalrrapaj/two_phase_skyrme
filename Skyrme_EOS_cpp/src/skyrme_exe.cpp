#include <iostream> 
#include <math.h> 

#include "EOSData.hpp" 
#include "EOSSkyrme.hpp" 

int main() {

  EOSSkyrme eos;
  for (double lrho = -5.0; lrho < 0.0; lrho += 0.01) {
    EOSData eosDat = EOSData::InputFromTMunNp(5.0/197.3, 0.0, pow(10.0, lrho)); 
     
    EOSData eosOut = eos.FromNpMunAndT(eosDat);
    
    std::cout << eosOut.Np() << " " << eosOut.Ye() << " " << eosOut.P() << std::endl;  
  }
  return 0;
}

