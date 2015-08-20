#include <iostream> 
#include <math.h> 

#include "EOSData.hpp" 
#include "EOSSkyrme.hpp" 
#include "GibbsPhaseConstruct.hpp"

int main() {
  EOSSkyrme eos;
  GibbsPhaseConstruct gibbs(eos);
  
  std::vector<EOSData> outDat = gibbs.FindPhasePoint(2.0/197.3, 0.0); 
  std::cout << std::endl;
  std::cout << "Nb:  " << outDat[0].Nb()  << " " << outDat[1].Nb()  << std::endl;
  std::cout << "Ye:  " << outDat[0].Ye()  << " " << outDat[1].Ye()  << std::endl;
  std::cout << "P:   " << outDat[0].P()   << " " << outDat[1].P()   << std::endl;
  std::cout << "Mun: " << outDat[0].Mun() << " " << outDat[1].Mun() << std::endl;
  std::cout << "Mup: " << outDat[0].Mup() << " " << outDat[1].Mup() << std::endl;
  return 0;
}

