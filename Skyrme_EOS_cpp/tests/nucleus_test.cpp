#include <iostream> 
#include <math.h> 
#include <vector>

#include "Util/Constants.hpp"
#include "EquationsOfState/EOSData.hpp" 
#include "EquationsOfState/EOSSkyrme.hpp" 
#include "EquationsOfState/GibbsPhaseConstruct.hpp"
#include "EquationsOfState/SkyrmeParameters.hpp"
#include "EquationsOfState/NucleusBase.hpp"
#include "EquationsOfState/LDNucleus.hpp"

int main() {
  const double HBC = Constants::HBCFmMeV;
  
  EOSSkyrme eos;

  double Zi = 28;
  double Ai = 56;
  double Ni = Ai-Zi;
  LDNucleus nuc(Zi, Ai, eos); 
  //StaticLDNucleus nucS(Zi, Ai, eos); 
  

  return 0;
  for (double ln = -3.0; ln < log10(0.2); ln += 0.01) { 
	double T = 0.1/HBC;
	double nB = pow(10.0, ln);
	double ye = 0.01;
	double nn = nb*(1.0-ye);
	double np = nb*ye;
    EOSData extState = eos.FromNAndT(EOSData::InputFromTNbYe(T, nB, ye));
  //  np0 = exState.
   // double den = Ni * np0 - nn0 * Zi;
//	double ni = (-nn0 * np + nn * np0)/den;
//	double uo = (Ni * np - nn * Zi)/den;
//    std::cout << pow(10.0, ln) << " " 
 //   << HBC/56.0 * nuc.GetBindingEnergy(extState, 0.1, uo) << std::endl; 
  }
  return 0;
}

