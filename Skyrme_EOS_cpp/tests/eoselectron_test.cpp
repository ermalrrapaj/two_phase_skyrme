#include <iostream> 
#include <math.h> 
#include <vector>

#include "Util/Constants.hpp"
#include "EquationsOfState/EOSData.hpp" 
#include "EquationsOfState/EOSElectron.hpp" 

int main() {
  const double HBC = Constants::HBCFmMeV;
  
  EOSElectron eose; 
  
  // Check that we agree with the high temperature, zero electron fraction limit 
  EOSData eosIn = EOSData::InputFromTNnNp(10000.0/HBC, 1.e-3, 0.0); 
	EOSData state = eose.FromNAndT(eosIn); 
  std::cout << state.Np() << " " << state.Nn() << " " << state.Ye() << " "
  << state.T()*HBC << " " << state.E()*HBC << " " << state.S()  << " "
  << state.Mue()*HBC << " " << state.P()*HBC << " " << std::endl;
  double err = 2.0 * 0.191909 * pow(eosIn.T(), 4)/state.P() - 1.0;
  std::cout << err << std::endl; 
  if (fabs(err)>1.e-6) return 1;

  // Check that we agree in the degenerate, relativistic limit
  eosIn = EOSData::InputFromTNnNp(1.e-10/HBC, 1.e-2, 1.e1);
	state = eose.FromNAndT(eosIn); 
  double mudeg = pow(3.0 * Constants::Pi * Constants::Pi * state.Np(), 1.0/3.0);
  double Pdeg = 1.0/(12.0 * Constants::Pi * Constants::Pi) * pow(mudeg, 4.0);
  
  err =  state.Mue()/mudeg - 1.0; 
  std::cout << state.Mue() << " " << mudeg << " " << err << std::endl;
  if (fabs(err)>1.e-6) return 1;
  err = state.P()/Pdeg - 1.0;
  std::cout << state.P() << " " << Pdeg << " " << err << std::endl;
  if (fabs(err)>1.e-6) return 1;
  
  // Check that we agree in the non-degenerate, non-relativistic limit
  eosIn = EOSData::InputFromTNnNp(1.e-3/HBC, 1.e-2, 1.e-18);
	state = eose.FromNAndT(eosIn); 
  double nQ = pow(Constants::ElectronMassInFm * state.T() 
      / (2.0 * Constants::Pi),1.5);
  double muST = Constants::ElectronMassInFm  + state.T() * log(state.Np()/(2.0*nQ));
  double pST = state.Np() * state.T();
  
  std::cout << state.Np() << " " << nQ <<std::endl; 
  err =  state.Mue()/muST - 1.0; 
  std::cout << state.Mue() << " " << muST << " " << err << std::endl;
  if (fabs(err)>1.e-5) return 1;
  err = state.P()/pST - 1.0;
  std::cout << state.P() << " " << pST << " " << err << std::endl;
  if (fabs(err)>1.e-5) return 1;
  
  //double Pdeg = pow(3.0 * Constants::Pi * Constants::Pi, 2.0/3.0) 
  //    / (5.0 * Constants::ElectronMassInFm) * pow(state.Np(), 5.0/3.0);  

  return 0;
}

