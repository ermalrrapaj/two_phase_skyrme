#include <iostream> 
#include <math.h> 
#include <vector>

#include "Util/Constants.hpp"
#include "EquationsOfState/EOSData.hpp" 
#include "EquationsOfState/EOSElectron.hpp" 
#include "EquationsOfState/EOSTestSuite.hpp"

extern"C" {
  void __fermi_dirac_MOD_dfermi(double *dk, double *denom, double *deta, 
      double *dtheta, double fd[4][5]); 
}


int FermiDiracEtaTest(double eta, double theta) {

  // First, test that we are taking derivatives of the fermi function correctly
  double dk,denom; 
  double fdeta,fdtheta,fdeta2,fdtheta2,fdeta3,fdtheta3;
  double fdetadtheta,fdeta2dtheta,fdetadtheta2;
  dk = 3.0; 
  denom = 1.0; 
  double delta = 1.e-4; 
   
  std::vector<std::vector<double>> fd(5, std::vector<double>(5, 0.0));
  
  for (int i=0; i<5; i++) {
    for (int j=0; j<5; j++) {
      double deta = eta * (1.0 + delta * (double)(i-2));
      double dtheta = theta * (1.0 + delta * (double)(j-2));
      double ft;
      double farr[4][5];
      __fermi_dirac_MOD_dfermi(&dk, &denom, &deta, &dtheta, farr); 
      fd[i][j] = farr[0][0];
    }
  } 

  double ft; 
  double farr[4][5];
  __fermi_dirac_MOD_dfermi(&dk, &denom, &eta, &theta, farr); 
  ft = farr[0][0]; 
  fdeta = farr[0][1]; 
  fdeta2 = farr[0][2]; 
  fdeta3 = farr[0][3]; 
  fdtheta = farr[1][0];   
  fdtheta2 = farr[2][0];   
  fdtheta3 = farr[3][0];   
  fdetadtheta = farr[1][1]; 
  fdeta2dtheta = farr[1][2]; 
  fdetadtheta2 = farr[2][1]; 
   
  int ierr = 0; 
  {
    double h = delta * eta; 
    double nd  = (fd[3][2] - fd[1][2])/(2.0*h); 
    double nd2 = (fd[3][2] - 2.0*fd[2][2] + fd[1][2])/pow(h, 2); 
    double nd3 = 0.5*(fd[4][2] - 2.0*fd[3][2] + 2.0*fd[1][2] - fd[0][2])/pow(h, 3); 
    
    std::cout << fdeta  << " " << nd/fdeta   - 1.0 << std::endl; 
    std::cout << fdeta2 << " " << nd2/fdeta2 - 1.0 << std::endl; 
    std::cout << fdeta3 << " " << nd3/fdeta3 - 1.0 << std::endl; 
    if (fabs(nd /fdeta  - 1.0)  > 1.e-7) ++ierr;
    if (fabs(nd2/fdeta2 - 1.0) > 1.e-5) ++ierr;
    if (fabs(nd3/fdeta3 - 1.0) > 1.e-3) ++ierr;
  }
   
  {
    double h = delta * theta; 
    double nd  = (fd[2][3] - fd[2][1])/(2.0*h); 
    double nd2 = (fd[2][3] - 2.0*fd[2][2] + fd[2][1])/pow(h, 2); 
    double nd3 = 0.5*(fd[2][4] - 2.0*fd[2][3] + 2.0*fd[2][1] - fd[2][0])/pow(h, 3); 
    
    std::cout << fdtheta  << " " << nd /fdtheta  - 1.0 << std::endl; 
    std::cout << fdtheta2 << " " << nd2/fdtheta2 - 1.0 << std::endl; 
    std::cout << fdtheta3 << " " << nd3/fdtheta3 - 1.0 << std::endl; 
    if (fabs(nd /fdtheta  - 1.0)  > 1.e-7) ++ierr;
    if (fabs(nd2/fdtheta2 - 1.0) > 1.e-5) ++ierr;
    if (fabs(nd3/fdtheta3 - 1.0) > 1.e-3) ++ierr;
  }
   
  double heta = delta*eta; 
  double htheta = delta*theta; 
  
  double nd11 = (fd[3][3] - fd[1][3] - fd[3][1] + fd[1][1])/(4.0*heta*htheta);
  std::cout << fdetadtheta << " " << nd11/fdetadtheta - 1.0 << std::endl; 
  if (fabs(nd11/fdetadtheta - 1.0) > 1.e-3) ++ierr;
   
  double nd21 = (fd[3][3] - 2.0*fd[2][3] + fd[1][3] - fd[3][1] + 2.0*fd[2][1] - fd[1][1])/(2.0*heta*htheta*heta);
  std::cout << fdeta2dtheta << " " << nd21/fdeta2dtheta - 1.0 << std::endl; 
  if (fabs(nd21/fdeta2dtheta - 1.0) > 1.e-3) ++ierr;
  
  double nd12 = (fd[3][3] - 2.0*fd[3][2] + fd[3][1] - fd[1][3] + 2.0*fd[1][2] - fd[1][1])/(2.0*htheta*htheta*heta);
  std::cout << fdetadtheta2 << " " << nd12/fdetadtheta2 - 1.0 << std::endl; 
  if (fabs(nd12/fdetadtheta2 - 1.0) > 1.e-3) ++ierr;
  
  return ierr; 
}

int CheckLimits() {
  

  const double HBC = Constants::HBCFmMeV;
  EOSElectron eose; 
  
  // Check that we agree with the high temperature, zero electron fraction limit 
  EOSData eosIn = EOSData::InputFromTNnNp(10000.0/HBC, 1.e-3, 0.0); 
	EOSData state = eose.FromNAndT(eosIn); 
   
  double Prad = 2.0 * 0.191909 * pow(eosIn.T(), 4);
  double err = Prad/state.P() - 1.0;
  std::cout << "rel, nondeg : " << state.P() << " " << Prad << " "<< err << std::endl; 
  if (fabs(err)>1.e-6) return 1;
  if (fabs(state.Mue()/eosIn.T())>1.e-10) return 1;
  err = 4.0*Prad/eosIn.T()/state.dPdT() - 1.0;
  std::cout << "rel, nondeg : " << state.dPdT() << " " << 4.0*Prad/eosIn.T() << " "<< err << std::endl; 
  if (fabs(err)>1.e-6) return 1;

  // Check that we agree in the degenerate, relativistic limit
  eosIn = EOSData::InputFromTNnNp(1.e-10/HBC, 1.e-2, 1.e1);
	state = eose.FromNAndT(eosIn); 
  
  double mudeg = pow(3.0 * Constants::Pi * Constants::Pi * state.Np(), 1.0/3.0);
  double Pdeg = 1.0/(12.0 * Constants::Pi * Constants::Pi) * pow(mudeg, 4.0);
  double dmuDeg = 1.0/3.0*mudeg/state.Np(); 
  double dPdeg = 4.0*Pdeg/mudeg*dmuDeg; 
   
  err =  state.Mue()/mudeg - 1.0; 
  std::cout << "rel, deg : " << state.Mue() << " " << mudeg << " " << err << std::endl;
  if (fabs(err)>1.e-6) return 1;
  err = state.P()/Pdeg - 1.0;
  std::cout << "rel, deg : " << state.P() << " " << Pdeg << " " << err << std::endl;
  if (fabs(err)>1.e-6) return 1;
  err = state.dPdNp()/dPdeg - 1.0;
  std::cout << "rel, deg : " << state.dPdNp() << " " << dPdeg << " " << err << std::endl;
  if (fabs(err)>1.e-6) return 1;
  
  // Check that we agree in the non-degenerate, non-relativistic limit
  eosIn = EOSData::InputFromTNnNp(1.e-3/HBC, 1.e-2, 1.e-18);
	state = eose.FromNAndT(eosIn); 
  
  double nQ = pow(Constants::ElectronMassInFm * state.T() 
      / (2.0 * Constants::Pi),1.5);
  double muST = Constants::ElectronMassInFm  + state.T() * log(state.Np()/(2.0*nQ));
  double pST = state.Np() * state.T();
  
  err =  state.Mue()/muST - 1.0; 
  std::cout << "Nonrel, nondeg : " << state.Mue() << " " << muST << " " << err << std::endl;
  if (fabs(err)>1.e-5) return 1;
  err = state.P()/pST - 1.0;
  std::cout << "Nonrel, nondeg : "<< state.P() << " " << pST << " " << err << std::endl;
  if (fabs(err)>1.e-5) return 1;
  err = state.dPdNp()/state.T() - 1.0;
  std::cout << "Nonrel, nondeg : "<< state.dPdNp() << " " << state.T() << " " << err << std::endl;
  if (fabs(err)>1.e-5) return 1;
  err = state.dPdT()/state.Np() - 1.0;
  std::cout << "Nonrel, nondeg : "<< state.dPdT() << " " << state.Np() << " " << err << std::endl;
  if (fabs(err)>1.e-5) return 1;
  
  /// non-relativistic, degenerate limit test 
  eosIn = EOSData::InputFromTNnNp(1.e-10/HBC, 1.e-2, 1.e-16);
	state = eose.FromNAndT(eosIn); 
  
  mudeg = Constants::ElectronMassInFm 
      + pow(3.0*Constants::Pi*Constants::Pi*state.Np(), 2.0/3.0) / 
      (2.0*Constants::ElectronMassInFm);
  dmuDeg = 2.0/3.0*mudeg/state.Np(); 
  Pdeg = pow(3.0 * Constants::Pi*Constants::Pi, 2.0/3.0)
      / (5.0*Constants::ElectronMassInFm) 
      * pow(state.Np(), 5.0/3.0);
  dPdeg = 5.0/3.0/state.Np()*Pdeg;
  
  err =  state.Mue()/mudeg - 1.0; 
  std::cout << "Nonrel, deg : " << state.Mue() << " " << mudeg << " " << err << std::endl;
  if (fabs(err)>1.e-6) return 1;
  err = state.P()/Pdeg - 1.0;
  std::cout << "Nonrel, deg : "<< state.P() << " " << Pdeg << " " << err << std::endl;
  if (fabs(err)>1.e-4) return 1;
  err = state.dPdNp()/dPdeg - 1.0;
  std::cout << "Nonrel, deg : "<< state.dPdNp() << " " << dPdeg << " " << err << std::endl;
  if (fabs(err)>1.e-4) return 1;
  
  return 0;
}

int main() {
  // Do tests that are specific to the electron Eos
  int ierr = FermiDiracEtaTest(10.0, 1.2);
  ierr += CheckLimits(); 
  ierr = 0; 
  // Do general EOS tests over temperature density grid 
  const double HBC = Constants::HBCFmMeV;
  EOSElectron eos; 
  EOSTestSuite test(eos, 1.e-5, true);
  for (double lT = -2.0; lT < 4.0; lT += 0.5) {
    for (double ln = -3.0; ln < 2.0; ln += 0.5) {
      double T = pow(10.0,lT)/HBC;
      double ne = pow(10.0,ln);
      
      ierr += test.CheckAnalyticDerivatives(T, 0.1, ne);   
      ierr += test.CheckThermodynamicConsistency(T, 0.1, ne);   
      if (ierr>0) return ierr;  
    }
  }  
  
  return ierr; 

}
