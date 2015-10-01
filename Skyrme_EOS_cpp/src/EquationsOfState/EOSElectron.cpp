/// \file EOSElectron.hpp
/// \author lroberts
/// \since Sep 23, 2015
///
/// \brief
///
///

#include "Util/Constants.hpp"
#include "EquationsOfState/EOSElectron.hpp"

extern"C" {
  void __fermi_dirac_MOD_two_fermion_eos(double *g, double *M, double *T, 
      double *n, double mu[][2], double p[][2], double e[][2], double s[][2]);
}

EOSData EOSElectron::FromNAndT(const EOSData& eosIn) {
  double g, mu, p, e, s, n;
  double dmun, dpn, den, dsn;
  double dmut, dpt, det, dst;
  double me = Constants::ElectronMassInFm;
  double np = eosIn.Np();
  double T = eosIn.T();
  
  double muarr[2][2]; 
  double parr[2][2]; 
  double earr[2][2]; 
  double sarr[2][2]; 
  g = 2.0;
  __fermi_dirac_MOD_two_fermion_eos(&g, &me, &T, &np, muarr, parr, earr, sarr);
  mu = muarr[0][0]; 
  p = parr[0][0]; 
  e = earr[0][0]; 
  s = sarr[0][0]; 
   
  EOSData eosOut;
  eosOut.Set("T", eosIn.T()); 
  eosOut.Set("Np", eosIn.Np()); 
  eosOut.Set("Nn", eosIn.Nn()); 
  
  eosOut.Set("P", p); 

  eosOut.SetDNp("P", parr[0][1]);  
  eosOut.SetDT( "P", parr[1][0]);  
  eosOut.SetDNn("P", 0.0);  
  
  double nb = eosIn.Nb();
  eosOut.Set("S", sarr[0][0]/nb);
  eosOut.SetDNp("S", sarr[0][1]/nb - sarr[0][0]/(nb*nb));  
  eosOut.SetDT( "S", sarr[1][0]/nb);  
  eosOut.SetDNn("S", -sarr[0][0]/(nb*nb));  
  
  eosOut.Set("E", e/nb); 
  eosOut.SetDNp("E", earr[0][1]/nb - earr[0][0]/(nb*nb));  
  eosOut.SetDT( "E", earr[1][0]/nb);  
  eosOut.SetDNn("E", -earr[0][0]/(nb*nb));  

  
  eosOut.Set("Mue", muarr[0][0]); 
  eosOut.SetDNp("Mue", muarr[0][1]);  
  eosOut.SetDT( "Mue", muarr[1][0]);  
  eosOut.SetDNn("Mue", 0.0);  
  
  
  return eosOut; 

}

