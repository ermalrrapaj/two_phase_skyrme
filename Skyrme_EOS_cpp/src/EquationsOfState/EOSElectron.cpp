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
  eosOut.Set("Mue", mu); 
  eosOut.Set("P", p); 
  eosOut.Set("E", e/eosIn.Nb()); 
  eosOut.Set("S", s/eosIn.Nb());
  eosOut.Set("dPdNp", parr[0][1]);  
  eosOut.Set("dPdT", parr[1][0]);  
  eosOut.Set("dPdNn", 0.0);  

  return eosOut; 

}

