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
  void __electron_eos_mod_MOD_electron_eos(double *M, double *n, double *T, 
      double *mu, double *p, double *e, double *s,
      double *dmut, double *dpt, double *det, double *dst,
      double *dmun, double *dpn, double *den, double *dsn);
}

EOSData EOSElectron::FromNAndT(const EOSData& eosIn) {
  double mu, p, e, s;
  double dmun, dpn, den, dsn;
  double dmut, dpt, det, dst;
  double me = Constants::ElectronMassInFm;
  double np = eosIn.Np();
  double T = eosIn.T();
  __electron_eos_mod_MOD_electron_eos(&me, &np, &T, &mu, &p, &e, &s, 
      &dmut, &dpt, &det, &dst, &dmun, &dpn, &den, &dsn);

  EOSData eosOut;
  eosOut.Set("T", eosIn.T()); 
  eosOut.Set("Np", eosIn.Np()); 
  eosOut.Set("Nn", eosIn.Nn()); 
  eosOut.Set("Mue", mu); 
  eosOut.Set("P", p); 
  eosOut.Set("E", e/eosIn.Nb()); 
  eosOut.Set("S", s/eosIn.Nb());
  return eosOut; 

}

