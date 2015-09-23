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
  void electron_eos_(double *M, double *n, double *T, 
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
  electron_eos_(&me, &np, &T, &mu, &p, &e, &s, &dmut, &dpt, &det, &dst, 
      &dmun, &dpn, &den, &dsn);

  EOSData eosOut;
  eosOut.Set("T", eosIn.T()); 
  eosOut.Set("Np", eosIn.Np()); 
  eosOut.Set("Np", eosIn.Nn()); 
  eosOut.Set("Mue", mu); 
  eosOut.Set("P", p); 
  eosOut.Set("E", e * eosIn.Ye()); 
  eosOut.Set("S", s * eosIn.Ye());
  return eosOut; 

}

