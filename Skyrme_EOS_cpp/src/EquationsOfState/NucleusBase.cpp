/// \file NucleusBase.cpp
/// \author lroberts
/// \since Sep 14, 2015
///
/// \brief
///
///

#include "EquationsOfState/EOSData.hpp"
#include "EquationsOfState/NucleusBase.hpp"
#include "Util/Constants.hpp"
#include "Util/OneDimensionalRoot.hpp"
#include <iostream> 
#include <math.h> 

double LDNucleus::GetBindingEnergy(const EOSData& eosIn, double ne) const {
  double Z = (double) NucleusBase::mZ; 
  double N = (double) NucleusBase::mN;
  double T = eosIn.T(); 
  auto pFunc = [this, ne, T, Z, N, &eosIn](double v) -> double {
    EOSData eosBulk = mpEos->FromNAndT(
        EOSData::InputFromTNnNp(T, N/v, Z/v)); 
    double Pb = eosBulk.P();
    double Ps = SurfacePressure(v);
    double Pc = CoulombPressure(v, eosIn.Np(), ne);
    return (Pb + Ps + Pc)/eosIn.P() - 1.0;  
  };
  
  OneDimensionalRoot rootFinder(1.e-12); 
  // Just choose an extremely small volume as the lower bound
  double vlo = 1.e-6 * NucleusBase::mA;
  // The maximum volume allowed for the given background electron density 
  // if the nuclear volume were larger, nuclei of this species would 
  // completely fill the space. 
  double vhi = 0.9999 * NucleusBase::mZ / ne;
 
  double v = rootFinder(pFunc, vlo, vhi);
  EOSData eosBulk = mpEos->FromNAndT(
      EOSData::InputFromTNnNp(T, N/v, Z/v)); 
  return -(eosBulk.E() - eosBulk.T()*eosBulk.S()) * NucleusBase::mA 
         - eosIn.P() * v - SurfaceEnergy(v) - CoulombEnergy(v, eosIn.Np(), ne);

}

double LDNucleus::SurfacePressure(double v) const {
  return -2.0 / 3.0 * pow(36.0 * Constants::Pi, 1.0/3.0) * mSigma0 
      * pow(v, -1.0/3.0); 
}

double LDNucleus::SurfaceEnergy(double v) const {
  return pow(36.0 * Constants::Pi, 1.0/3.0) * mSigma0 * pow(v, 2.0/3.0); 
}

namespace {
inline double D(double u) {
  return 1.0 - 1.5 * pow(u, 1.0/3.0) + 0.5*u;
}
inline double DpoD(double u) {
  return -0.5 * pow(u, -2.0/3.0) + 0.5/(D(u) +  1.e-40);
}
}

double LDNucleus::CoulombEnergy(double v, double npo, double ne) const {
  double u = (ne - npo) / ((double)NucleusBase::mZ/v - npo);
  return 3.0 * Constants::ElementaryChargeSquared 
      / (5.0 * pow(3.0 * v / (4.0 * Constants::Pi), 1.0/3.0))
      * pow((double)NucleusBase::mZ - v*npo, 2) * D(u);
}

double LDNucleus::CoulombPressure(double v, double npo, double ne) const {
  double Z = (double) NucleusBase::mZ;
  double denom = 1.0 / (Z/v - npo);
  double u = (ne - npo) * denom;
  double dudlnne = ne / (Z/v - npo);
  return CoulombEnergy(v, npo, ne) / v  
      * (DpoD(u)*u*denom*Z / (v*v) - 2.0*npo/(Z/v - npo) - 1.0/3.0 ); 
}
