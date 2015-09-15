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
    double Pc = CoulombPressure(v, ne);
    return (Pb + Ps + Pc)/eosIn.P() - 1.0;  
  };
  
  OneDimensionalRoot rootFinder(1.e-12); 
  double vlo = 1.e-6 * NucleusBase::mA; 
  double vhi = 1.e4 * NucleusBase::mA;
 
  double v = rootFinder(pFunc, vlo, vhi);
  EOSData eosBulk = mpEos->FromNAndT(
      EOSData::InputFromTNnNp(T, N/v, Z/v)); 
  return -(eosBulk.E() - eosBulk.T()*eosBulk.S()) * NucleusBase::mA 
         - eosIn.P() * v - SurfaceEnergy(v);

}

double LDNucleus::SurfacePressure(double v) const {
  return -2.0 / 3.0 * pow(36.0 * Constants::Pi, 1.0/3.0) * mSigma0 * pow(v, -1.0/3.0); 
}

double LDNucleus::SurfaceEnergy(double v) const {
  return pow(36.0 * Constants::Pi, 1.0/3.0) * mSigma0 * pow(v, 2.0/3.0); 
}

double LDNucleus::CoulombPressure(double v, double ne) const {
  return 0.0 * v * ne; 
}
