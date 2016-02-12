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

static double HBC = Constants::HBCFmMeV; 
static double MNUC = 938.918/HBC;

std::vector<double> StaticNucleus::CoulombEnergy(double v, double npo, 
    double ne) const {
  double Z = (double) NucleusBase::mZ; 
  double u = v*(ne - npo) / (Z - v*npo);
  double D = 0.5*u - 1.5*pow(u,1.0/3.0);
  double udDduoD = (0.5*u- 0.5*pow(u,1.0/3.0))/(D + 1.e-40);
  double dludv = 1.0/ (v + 1.e-40) + npo / (Z - v*npo);
  double dludne = v / (ne - npo + 1.e-40);
  double dludnpo = (u - 1.0) / (ne - npo + 1.e-40);
  double vn1o3 = pow(v, -1.0/3.0);
  double fac = 3.0 * Constants::ElementaryChargeSquared 
      / (5.0 * pow(3.0 / (4.0 * Constants::Pi), 1.0/3.0));

  double Ec      = fac * vn1o3 * D * pow(Z - npo*v, 2);
  double dEcdv   = Ec * (udDduoD * dludv - 1.0/(3.0 * v) - 2.0 * npo / (Z-npo*v));
  double dEcdnpo = Ec * (udDduoD * dludnpo - 2.0 * v / (Z-npo*v)); 
  double dEcdne  = Ec *  udDduoD * dludne; 
  
  return {Ec, dEcdv, dEcdnpo, dEcdne};
}

double StaticNucleus::FreeEnergy(const EOSData& eosIn, double ne, double ni) const {
	double T  = eosIn.T();
	double BE = GetBindingEnergy(eosIn, ne);
	double nQ = pow(MNUC*T/2/Constants::Pi,1.5); 
	double Fk = T*log((ni+1.e-100)/nQ/pow(NucleusBase::mA,1.5)) - T;
	return Fk - BE;
}

double StaticNucleus::Entropy (const EOSData& eosIn, double ne, double ni) const {
	return 0.0;
}

double StaticNucleus::NucleusPressure (const EOSData& eosIn, double ne, double uo) const{
  auto Ec = CoulombEnergy(GetVolume(eosIn, ne), eosIn.Np(), ne);
  return eosIn.T() + Ec[3] + eosIn.Np()/(uo+1.e-80)*Ec[2]; 
}

double StaticNucleus::Nucleusmup (const EOSData& eosIn, double ne, double uo, double ni) const {
  auto Ec = CoulombEnergy(GetVolume(eosIn, ne), eosIn.Np(), ne);
	return ni/uo * Ec[2];
}

double StaticNucleus::Nucleusmun (const EOSData& eosIn, double ne, double uo, double ni) const {
	return 0.0;
}

