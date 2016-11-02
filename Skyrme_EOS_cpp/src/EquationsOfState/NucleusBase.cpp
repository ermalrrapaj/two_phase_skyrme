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
  double u = v*(ne - npo) / (Z - v*npo + 1.e-80);
  double D = 0.5*u - 1.5*pow(u,1.0/3.0) + 1.0;
  double DprimeD = 0.5*(1.0- pow(u,-2.0/3.0))/(D + 1.e-80);
  double dudlogne = ne*v / (Z - v*npo + 1.e-80);
  double dudlogv = u*Z/(Z-v*npo+ 1.e-80);
  double dudlognpo = u*npo*(1.0/(npo - ne + 1.e-80) + v/(Z-v*npo+1.e-80));
  double vn1o3 = pow(v, -1.0/3.0);
  double fac = 3.0 * Constants::ElementaryChargeSquared 
      / (5.0 * pow(3.0 / (4.0 * Constants::Pi), 1.0/3.0));

  double Ec      = fac * vn1o3 * D * pow(Z - npo*v, 2);
  double dEcdv   = Ec/v * (DprimeD*dudlogv - 1.0/3.0 - 2.0*npo*v/(Z-v*npo+1.e-80));
  double dEcdlognpo = Ec * (DprimeD * dudlognpo - 2.0 * v*npo/(Z-npo*v+1.e-80)); 
  double dEcdlogne  = Ec *  DprimeD * dudlogne; 
  
  std::vector<double> result(4);
  result[0] = Ec;
  result[1] = dEcdv;
  result[2] = dEcdlognpo;
  result[3] = dEcdlogne;
  
  return result;
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

///
/// Calculate the pressure from a single nucleus,
/// to get the total pressure contribution -> ni * NucleusPressure
///
double StaticNucleus::NucleusPressure (const EOSData& eosIn, double ne, double uo) const{
  auto Ec = CoulombEnergy(GetVolume(eosIn, ne), eosIn.Np(), ne);
  return eosIn.T() + Ec.at(3) + eosIn.Np()/(uo+1.e-80)*Ec[2]; 
}

double StaticNucleus::Nucleusmup (const EOSData& eosIn, double ne, double uo, double ni) const {
  auto Ec = CoulombEnergy(GetVolume(eosIn, ne), eosIn.Np(), ne);
	return ni/uo * Ec[2];
}

double StaticNucleus::Nucleusmun (const EOSData& eosIn, double ne, double uo, double ni) const {
	return 0.0;
}

double StaticNucleus::GetDensity(const EOSData& eosIn, double ne, double uo, 
    double v) const {
  double nQ = pow(Constants::NeutronMassInFm 
      * eosIn.T() / (2.0 * Constants::Pi), 1.5);
  double aa = (GetN()*eosIn.Mun() + GetZ()*eosIn.Mup() 
      + GetBindingEnergy(eosIn, ne, v) - v*eosIn.P())/eosIn.T();
  return std::min(nQ * pow((double) GetA(), 1.5) * exp(aa), 1.e200);
}
