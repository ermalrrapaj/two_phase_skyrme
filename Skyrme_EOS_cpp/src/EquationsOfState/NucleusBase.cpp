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


double NucleusBase::CoulombEnergy(double v, double npo, double ne) const {
  double u = v*(ne - npo) / ((double)NucleusBase::mZ - v*npo);
  double Ec = 3.0 * Constants::ElementaryChargeSquared 
      / (5.0 * pow(3.0 * v / (4.0 * Constants::Pi), 1.0/3.0))
      * pow((double)NucleusBase::mZ - npo*v, 2) * (0.5*u - 1.5*pow(u,1.0/3.0)) 
      * pow(ne/(double)NucleusBase::mZ,1.0/3.0);
  return Ec;
}

double NucleusBase::CoulombPressureExternal(double v, double npo, double ne) const {
  double u = v*(ne - npo) / ((double)NucleusBase::mZ - v*npo);
  return 3.0 * Constants::ElementaryChargeSquared 
      / (5.0 * pow(3.0 / (4.0 * Constants::Pi), 1.0/3.0))
      * pow((double)NucleusBase::mZ - v*npo, 2) * pow(u,-1.0/3.0)/3.0
      * v/((double)NucleusBase::mZ - v*npo);
}

double NucleusBase::FreeEnergy(const EOSData& eosIn, double ne, double ni) const {
	double T  = eosIn.T();
	double BE = GetBindingEnergy(eosIn, ne);
	double nQ = pow(MNUC*T/2/Constants::Pi,1.5); 
	double Fk = T*log((ni+1.e-100)/nQ/pow(NucleusBase::mA,1.5)) - T;
	return Fk - BE;
}

double NucleusBase::Entropy (const EOSData& eosIn, double ne, double ni) const {
	return 0.0;
}

double NucleusBase::NucleusPressure (const EOSData& eosIn, double ne, double uo) const{
	double npo = eosIn.Np();
	double T = eosIn.T();
	double v = GetVolume(eosIn, ne);
	double u = v*(ne - npo) / ((double)NucleusBase::mZ - v*npo +  1.e-100);
	double EC = CoulombEnergy(v, npo, ne);
  double DodD = (0.5*pow(u,5.0/3.0) - 1.5*u)/(0.5*pow(u, 2.0/3.0) - 0.5 + 1.e-40); 
  double dDdu = 0.5 - 0.5*pow(u, -2.0/3.0);
  
  double dECdlogne = EC*DodD * v * ne / ((double)NucleusBase::mZ - v*npo +  1.e-40);
  double dECdlognpo = -EC*DodD* v * npo / ((double)NucleusBase::mZ - v*npo +  1.e-40)
    * (1.0 - u);

  return T + dECdlogne + uo*dECdlognpo; 
}

double NucleusBase::Nucleusmup (const EOSData& eosIn, double ne, double uo, double ni) const {
	double npo = eosIn.Np();
	double T = eosIn.T();
	double v = GetVolume(eosIn, ne);
	double u = v*(ne - npo) / ((double)NucleusBase::mZ - v*npo);
	double EC = CoulombEnergy(v, npo, ne);
	double DuDnpo = v / ((double)NucleusBase::mZ - v*npo)*(u-1.0);
	double DDu = 0.5*(1 - pow(u, -2.0/3.0))/(0.5*pow(u,5.0/3.0) - 1.5*u + 1.e-40);
	double DFDnpo = EC*(DDu*DuDnpo - v/((double)NucleusBase::mZ - v*npo)/npo);
	return ni/uo * DFDnpo;
}

double NucleusBase::Nucleusmun (const EOSData& eosIn, double ne, double uo, double ni) const {
	return 0.0;
}

