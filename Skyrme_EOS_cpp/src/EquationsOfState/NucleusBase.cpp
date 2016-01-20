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


namespace {
inline double D(double u) {
  return 1.0 - 1.5 * pow(u, 1.0/3.0) + 0.5*u;
}
inline double uDpoD(double u) {
  return (-0.5 * pow(u, 1.0/3.0) + 0.5*u)/(D(u) +  1.e-40);
}
}

double NucleusBase::CoulombEnergy(double v, double npo, double ne) const {
  double u = v*ne / (double)NucleusBase::mZ;
  double Ec = 3.0 * Constants::ElementaryChargeSquared 
      / (5.0 * pow(3.0 / (4.0 * Constants::Pi), 1.0/3.0))
      * pow((double)NucleusBase::mZ, 2) * (0.5*pow(u,2.0/3.0) - 1.5) 
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


StaticNucleus LDNucleus::GetStaticNucleus() const {
  double n0 = 1.e-10;
  EOSData eosIn = mpEos->FromNAndT(
      EOSData::InputFromTNbYe(1.e-3/197.3, n0, 0.5));
  double v = GetVolume(eosIn, n0); 
  double BE = GetBindingEnergy(eosIn, n0, v);
  int A = NucleusBase::mZ + NucleusBase::mN;
  std::cerr << NucleusBase::mZ << " " << NucleusBase::mN << " " << (double) A/v
      << " " << BE/A*197.3 << std::endl;
  return StaticNucleus(NucleusBase::mZ, A, BE, {}, {}, v);
}

double StaticNucleus::FreeEnergy(const EOSData& eosIn, double ne, double ni) const {
	double T = eosIn.T();
	double BE = GetBindingEnergy(eosIn, ne);
	double nQ = pow(MNUC*T/2/Constants::Pi,1.5); 
	double Fk = T*log((ni+1.e-100)/nQ/pow(NucleusBase::mA,1.5))-T;
	double FE= Fk-BE;
	return FE;	
}

double StaticNucleus::Entropy (const EOSData& eosIn, double ne, double ni) const {
	return 0.0;
}

double StaticNucleus::NucleusPressure (const EOSData& eosIn, double ne, double uo) const{
	double npo = 0.0;//eosIn.Np();
	double T = eosIn.T();
	double v = GetVolume(eosIn, ne);
	double u = v*(ne - npo) / ((double)NucleusBase::mZ - v*npo+  1.e-40);
	double EC = CoulombEnergy(v, npo, ne);
	double DuDlogne = v*ne / ((double)NucleusBase::mZ - v*npo+  1.e-40);
	double DuDlognpo = npo * ( v*u/((double)NucleusBase::mZ - v*npo+  1.e-40)
	                     -v/((double)NucleusBase::mZ - v*npo+  1.e-40) );
	double DDu = 0.5*(1 - pow(u, -2.0/3.0))/(D(u) +  1.e-40);
	
	double DFDlogne = EC*DDu*DuDlogne;
	double DFDlognpo = EC*(DDu*DuDlognpo - v/((double)NucleusBase::mZ - v*npo+  1.e-40));
	double pressure = T +uo *DFDlognpo + DFDlogne;
	return pressure;
}	

double StaticNucleus::Nucleusmup (const EOSData& eosIn, double ne, double uo, double ni) const {
	double npo = eosIn.Np();
	double T = eosIn.T();
	double v = GetVolume(eosIn, ne);
	double u = v*(ne - npo) / ((double)NucleusBase::mZ - v*npo);
	double EC = CoulombEnergy(v, npo, ne);
	double DuDnpo = v*u/((double)NucleusBase::mZ - v*npo)
	                     -v/((double)NucleusBase::mZ - v*npo) ;
	double DDu = 0.5*(1 - pow(u, -2.0/3.0))/D(u);
	double DFDnpo = EC*(DDu*DuDnpo - v/((double)NucleusBase::mZ - v*npo)/npo);
	double mup = ni/uo * DFDnpo;
	return mup;
}

double StaticNucleus::Nucleusmun (const EOSData& eosIn, double ne, double uo, double ni) const {
	return 0.0;
}

double LDNucleus::GetVolume(const EOSData& eosIn, double ne) const {
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
  
  OneDimensionalRoot rootFinder(1.e-8); 
  double vlo, vhi;
  if (eosIn.Np() <= ne) {
    // Just choose an extremely small volume as the lower bound
    vlo = 1.e-8 * (N + Z);
    // The maximum volume allowed for the given background electron density 
    // if the nuclear volume were larger, nuclei of this species would 
    // completely fill the space. 
    vhi = 0.9999 * Z / ne;
  } else {
   vlo = 1.000001* Z / ne; 
   vhi = 1.e12 * (N+Z); 
  }
  try {
    return rootFinder(pFunc, vlo, vhi);
  } catch(...) { 
    vlo = (N+Z)/0.2;
    vhi = (N+Z)/1.e-3;
    return rootFinder(pFunc, vlo, vhi);
  }
}

double LDNucleus::GetBindingEnergy(const EOSData& eosIn, double ne) const {
  return GetBindingEnergy(eosIn, ne, GetVolume(eosIn, ne));
}

double LDNucleus::GetBindingEnergy(const EOSData& eosIn, 
    double ne, double v) const {
  double Z = (double) NucleusBase::mZ; 
  double N = (double) NucleusBase::mN;
  double A = (double) NucleusBase::mA;
  double T = eosIn.T(); 
   
  EOSData eosBulk = mpEos->FromNAndT(EOSData::InputFromTNnNp(T, N/v, Z/v)); 
  return -(eosBulk.E() - eosBulk.T()*eosBulk.S()) * A 
         - SurfaceEnergy(v) - CoulombEnergy(v, eosIn.Np(), ne);

}

double LDNucleus::FreeEnergy(const EOSData& eosIn, double ne, double ni) const {
	double T = eosIn.T();
	double A = (double) NucleusBase::mA;
	double BE = GetBindingEnergy(eosIn,ne);
	double nQ = pow(MNUC*T/2/Constants::Pi,1.5); 
	double Fk = T*log(ni/nQ/pow(A,1.5))-T;
	double FE= Fk-BE;
	return FE;	
}

double LDNucleus::Entropy (const EOSData& eosIn, double ne, double ni) const {
	double v = GetVolume(eosIn, ne);
	double Z = (double) NucleusBase::mZ; 
	double N = (double) NucleusBase::mN;
	double A = (double) NucleusBase::mA;
	double T = eosIn.T();  
	double nQ = pow(MNUC*T/2/Constants::Pi,1.5);
	EOSData eosBulk = mpEos->FromNAndT(EOSData::InputFromTNnNp(T, N/v, Z/v));
	double Sk = 2.5-log(pow(A,-1.5)*ni/nQ);
	double Stot = Sk + A*eosBulk.S();
	return Stot;
}

double LDNucleus::SurfacePressure(double v) const {
  return -2.0 / 3.0 * pow(36.0 * Constants::Pi, 1.0/3.0) * mSigma0 
      * pow(v, -1.0/3.0); 
}

double LDNucleus::SurfaceEnergy(double v) const {
  return pow(36.0 * Constants::Pi, 1.0/3.0) * mSigma0 * pow(v, 2.0/3.0); 
}


double LDNucleus::CoulombEnergy(double v, double npo, double ne) const {
  double u = v*(ne - npo) / ((double)NucleusBase::mZ - v*npo);
  double Ec = 3.0 * Constants::ElementaryChargeSquared 
      / (5.0 * pow(3.0 * v / (4.0 * Constants::Pi), 1.0/3.0))
      * pow((double)NucleusBase::mZ - v*npo, 2) * D(u);
  return Ec;
}

double LDNucleus::CoulombPressure(double v, double npo, double ne) const {
  double Z = (double) NucleusBase::mZ;
  double denom = 1.0 / (Z/v - npo);
  double u = (ne - npo) * denom;
  return CoulombEnergy(v, npo, ne) / v  
      * (uDpoD(u)*denom*Z / (v*v) - 2.0*npo/(Z/v - npo) - 1.0/3.0 ); 
}

double LDNucleus::NucleusPressure (const EOSData& eosIn, double ne, double uo) const {
	double npo = eosIn.Np();
	double T = eosIn.T();
	double v = GetVolume(eosIn, ne);
	double u = v*(ne - npo) / ((double)NucleusBase::mZ - v*npo);
	double EC = CoulombEnergy(v, npo, ne);
	double DuDlogne = v*ne / ((double)NucleusBase::mZ - v*npo);
	double DuDlognpo = npo * ( v*u/((double)NucleusBase::mZ - v*npo)
	                     -v/((double)NucleusBase::mZ - v*npo) );
	double DDu = 0.5*(1 - pow(u, -2.0/3.0))/(D(u) +  1.e-40);
	
	double DFDlogne = EC*DDu*DuDlogne;
	double DFDlognpo = EC*(DDu*DuDlognpo - v/((double)NucleusBase::mZ - v*npo));
	double pressure = T +uo *DFDlognpo + DFDlogne;
	return pressure;
}

double LDNucleus::Nucleusmup (const EOSData& eosIn, double ne, double uo, double ni) const {
	double npo = eosIn.Np();
	double T = eosIn.T();
	double v = GetVolume(eosIn, ne);
	double u = v*(ne - npo) / ((double)NucleusBase::mZ - v*npo);
	double EC = CoulombEnergy(v, npo, ne);
	double DuDnpo = v*u/((double)NucleusBase::mZ - v*npo)
	                     -v/((double)NucleusBase::mZ - v*npo) ;
	double DDu = 0.5*(1 - pow(u, -2.0/3.0))/D(u);
	double DFDnpo = EC*(DDu*DuDnpo - v/((double)NucleusBase::mZ - v*npo)/npo);
	double mup = ni/uo * DFDnpo;
	return mup;
}

double LDNucleus::Nucleusmun (const EOSData& eosIn, double ne, double uo, double ni) const {
	return 0.0;
}
