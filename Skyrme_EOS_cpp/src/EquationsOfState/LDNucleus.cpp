/// \file NucleusBase.cpp
/// \author lroberts
/// \since Sep 14, 2015
///
/// \brief
///
///

#include "EquationsOfState/EOSData.hpp"
#include "EquationsOfState/NucleusBase.hpp"
#include "EquationsOfState/LDNucleus.hpp"
#include "Util/Constants.hpp"
#include "Util/OneDimensionalRoot.hpp"
#include <iostream> 
#include <math.h> 

static double HBC = Constants::HBCFmMeV; 
static double MNUC = 938.918/HBC;


StaticNucleus LDNucleus::GetStaticNucleus() const {
  double n0 = 1.e-10;
  EOSData eosIn = mpEos->FromNAndT(
      EOSData::InputFromTNbYe(1.e-3/197.3, n0, 0.5));
  double v = GetVolume(eosIn, n0); 
  double BE = GetBindingEnergy(eosIn, n0, v);
  int A = NucleusBase::mZ + NucleusBase::mN;
  return StaticNucleus(NucleusBase::mZ, A, BE, {}, {}, v);
}

LDNucleus::LDNucleus(int Z, int A, const EOSBase& eos) : NucleusBase(Z, A),
      mpEos(eos.MakeUniquePtr()),
      mSigma0(1.15/Constants::HBCFmMeV),
      mSs0(45.8/Constants::HBCFmMeV) {
    
  double n0 = 1.e-10;
  EOSData eosIn = mpEos->FromNAndT(EOSData::InputFromTNbYe(1.e-3/197.3, n0, 0.5));
  mV0 = -1.0;
  mV0 = LDNucleus::GetVolume(eosIn, n0); 
  
  double vmax = (double) A / 0.5; 
  double v = (double) A / 1.e-4;
  while (v>vmax) {
    EOSData eosBulk = mpEos->FromNAndT(
        EOSData::InputFromTNnNp(1.e-3/197.3, (double)(A-Z)/v, (double)Z/v)); 
    double Pb = eosBulk.P();
    double Ps = -SurfaceEnergy(v, 0.0, 0.0, 0.0)[1];
    double Pc = -CoulombEnergy(v, 0.0, 0.0, 0.0)[1];
    std::cout << A-Z << " " << v << " " << Pb << " " << Ps << " " << Pc << std::endl;
    PvsV.push_back(std::pair<double, double>(Pb+Ps+Pc, v)); 
    v /= 1.3;
  } 

} 

double LDNucleus::GetVolume(const EOSData& eosIn, double ne) const {
  double Z = (double) NucleusBase::mZ; 
  double N = (double) NucleusBase::mN;
  double T = eosIn.T(); 
  auto pFunc = [this, ne, T, Z, N, &eosIn](double v) -> double {
    EOSData eosBulk = mpEos->FromNAndT(
        EOSData::InputFromTNnNp(T, N/v, Z/v)); 
    double Pb = eosBulk.P();
    double Ps = -SurfaceEnergy(v, eosIn.Np(), eosIn.Nn(), ne)[1];
    double Pc = -CoulombEnergy(v, eosIn.Np(), eosIn.Nn(), ne)[1];
    return (Pb + Ps + Pc)/eosIn.P() - 1.0;  
  };
  
  OneDimensionalRoot rootFinder(1.e-8); 
  double vlo, vhi;
  //if (PvsV.size()>1) {
  //int idxl = 0;
  //int idxu = PvsV.size()-1;
  //double Pb = eosIn.P(); 
  //while (idxu-idxl>1) {
  //  int mid = (int)(idxl/2 + idxu/2);
  //  if (PvsV[mid].first > Pb) {
  //    idxl = mid; 
  //  } else { 
  //    idxu = mid; 
  //  } 
  //}

  //vlo = PvsV[idxl].second; 
  //vhi = PvsV[idxu].second;
  //
  //} else {
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
  //}
  try {
    return rootFinder(pFunc, vlo, vhi);
  } catch(...) { 
    return rootFinder(pFunc, (double) (N+Z)/0.2, (double) (N+Z)/1.e-3);
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
   
  EOSData eosBulk = GetBulk(T, v); 
  return -(eosBulk.E() - eosBulk.T()*eosBulk.S()) * A 
         - SurfaceEnergy(v, eosIn.Nn(), eosIn.Np(), ne)[0] 
         - CoulombEnergy(v, eosIn.Nn(), eosIn.Nn(), ne)[0];
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

std::vector<double> LDNucleus::SurfaceEnergy(double v, double nno, 
    double npo, double /*ne*/) const {
  double A = (double) NucleusBase::mN + (double) NucleusBase::mZ;
  double Es = pow(36.0 * Constants::Pi, 1.0/3.0) * mSigma0 * pow(v, 2.0/3.0);
  double dEsdv = 2.0/3.0*Es/v; 
  double dEsdnno = 0.0; 
  double dEsdnpo = 0.0;
  
  if (mV0>0.0) {
    Es *= pow(nno/A + npo/A - 1.0/v, 2) * mV0 * mV0;
    dEsdv   += 2.0*Es/(v*v*(nno/A + npo/A - 1.0/v)); 
    dEsdnno += 2.0*Es/(A*(nno/A + npo/A - 1.0/v)); 
    dEsdnpo += 2.0*Es/(A*(nno/A + npo/A - 1.0/v));
    
  } 
  // First value is surface energy 
  // second is derivative wrt v 
  // second is derivative wrt nno
  // second is derivative wrt npo
  // second is derivative wrt ne
  return {Es, dEsdv, dEsdnno, dEsdnpo, 0.0}; 
}

std::vector<double> LDNucleus::CoulombEnergy(double v, double /*nno*/, 
    double npo, double ne) const {
  double Z = (double) NucleusBase::mZ; 
  double u = v*(ne - npo) / (Z - v*npo);
  double D = 1.0 + 0.5*u - 1.5*pow(u,1.0/3.0); // Note the extra 1.0
  double udDduoD = (0.5*u - 0.5*pow(u,1.0/3.0))/(D + 1.e-40);
  double dludv   = 1.0 / (v + 1.e-40) + npo / (Z - v*npo);
  double dludne  = 1.0 / (ne - npo + 1.e-40);
  double dludnpo = (u-1.0) / (ne - npo + 1.e-40);
  double vn1o3 = pow(v, -1.0/3.0);
  double fac = 3.0 * Constants::ElementaryChargeSquared 
      / (5.0 * pow(3.0 / (4.0 * Constants::Pi), 1.0/3.0));

  double Ec      = fac * vn1o3 * D * pow(Z - npo*v, 2);
  double dEcdv   = Ec * (udDduoD*dludv - 1.0/(3.0 * v) - 2.0 * npo / (Z-npo*v));
  double dEcdnpo = Ec * (udDduoD*dludnpo - 2.0 * v / (Z-npo*v)); 
  double dEcdne  = Ec *  udDduoD*dludne; 
  
  // First value is coulomb energy 
  // second is derivative wrt v 
  // second is derivative wrt nno
  // second is derivative wrt npo
  // second is derivative wrt ne
  return {Ec, dEcdv, 0.0, dEcdnpo, dEcdne};
}

double LDNucleus::NucleusPressure (const EOSData& eosIn, double ne, double uo) const {
	double v = GetVolume(eosIn, ne);
	auto Ec = CoulombEnergy(v, eosIn.Nn(), eosIn.Np(), ne); 
	auto Es = SurfaceEnergy(v, eosIn.Nn(), eosIn.Np(), ne);
  uo += 1.e-20;
  return eosIn.T() + eosIn.Np()/uo * (Ec[3] + Es[3]) + eosIn.Nn()/uo * Es[2] 
      + ne * Ec[4];
}

double LDNucleus::Entropy (const EOSData& eosIn, double ne, double ni) const {
	double v = GetVolume(eosIn, ne);
	double Z = (double) NucleusBase::mZ; 
	double N = (double) NucleusBase::mN;
	double A = (double) NucleusBase::mA;
	double T = eosIn.T();  
	double nQ = pow(MNUC*T/2/Constants::Pi,1.5);
	EOSData eosBulk = GetBulk(T, v);
	double Sk = 2.5-log(pow(A,-1.5)*ni/nQ);
	double Stot = Sk + A*eosBulk.S();
	return Stot;
}

double LDNucleus::Nucleusmup (const EOSData& eosIn, double ne, double uo, 
    double ni) const {
	double v = GetVolume(eosIn, ne);
	auto Ec = CoulombEnergy(v, eosIn.Nn(), eosIn.Np(), ne); 
	auto Es = SurfaceEnergy(v, eosIn.Nn(), eosIn.Np(), ne);
  uo += 1.e-20;
	return ni/uo * (Ec[3] + Es[3]);
}

double LDNucleus::Nucleusmun (const EOSData& eosIn, double ne, double uo, double ni) const {
	double v = GetVolume(eosIn, ne);
	auto Ec = CoulombEnergy(v, eosIn.Nn(), eosIn.Np(), ne); 
	auto Es = SurfaceEnergy(v, eosIn.Nn(), eosIn.Np(), ne);
  uo += 1.e-20;
	return ni/uo * (Ec[2] + Es[2]);
}
