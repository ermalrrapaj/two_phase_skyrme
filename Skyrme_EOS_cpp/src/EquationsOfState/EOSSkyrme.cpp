/// \file EOSSkyrme.cpp
/// \author lroberts
/// \since Aug 18, 2015
///
/// \brief
///
///

#include <vector>
#include <math.h>

#include "EquationsOfState/EOSSkyrme.hpp" 
#include "Util/OneDimensionalRoot.hpp"
#include "Util/MultiDimensionalRoot.hpp"
#include "Util/Constants.hpp" 

static double HBC = Constants::HBCFmMeV; 
static double MNUC = 938.918/HBC;
static double PI = Constants::Pi; 
static double ALPHA = 3.0*HBC/(10.0*MNUC)*pow(3.0/2.0*PI*PI,2.0/3.0);
static double e_ele = sqrt(1.4299764/HBC);

extern "C" {
  double ifermi12_(double* scale_density, double * deriv);
  double zfermi32_(double* eta, double * deriv);
}

EOSSkyrme::EOSSkyrme() :
    mA(-286.1/HBC), 
    mB(-107.1/HBC), 
    mC(968.0/HBC),
    mD(0.0), 
    mF(0.0), 
    mG(0.0),
    mDelta(2.002) {}

std::vector<EOSData> EOSSkyrme::FromMuAndT(const EOSData& eosIn) const {
  
  MultiDimensionalRoot rootFinder(1.e-10);
  auto root_func = [&eosIn, this] (std::vector<double> logN) -> 
      std::vector<double> { 
    EOSData out = BaseEOSCall(eosIn.T(), exp(logN[0]), exp(logN[1]));
    return {(out.Mun() - eosIn.Mun())/(eosIn.Mun() + 1.e-7), 
            (out.Mup() - eosIn.Mup())/(eosIn.Mup() + 1.e-7)};
  };
  
  std::vector<EOSData> eosOut; 
   
  // Look for a low density solution first   
  try {
    std::vector<double> logN = rootFinder(root_func, {1.e-8, 1.e-8}, 2); 
    EOSData low = BaseEOSCall(eosIn.T(), exp(logN[0]), exp(logN[1]));
    eosOut.push_back(low);
  } catch(...) {}
  
  // Look for a high density solution second
  try {
    std::vector<double> logN = rootFinder(root_func, {1.e1, 1.e1}, 2); 
    EOSData hi = BaseEOSCall(eosIn.T(), exp(logN[0]), exp(logN[1]));
    eosOut.push_back(hi);
  } catch(...) {}
  return eosOut;
}

EOSData EOSSkyrme::FromNpMunAndT(const EOSData& eosIn) const {
  
  auto root_func = [&eosIn, this](double logNn)->double {  
      EOSData out = BaseEOSCall(eosIn.T(), exp(logNn), eosIn.Np());
      return (out.Mun() - eosIn.Mun()) / (eosIn.Mun() + 1.e-10);
  }; 
  
  OneDimensionalRoot rootFinder(1.e-12);
  double nn_lo = log(1.e-120);
  double nn_hi = log(0.95/(fabs(mF + mG) + 1.e-5));
  double logNn = rootFinder(root_func, nn_lo, nn_hi);
  
  return BaseEOSCall(eosIn.T(), exp(logNn), eosIn.Np());  
}

EOSData EOSSkyrme::FromNnMupAndT(const EOSData& eosIn) const {
  
  auto root_func = [&eosIn, this](double logNp)->double {  
      EOSData out = BaseEOSCall(eosIn.T(), eosIn.Nn(), exp(logNp)); 
      return (out.Mup() - eosIn.Mup()) / (eosIn.Mup() + 1.e-10);
  }; 
  
  OneDimensionalRoot rootFinder(1.e-12);
  double nn_lo = log(1.e-120);
  double nn_hi = log(0.95/(fabs(mF + mG) + 1.e-5));
  double logNp = rootFinder(root_func, nn_lo, nn_hi);
  
  return BaseEOSCall(eosIn.T(), eosIn.Nn(), exp(logNp));  
}

EOSData EOSSkyrme::FromNAndT(const EOSData& eosIn) {
  return BaseEOSCall(eosIn.T(), eosIn.Nn(), eosIn.Np()); 
} 

EOSSkyrme EOSSkyrme::FreeGas() {
  return EOSSkyrme(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
}

EOSSkyrme EOSSkyrme::FromErmalSkyrme(const std::array<const double, 7>& param) {
	double a= param[0], b = param[1], t0 = param[2], x0 = param[3];
	double t3=param[4], x3=param[5], alpha = param[6];
    double A, B, C, D, F, G, delta;
    F = (param[0]+param[1])*MNUC/4.0/HBC;
	G = -param[1]*MNUC/4.0/HBC;
	A = 0.25*param[2]*(1.0-param[3])/HBC; 
	B = 0.125*param[2]*(1.0+2.0*param[3])/HBC;
	C = param[4]*(1.0-param[5])/24.0/HBC;
	D = param[4]*(1.0+2.0*param[5])/48.0/HBC;
	delta = param[6]+1.0;
	// If needed to verify conversion!
	//std::cout<<A<<", "<<B<<", "<<C<<", "<<D<<", "<<F<<", "<<G<<", "<<delta<<"\n";
	EOSSkyrme out(A,B,C,D,F,G,delta);
	return out;
}

EOSSkyrme EOSSkyrme::FromSaturation(const std::array<const double, 7>& param){
	double rho_s = param[0], BE = param[1], MsoM = param[2], KM = param[3];
	double S = param[4], L = param[5], KS = param[6]; 
  double A, B, C, D, F, G, a, b, c, d, f, g, delta;
  f = (1.0- MsoM )/(MsoM *rho_s);
	delta = ( KM+2.0*pow(rho_s,2.0/3.0)*(1.0-5.0*f*rho_s)*ALPHA )
	        /( 3.0*pow(rho_s,2.0/3.0)*ALPHA*(1.0-2.0*f*rho_s)-9.0*BE );
	g = ( 9.0*KS-27.0*(L-3.0*S)*delta +5.0*pow(rho_s,2.0/3.0)*ALPHA
	    *(2.0-3.0*delta+2.0*F*rho_s*(3.0*delta-5.0)))/(30.0*pow(rho_s,5.0/3.0)*ALPHA*(3.0*delta-5.0)); 
	d = (5.0*(3.0*L-9.0*S+pow(rho_s,2.0/3.0)*ALPHA)-3.0*KS)/(9.0*(5.0-8.0*delta+3.0*pow(delta,2.0))*pow(rho_s,delta));
	c = (pow(rho_s,2.0/3.0)*(1.0-2.0*F*rho_s)*ALPHA-3.0*BE)/(3.0*(delta-1)*pow(rho_s,delta))-d;
	b = (L*(6.0+9.0*delta)+5.0*(ALPHA*pow(rho_s,2.0/3.0)*(3.0*delta-2.0)-9.0*S*delta)-3.0*KS)/(18.0*rho_s*(delta-1.0));
	a = -(2.0/3.0*ALPHA*pow(rho_s,-1.0/3.0)+b+5.0/3.0*F*ALPHA*pow(rho_s,2.0/3.0)
	    +(c+d)*delta*pow(rho_s,delta-1.0));
	A=a/HBC;
	B=b/HBC;
	C=c/HBC;
	D=d/HBC;
	F=f;
	G=g;
	// If needed to verify solution!
	//std::cout<<A<<", "<<B<<", "<<C<<", "<<D<<", "<<F<<", "<<G<<", "<<delta<<"\n";
	EOSSkyrme out(A,B,C,D,F,G,delta);
	return out;
} 

EOSData EOSSkyrme::BaseEOSCall(const double T, const double nn, 
    const double np) const {
  
  const double nt = nn + np; 
  const double xp = np/(nt+1.e-40);
  
  double momsp = 1.0 + mF*(nn + np) + mG*(nn - np);
  double momsn = 1.0 + mF*(nn + np) - mG*(nn - np);
  
  double invetan = 2.0*PI*PI * nn * pow(2.0*MNUC/momsn*T, -1.5); 
  double invetap = 2.0*PI*PI * np * pow(2.0*MNUC/momsp*T, -1.5); 
  double detan, detap;
  double etan = ifermi12_(&invetan, &detan); 
  double etap = ifermi12_(&invetap, &detap);
  double dtaun, dtaup; 
  double taup = pow(2.0*MNUC/momsp*T, 2.5)/(2.0*PI*PI) * zfermi32_(&etap, &dtaun);
  double taun = pow(2.0*MNUC/momsn*T, 2.5)/(2.0*PI*PI) * zfermi32_(&etan, &dtaup);
  
  double Up = (taup*(mF - mG) + taun*(mF + mG)) / (2.0*MNUC) 
      + 2.0*mA*nt + 4.0*mB*nn + mC*(1.0 + mDelta)*pow(nt, mDelta) 
      + 4.0*mD*pow(nt, mDelta-2.0)*(nn*nt + (mDelta-1.0)*nn*np);
  double Un = (taup*(mF + mG) + taun*(mF - mG)) / (2.0*MNUC) 
      + 2.0*mA*nt + 4.0*mB*np + mC*(1.0 + mDelta)*pow(nt, mDelta) 
      + 4.0*mD*pow(nt, mDelta-2.0)*(np*nt + (mDelta-1.0)*nn*np);
  
  double mup = etap*T + Up; 
  double mun = etan*T + Un; 
  
  double ee = taup*momsp/(2.0*MNUC) + taun*momsn/(2.0*MNUC) 
      + (mA + 4.0*mB*xp*(1.0-xp))*nt*nt 
      + (mC + 4.0*mD*xp*(1.0-xp))*pow(nt, mDelta + 1.0);
  double pp = (5.0/6.0*momsp - 0.5)/MNUC*taup + (5.0/6.0*momsn - 0.5)/MNUC*taun
      + (mA + 4.0*mB*xp*(1.0-xp))*nt*nt 
      + (mC + 4.0*mD*xp*(1.0-xp))*mDelta*pow(nt, mDelta + 1.0);
  double ss = 5.0/(6.0*MNUC*T)*(taup*momsp + taun*momsn) - np*etap - nn*etan; 
   
  return EOSData::Output(T, nn, np, mun, mup, pp, ss/nt, ee/nt);  
}



