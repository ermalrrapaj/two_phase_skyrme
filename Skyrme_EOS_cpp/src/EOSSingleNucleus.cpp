/// \file EOSSingleNucleus.cpp
/// \authorr lroberts
/// \since Aug 30, 2015
///
/// \brief
///
///

#include <math.h>
#include <iostream> 
#include <algorithm> 
#include <iterator> 

#include "Constants.hpp"
#include "EOSSingleNucleus.hpp"
#include "MultiDimensionalRoot.hpp"
#include "OneDimensionalRoot.hpp" 

EOSData EOSSingleNucleus::FromNAndT(const EOSData& eosIn) {
  
  EOSData gibbsState = GibbsPhaseConstruct::FromNAndT(eosIn);
  
  if (gibbsState.Phases().size()>1) {
    
    auto snState = EquilibriumConditions(eosIn, gibbsState.Phases()[0],
        gibbsState.Phases()[1]);
    
    if (snState.size()>1) {
      double u = (eosIn.Nn() - snState[0].Nn()) 
          / (snState[1].Nn() - snState[0].Np());
      double h = 1.0; 
      double T = snState[1].T();
      double muh = GetNQ(u, snState[1].Nb(), T); 
      
      std::array<double, 2> D = DSurf(u); 
      std::array<double, 2> sig = Sigma(snState[1].Ye()); 
      double beta = 4.5 * pow(8.0*Constants::Pi/15.0, 1.0/3.0)
        * pow(Constants::ElementaryCharge*snState[1].Np()*sig[0], 2.0/3.0);
      
      double P = snState[0].P() + u*snState[1].Nb()/mA0*h*(T*(1-u) - u*muh);
      P -= beta*u*(D[0] - D[1]); 
      double S = u*snState[1].Nb()*snState[1].S() 
          + (1.0-u)*snState[0].Nb()*snState[0].S() 
          + u*(1.0-u)*snState[1].Nb()*h/mA0*(2.5 - muh/T);

      /// \todo Include h factor correction to Entropy 
       
      /// \todo Implement calculation of the internal energy 

      EOSData out = EOSData::Output(T, eosIn.Nn(), eosIn.Np(), snState[0].Mun(),
          snState[0].Mup(), P, S/eosIn.Nb(), 0.0);
      out.SetPhases(snState);
      return out;
    } else {
      return gibbsState;
    }
  } else {
    return gibbsState;
  } 
}

double EOSSingleNucleus::GetNQ(double u, double nn, double T) {
  double nQ = pow((Constants::NeutronMassInFm*T/(2.0*Constants::Pi)), 1.5);
  return T*log(u*(1.0-u)*nn/(nQ*pow(mA0, 2.5))); 
}

std::vector<EOSData> EOSSingleNucleus::EquilibriumConditions( 
    const EOSData& eosIn, const EOSData& eosLo, const EOSData& eosHi) {
  
  double T = eosIn.T();
  bool linear = true;
  double PScale = 1.0; 
  double MunScale = 1.0; 
  double MupScale = 1.0;
  double lamV = 0.0;
  double lamS = 0.0;

  // Define the equations to solve  
  auto equil_f = [this, T, &eosIn, &lamV, &lamS, &PScale, &MunScale, &MupScale]
      (std::vector<double> xx) -> std::vector<double> {
    double u = xx[4];
    double nnl = exp(xx[0]);
    double npl = exp(xx[1]);
    double nnu = exp(xx[0]) + exp(xx[2]);
    double npu = exp(xx[1]) + exp(xx[3]);

    EOSData eLo = mpEos->FromNAndT(EOSData::InputFromTNnNp(T, nnl, npl)); 
    EOSData eHi = mpEos->FromNAndT(EOSData::InputFromTNnNp(T, nnu, npu));
    
    // Calculate the Gibbs portion of the functions   
    std::vector<double> f = {eHi.P() - eLo.P(),
         eHi.Mun() - eLo.Mun(),
         eHi.Mup() - eLo.Mup(),
         ((1.0 - u)*eLo.Np() + u*eHi.Np())/eosIn.Np() - 1.0,
         ((1.0 - u)*eLo.Nn() + u*eHi.Nn())/eosIn.Nn() - 1.0};
    
    // Factor for approaching critical temperature  
    /// \todo Implement actual critical T factor
    double h = 1.0;
    double hp = 0.0; // Derivative of h wrt T
      
    // Calculate velocity corrections
    double muh = GetNQ(u, eHi.Nb(), T); 
    f[0] += lamV*h*muh/mA0*u*eHi.Nb();  
    f[1] += lamV*h*muh/mA0*(1.0-u); 
    f[2] += lamV*h*muh/mA0*(1.0-u); 

    // Calculate surface corrections
    std::array<double, 2> D = DSurf(u); 
    std::array<double, 2> sig = Sigma(eHi.Ye()); 
    double beta = 4.5 * pow(8.0*Constants::Pi/15.0, 1.0/3.0)
        * pow(Constants::ElementaryCharge*eHi.Np()*sig[0], 2.0/3.0);
    
    f[0] += lamS*beta*(2.0*D[0]/3.0 - D[1]); 
    f[1] -= lamS*2.0*beta*D[0]*eHi.Ye()*sig[1]/(3.0*eHi.Nb()*sig[0]); 
    f[2] -= lamS*2.0*beta*D[0]*(1.0-eHi.Ye())*sig[1]/(3.0*eHi.Nb()*sig[0]); 
    f[2] -= lamS*2.0*beta*D[0]/(3.0*eHi.Np()); 
       
    //Rescale the equations 
    f[0] = f[0]/PScale; 
    f[1] = f[1]/MunScale; 
    f[2] = f[2]/MupScale; 

    return f;
  };
  
  PScale = std::max(eosHi.P(), 1.e-8); 
  MunScale = fabs(eosHi.Mun()); 
  MupScale = fabs(eosHi.Mup()); 
  MultiDimensionalRoot rootFinder = MultiDimensionalRoot(1.e-12, 200);
  
  std::vector<double> pars = {log(eosLo.Nn()), log(eosLo.Np()),
      log(eosHi.Nn() - eosLo.Nn()), log(eosHi.Np() - eosLo.Np()), 
      (eosIn.Np()-eosLo.Np())/(eosHi.Np() - eosLo.Np())};
       
  // Try to slowly increase the effect of the finite size terms
  double llam = -6.0;
  double iter = 0;
  while (llam<=0.0 && iter < 100) {
    ++iter;
    lamS = std::min(pow(10.0, llam), 1.0); 
    lamV = 0.0; 
    try {
      pars = rootFinder(equil_f, pars, 5);
      llam = llam + log10(2.0);
    } catch (std::exception& e) {
      llam = llam - log10(0.8);
    }
    if (lamS >= 1.0) break;
  }

  if (iter <= 98) { 
    llam = -6.0;
    iter = 0;
    while (llam<=0.0 && iter < 100) {
      ++iter;
      lamS = 1.0; 
      lamV = std::min(pow(10.0, llam), 1.0);
      try {
        pars = rootFinder(equil_f, pars, 5);
        llam = llam + log10(2.0);
      } catch (std::exception& e) {
        llam = llam - log10(0.8);
      }
      if (lamV >= 1.0) break;
    }
  }

  if (iter>98) { 
      std::cerr << "Didn't converge. " << std::to_string(lamS) << std::endl;
      return {mpEos->FromNAndT(
          EOSData::InputFromTNnNp(T, eosIn.Nn(), eosIn.Np()))}; 
  }

  EOSData eL = mpEos->FromNAndT(
      EOSData::InputFromTNnNp(T, exp(pars[0]), exp(pars[1])));
  EOSData eH = mpEos->FromNAndT(
      EOSData::InputFromTNnNp(T, exp(pars[0]) + exp(pars[2]), 
      exp(pars[1]) + exp(pars[3])));
   
  return {eL, eH};
}

std::array<double, 2> EOSSingleNucleus::DSurf(double u) {
  
  double a  = 1.0 - 1.5*pow(u, 1.0/3.0) + 0.5*u + 1.e-50;
  double uda = -0.5*(u, 1.0/3.0) + 0.5*u;
  double b = std::max(1.0 - 1.5*pow(1.0-u, 1.0/3.0) + 0.5*(1.0-u) + 1.e-50,0.0);
  double mudb = 0.5*pow(1.0-u, 1.0/3.0) - 0.5*(1.0-u);
  double denom = u*u + (1.0-u)*(1.0-u) + 0.6*u*u*pow(1.0-u, 2.0); 
  double ddeno = 2.0*u - 2.0*(1.0 - u) + 1.2*u*pow(1.0-u, 2) - 1.2*u*u*(1.0-u);  
  
  // This is actually D/u 
  double D  = (1.0-u)*((1.0-u)*pow(a, 1.0/3.0) + u*pow(b, 1.0/3.0))/denom;
  
  double Dp = (1.0-2.0*u)*((1.0-u)*pow(a, 1.0/3.0) + u*pow(b, 1.0/3.0))/denom; 
  Dp +=(pow(1.0-u, 2)*pow(a, -2.0/3.0)/3.0*uda 
      + u*u*pow(b, -2.0/3.0)/3.0*mudb)/denom;
  Dp += u*(1.0-u)*(-pow(a, 1.0/3.0) + pow(b, 1.0/3.0))/denom;
  Dp -= u*(1.0-u)*((1.0-u)*pow(a, 1.0/3.0) 
      + u*pow(b, 1.0/3.0))/(denom*denom)*ddeno;
  
  return {D, Dp}; 

}

std::array<double, 2> EOSSingleNucleus::Sigma(double xp) {
  /// \todo Include h factor from L&S
  double r0 = pow(3.0/(4.0*Constants::Pi*0.155), 1.0/3.0); 
  double q = 384.0*Constants::Pi*r0*r0*mSigma0/mSs0 - 16.0; 
  double denom = q + pow(xp, -3.0) + pow(1.0 - xp, -3.0);
  double sig = mSigma0*(16.0 + q) / denom; 
  double sigp = sig / denom * (-3.0*pow(xp, -4.0) + 3.0*pow(xp, -4.0));
  return {sig, sigp};
}


