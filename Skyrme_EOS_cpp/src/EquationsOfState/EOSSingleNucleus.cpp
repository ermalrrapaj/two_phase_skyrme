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

#include "EquationsOfState/EOSSingleNucleus.hpp"
#include "Util/Constants.hpp"
#include "Util/MultiDimensionalRoot.hpp"
#include "Util/OneDimensionalRoot.hpp" 

EOSData EOSSingleNucleus::FromNAndT(const EOSData& eosIn) {
  double llam = -5.0; 
  std::vector<EOSData> gibbsState = GibbsPhaseConstruct::FromNAndT(eosIn).Phases();
  std::vector<EOSData> snState;
  if (gibbsState.size()>1)
      snState = EquilibriumConditions(eosIn, gibbsState[0], gibbsState[1], llam);
  
  // The Gibbs phase boundary is not the same as the SN phase boundary
  // so try to work from higher and lower density to the phase boundary
  if (snState.size()<2) { 
    for (double nb = eosIn.Nb()*1.0e1; nb > eosIn.Nb()*1.05; nb /= 1.05) {
      EOSData tmpIn =  EOSData::InputFromTNbYe(eosIn.T(), 
          std::max(nb, eosIn.Nb()), eosIn.Ye());
      
      std::vector<EOSData> tmp;
      if (snState.size()<2) {  
        tmp = GibbsPhaseConstruct::FromNAndT(tmpIn).Phases();
        if (tmp.size()<2) break;
      } else { 
        tmp = snState;
      }

      tmp = EquilibriumConditions(tmpIn, tmp[0], tmp[1]);
      if (tmp.size()<2) break;
      //std::cout << "Searching down from higher density. " << std::to_string(eosIn.Nb()) << " " 
      //    << std::to_string(nb)  << std::endl; 
      snState = tmp;
      llam = 0.0;
      
    }
    if (snState.size()>1)
        snState = EquilibriumConditions(eosIn, snState[0], snState[1], llam);
  }  

  if (snState.size()<2) {
    double nb; 
    for (nb = eosIn.Nb()*1.0e-1; nb < eosIn.Nb()*1.05; nb *= 1.05) {
      EOSData tmpIn =  EOSData::InputFromTNbYe(eosIn.T(), 
          std::min(nb, eosIn.Nb()), eosIn.Ye());
      
      std::vector<EOSData> tmp;
      if (snState.size()<2) {  
        tmp = GibbsPhaseConstruct::FromNAndT(tmpIn).Phases();
        if (tmp.size()<2) break;
      } else { 
        tmp = snState;
      }

      tmp = EquilibriumConditions(tmpIn, tmp[0], tmp[1]);
      if (tmp.size()<2) break;
      //std::cout << "Searching up from lower density. " << std::to_string(eosIn.Nb()) << " " 
      //    << std::to_string(nb)  << std::endl; 
      snState = tmp;
      llam = 0.0;
    }
    if (snState.size()>1)
      snState = EquilibriumConditions(eosIn, snState[0], snState[1], llam);
  }  

  if (snState.size()>1) {
    // Calculate the thermodynamic quantities for this mixed phase 
    double u = (eosIn.Nn() - snState[0].Nn()) 
        / (snState[1].Nn() - snState[0].Nn());
    double T = snState[1].T();
    double muh = GetMuh(u, snState[1].Nb(), T); 
    
    std::array<double, 2> D = DSurf(u); 
    std::array<double, 2> sig = Sigma(snState[1].Ye()); 
    std::array<double, 3> h = HFunc(T, snState[1].Ye());
    sig[0] = sig[0]*h[0];
    double beta = 4.5 * pow(8.0*Constants::Pi/15.0 
        *Constants::ElementaryChargeSquared , 1.0/3.0)
        * pow(snState[1].Np()*sig[0], 2.0/3.0);
    
    // One way of calculating P  
    double P = snState[0].P(); 
    P += u*snState[1].Nb()/mA0*h[0]*(T*(1.0-u) - u*muh);
    P -= beta*u*(D[0] - D[1]);

    // Second way of calculating P
    //double P2 = snState[1].P(); 
    //P2 += beta*((2.0/3.0 - u)*D[0] - D[1]*(1.0-u));
    //P2 += u*(1.0-u)*snState[1].Nb()/mA0*h[0]*(muh + T);
    //P = P2;
     
    double S = u*snState[1].Nb()*snState[1].S() 
        + (1.0-u)*snState[0].Nb()*snState[0].S()
        + u*(1.0-u)*snState[1].Nb()*h[0]/mA0*(5.0/2.0 - muh/T)
        - (2*beta*u*D[0]/(3*h[0]) + u*(1.0-u)*snState[1].Nb()/mA0*(muh-T))*h[2]; 
    
    double F = (1.0-u)*snState[0].Nb()*(snState[0].E() - T*snState[0].S());
    F += u*snState[1].Nb()*(snState[1].E() - T*snState[1].S());
    F += beta*u*D[0];
    F += u*(1.0-u)*snState[1].Nb()/mA0*h[0]*(muh-T);
    
    double E = F - T*S;
    
    EOSData bulk = mpEos->FromNAndT(eosIn); 
    double Fbulk = (bulk.E()-T*bulk.S())*eosIn.Nb();
    if ( Fbulk < F) {
      std::cout << " Bulk energy is lower than single nucleus prediction ";
      std::cout << eosIn.Nb() << " "<< F << " " << Fbulk << std::endl;
      //return bulk;
    }

    EOSData out = EOSData::Output(T, eosIn.Nn(), eosIn.Np(), snState[0].Mun(),
        snState[0].Mup(), P, S/eosIn.Nb(), E/eosIn.Nb());
    out.SetPhases(snState);
    return out;
  } else {
    //EOSData out = GibbsPhaseConstruct::FromNAndT(eosIn);
    //if (out.Phases().size()>1) std::cout << "Returning two phase gibbs" << std::endl;
    EOSData out = mpEos->FromNAndT(eosIn);
    if (out.Phases().size()>1) std::cout << "Returning bulk Skyrme" << std::endl;
    return out;
  } 
}

double EOSSingleNucleus::GetMuh(double u, double nn, double T) {
  double nQ = pow((Constants::NeutronMassInFm*T/(2.0*Constants::Pi)), 1.5);
  return T*log(u*(1.0-u)*nn/(nQ*pow(mA0, 2.5))); 
}

std::vector<EOSData> EOSSingleNucleus::EquilibriumConditions( 
    const EOSData& eosIn, const EOSData& eosLo, const EOSData& eosHi, 
    double llam0) {
  
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
    auto h = HFunc(T, eHi.Ye());
      
    // Calculate velocity corrections
    double muh = GetMuh(u, eHi.Nb(), T); 
    f[0] += lamV*h[0]*muh/mA0*u*eHi.Nb();  
    f[1] += lamV*(h[0]*muh - h[1]*eHi.Ye()*(muh-T))/mA0*(1.0-u); 
    f[2] += lamV*(h[0]*muh + h[1]*(1.0-eHi.Ye())*(muh-T))/mA0*(1.0-u); 

    // Calculate surface corrections
    std::array<double, 2> D = DSurf(u); 
    std::array<double, 2> sig = Sigma(eHi.Ye());
    sig[1] = h[1]/h[0] + sig[1]; // sigma'/sigma
    sig[0] = sig[0]*h[0];
    double beta = 4.5 * pow(8.0*Constants::Pi/15.0
        * Constants::ElementaryChargeSquared, 1.0/3.0)
        * pow(eHi.Np()*sig[0], 2.0/3.0);
    
    f[0] += lamS*2.0*beta*(D[0]/3.0 - D[1]/2.0); 
    f[1] -= lamS*2.0*beta*D[0]*sig[1]/(3.0*eHi.Nb())*eHi.Ye(); 
    f[2] += lamS*2.0*beta*D[0]*sig[1]/(3.0*eHi.Nb())*(1.0-eHi.Ye()); 
    f[2] += lamS*2.0*beta*D[0]/(3.0*eHi.Np()); 
       
    // Rescale the equations 
    f[0] = f[0]/PScale; 
    f[1] = f[1]/MunScale; 
    f[2] = f[2]/MupScale; 

    return f;
  };
  
  PScale = std::max(fabs(eosHi.P()), 1.e-6); 
  MunScale = fabs(eosHi.Mun()); 
  MupScale = fabs(eosHi.Mup()); 
  MultiDimensionalRoot rootFinder = MultiDimensionalRoot(1.e-10, 250);
  
  std::vector<double> pars;
  bool lastSuccess = false;  
  if (mLastSet.size()==5) {
    try { 
      lamV = 1.0; 
      lamS = 1.0; 
      pars = rootFinder(equil_f, mLastSet, 5);
      lastSuccess = true;
    } catch (...) {
      lastSuccess = false;
    }
  }      
  
  if (!lastSuccess) { 
    if (eosIn.T() < GibbsPhaseConstruct::mTMin) T = GibbsPhaseConstruct::mTMin;
    pars = {log(eosLo.Nn()), log(eosLo.Np()),
        log(eosHi.Nn() - eosLo.Nn()), log(eosHi.Np() - eosLo.Np()), 
        (eosIn.Np()-eosLo.Np())/(eosHi.Np() - eosLo.Np())};
    
    // Try to slowly increase the effect of the finite size terms
    double dllam = 0.1;
    for (double llam = llam0; llam <0.01; llam += dllam) {
      lamV = std::min(pow(10.0, llam), 1.0); 
      lamS = std::min(pow(10.0, llam), 1.0); 
      try {
        pars = rootFinder(equil_f, pars, 5);
      } catch (std::exception& e) {
        //std::cerr << "Didn't converge. " << std::to_string(lamV) << " " 
        //    << std::to_string(lamV-1.0) << " " << std::to_string(eosIn.Nb())  
        //    << " " << pars[4] << std::endl;
        //llam -= dllam; 
        //dllam = 0.7*dllam;
        //if (dllam < 0.01) 
          return {mpEos->FromNAndT(
              EOSData::InputFromTNnNp(T, eosIn.Nn(), eosIn.Np()))}; 
      }
      if (T>eosIn.T()) { 
        for (T = GibbsPhaseConstruct::mTMin; T >= eosIn.T(); T *= 0.9) {
          try {
            pars = rootFinder(equil_f, pars, 5);
          } catch (std::exception& e) {
            std::cerr << "Unable to work down to desired temperature" << std::endl;
            return {mpEos->FromNAndT(
              EOSData::InputFromTNnNp(eosIn.T(), eosIn.Nn(), eosIn.Np()))}; 
          }
        }
        T = eosIn.T();
        try {
          pars = rootFinder(equil_f, pars, 5);
        } catch(...) {
          std::cerr << "Unable to work down to desired temperature" << std::endl;
          return {mpEos->FromNAndT(
            EOSData::InputFromTNnNp(eosIn.T(), eosIn.Nn(), eosIn.Np()))}; 
        }
      }
    } 
  }
   
  EOSData eL = mpEos->FromNAndT(
      EOSData::InputFromTNnNp(T, exp(pars[0]), exp(pars[1])));
  EOSData eH = mpEos->FromNAndT(
      EOSData::InputFromTNnNp(T, exp(pars[0]) + exp(pars[2]), 
      exp(pars[1]) + exp(pars[3])));
  
  mLastSet = pars; 
  
  return {eL, eH};
}

std::array<double, 2> EOSSingleNucleus::DSurf(double u) {
  
  double a  = 1.0 - 1.5*pow(u, 1.0/3.0) + 0.5*u;
  double adu = -0.5*pow(u,-2.0/3.0) + 0.5; 
  double b  = 1.0 - 1.5*pow(1.0-u, 1.0/3.0) + 0.5*(1.0-u);
  double bdu = 0.5*pow(1.0-u,-2.0/3.0) - 0.5; 

  double denom = u*u + (1.0-u)*(1.0-u) + 0.6*u*u*pow(1.0-u, 2); 
  double ddeno = 2.0*u - 2.0*(1.0 - u) + 1.2*u*pow(1.0-u, 2) - 1.2*u*u*(1.0-u);  
  
  // This is actually D/u 
  double D  = (1.0-u)*((1.0-u)*pow(a, 1.0/3.0) + u*pow(b, 1.0/3.0))/denom;
  
  // This is D' 
  double Dp = (1.0-2.0*u)*((1.0-u)*pow(a, 1.0/3.0) + u*pow(b, 1.0/3.0))/denom; 
  Dp -= u*D/denom*ddeno; 
   
  Dp += u*(1.0-u)/denom * (1.0-u)*pow(a, -2.0/3.0)/3.0 * adu; 
  Dp -= u*(1.0-u)/denom * pow(a, 1.0/3.0);
 
  Dp += u*(1.0-u)/denom * u*pow(b, -2.0/3.0)/3.0 * bdu;
  Dp += u*(1.0-u)/denom * pow(b, 1.0/3.0);
  
  return {D, Dp}; 

}

std::array<double, 3> EOSSingleNucleus::HFunc(double T, double xp) {
  double Tc = 87.76/Constants::HBCFmMeV
      *sqrt(mKs0*Constants::HBCFmMeV/375.0)*xp*(1.0-xp);
  double dTc = 87.76/Constants::HBCFmMeV
      *sqrt(mKs0*Constants::HBCFmMeV/375.0)*(1.0-2.0*xp);
  if (T>Tc) return {0.0, 0.0};
  double t = 1.0 - T*T/(Tc*Tc);
  return {t*t, 4.0*t*T*T/pow(Tc, 3)*dTc, 4.0*t*T/(Tc*Tc)};
  //return {1.0, 0.0, 0.0};
}

std::array<double, 2> EOSSingleNucleus::Sigma(double xp) {
  double r0 = pow(3.0/(4.0*Constants::Pi*0.155), 1.0/3.0); 
  double q = 384.0*Constants::Pi*r0*r0*mSigma0/mSs0 - 16.0; 
  double denom = q + pow(xp, -3.0) + pow(1.0 - xp, -3.0);
  double sig = mSigma0*(16.0 + q) / denom;
  // sigp is sig'/sig to prevent divide by zeros
  double sigp = -(-3.0*pow(xp, -4.0) + 3.0*pow(1.0-xp, -4.0)) / denom;
  return {sig, sigp};
}


