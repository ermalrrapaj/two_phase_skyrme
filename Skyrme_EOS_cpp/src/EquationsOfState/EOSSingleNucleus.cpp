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

EOSData EOSSingleNucleus::FromNAndT(const EOSData& eosSave) {
  
  double llam = -5.0; // Starting point for turning on interactions

  // Use starting point above minimum temperature and work down in temperature later
  double Tsearch = 0.1/Constants::HBCFmMeV;
  EOSData eosIn = EOSData::InputFromTNnNp(eosSave.T(), eosSave.Nn(), eosSave.Np()); 
  if (eosIn.T() < Tsearch && eosSave.Nb() < 1.e-6) {  
    eosIn = EOSData::InputFromTNbYe(Tsearch, eosSave.Nb()*7.5e0, eosSave.Ye()); 
  } 
  else if ( eosIn.T() < Tsearch) { 
    eosIn = EOSData::InputFromTNbYe(Tsearch, eosSave.Nb(), eosSave.Ye()); 
  }

  // Find the Gibbs phases and try to find a single nucleus state using these 
  // as guesses 
  std::vector<EOSData> gibbsState = GibbsPhaseConstruct::FromNAndT(eosIn).Phases();
  std::vector<EOSData> snState;
  if (gibbsState.size()>1)
      snState = EquilibriumConditions(eosIn, gibbsState[0], gibbsState[1], llam);
  
  // The Gibbs phase boundary is not the same as the SN phase boundary
  // so try to work from the Gibbs phase boundary to lower density
  // if we haven't succesfully found a single nucleus state.
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
      snState = tmp;
      llam = 0.0;
      
    }
    if (snState.size()>1)
        snState = EquilibriumConditions(eosIn, snState[0], snState[1], llam);
  }  

  // The Gibbs phase boundary is not the same as the SN phase boundary
  // so try to work from the Gibbs phase boundary to higher density
  // if we haven't succesfully found a single nucleus state.
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
      snState = tmp;
      llam = 0.0;
    }
    if (snState.size()>1)
      snState = EquilibriumConditions(eosIn, snState[0], snState[1], llam);
  }  
  // If we have found the two single nucleus phases, calculate EoS 
  if (snState.size()>1) {
    
    // First, if we calculated at higher temperature work down to the desired
    // temperature
    if (eosIn.T() > eosSave.T()) {
      for (double Titer = eosIn.T(); Titer >= eosSave.T(); Titer*=0.99) {
        eosIn = EOSData::InputFromTNbYe(Titer, eosIn.Nb(), eosSave.Ye()); 
        snState = EquilibriumConditions(eosIn, snState[0], snState[1], 1.0);
      }
      // We also started from a higher density, so decrease the density
      for (double nbIter = eosIn.Nb(); nbIter >= eosSave.Nb(); nbIter*=0.99) {
        eosIn = EOSData::InputFromTNbYe(eosSave.T(), nbIter, eosSave.Ye()); 
        snState = EquilibriumConditions(eosIn, snState[0], snState[1], 1.0);
      }

      eosIn = EOSData::InputFromTNbYe(eosSave.T(), eosSave.Nb(), eosSave.Ye()); 
      snState = EquilibriumConditions(eosIn, snState[0], snState[1], 1.0);
    } 

    // Calculate the thermodynamic quantities for this mixed phase 
    double u = std::max((eosIn.Nn() - snState[0].Nn()) 
        / (snState[1].Nn() - snState[0].Nn()), 0.0);
    if (u < 1.e-30) return mpEos->FromNAndT(eosSave);
    EOSData eosSn = GetStateFromPhases(u, snState[0], snState[1]); 
    
    // Compare to bulk phase just for kicks 
    EOSData bulk = mpEos->FromNAndT(eosIn); 
    double F = eosSn.E() - eosSn.T()*eosSn.S();  
    double Fbulk = bulk.E()-bulk.T()*bulk.S();
    if ( Fbulk < F) {
      std::cout << " Bulk energy is lower than single nucleus prediction ";
      std::cout << eosIn.Nb() << " "<< F << " " << Fbulk << std::endl;
      //return bulk;
    }
    return eosSn; 
  } else {
    // We were unsuccesful at finding SN phases, so return bulk 
    EOSData out = mpEos->FromNAndT(eosSave);
    std::cout << "Returning bulk " << out.P() << std::endl;
    if (out.Phases().size()>1) std::cout << "Returning bulk Skyrme" << std::endl;
    return out;
  } 
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
  
  // Try the succesful parameters from last time function was called, since
  // it is often called on nearby points 
  std::vector<double> pars;
  if (mLastSet.size()==5) {
    try { 
      lamV = 1.0; 
      lamS = 1.0; 
      pars = rootFinder(equil_f, mLastSet, 5);
      mLastSet = pars; 
      return {mpEos->FromNAndT(EOSData::InputFromTNnNp(T, exp(pars[0]), exp(pars[1]))), 
              mpEos->FromNAndT(EOSData::InputFromTNnNp(T, exp(pars[0]) + exp(pars[2]), 
              exp(pars[1]) + exp(pars[3])))};
    } catch (...) {}
  }      
  
  // Try to slowly increase the effect of the finite size terms
  pars = {log(eosLo.Nn()), log(eosLo.Np()),
      log(eosHi.Nn() - eosLo.Nn()), log(eosHi.Np() - eosLo.Np()), 
      (eosIn.Np()-eosLo.Np())/(eosHi.Np() - eosLo.Np())};
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
   
  mLastSet = pars; 
  return {mpEos->FromNAndT(EOSData::InputFromTNnNp(T, exp(pars[0]), exp(pars[1]))), 
      mpEos->FromNAndT(EOSData::InputFromTNnNp(T, exp(pars[0]) + exp(pars[2]), 
      exp(pars[1]) + exp(pars[3])))};
  
}

EOSData EOSSingleNucleus::GetStateFromPhases(double u, const EOSData& low, 
    const EOSData& high) {

  double T = high.T();
  double muh = GetMuh(u, high.Nb(), T); 
  
  std::array<double, 2> D = DSurf(u); 
  std::array<double, 2> sig = Sigma(high.Ye()); 
  std::array<double, 3> h = HFunc(T, high.Ye());
  sig[0] = sig[0]*h[0];
  double beta = 4.5 * pow(8.0*Constants::Pi/15.0 
      *Constants::ElementaryChargeSquared , 1.0/3.0)
      * pow(high.Np()*sig[0], 2.0/3.0);
  
  // One way of calculating P  
  double P = low.P(); 
  if (u>0.0) {
    P += u*high.Nb()/mA0*h[0]*(T*(1.0-u) - u*muh);
    P -= beta*u*(D[0] - D[1]);
  }
  // Second way of calculating P
  //double P2 = snState[1].P(); 
  //P2 += beta*((2.0/3.0 - u)*D[0] - D[1]*(1.0-u));
  //P2 += u*(1.0-u)*snState[1].Nb()/mA0*h[0]*(muh + T);
  //P = P2;
   
  double S = u*high.Nb()*high.S() 
      + (1.0-u)*low.Nb()*low.S()
      + u*(1.0-u)*high.Nb()*h[0]/mA0*(5.0/2.0 - muh/T)
      - (2*beta*u*D[0]/(3*h[0]) + u*(1.0-u)*high.Nb()/mA0*(muh-T))*h[2]; 
  
  double F = (1.0-u)*low.Nb()*(low.E() - T*low.S());
  F += u*high.Nb()*(high.E() - T*high.S());
  F += beta*u*D[0];
  F += u*(1.0-u)*high.Nb()/mA0*h[0]*(muh-T);
  
  double E = F - T*S;
  double nn = u*low.Nn() + (1.0-u)*high.Nn(); 
  double np = u*low.Np() + (1.0-u)*high.Np();
  double nb = nn + np; 
  EOSData out = EOSData::Output(T, nn, np, low.Mun(),
      low.Mup(), P, S/nb, E/nb);
  out.SetPhases({low, high});
  return out;

}

double EOSSingleNucleus::GetMuh(double u, double nn, double T) {
  double nQ = pow((Constants::NeutronMassInFm*T/(2.0*Constants::Pi)), 1.5);
  return T*log(u*(1.0-u)*nn/(nQ*pow(mA0, 2.5))); 
}

std::array<double, 2> EOSSingleNucleus::DSurf(double u) {
  // Use small u expansion to reduce cancellation errors in derivative
  if (u < 1.e-6) {
    double uot = pow(u, 1.0/3.0);
    double D = 1.0 - 0.5*uot - 0.25*uot*uot - u/24.0;
    double Dp = 1.0 - 2.0/3.0*uot - 5.0/12.0*uot*uot - u/12.0;
    return {D, Dp};
  } 
    
  double a  = 1.0 - 1.5*pow(u, 1.0/3.0) + 0.5*u;
  double adu = -0.5*pow(u,-2.0/3.0) + 0.5; 
  double b  = 1.0 - 1.5*pow(1.0-u, 1.0/3.0) + 0.5*(1.0-u);
  double bdu = 0.5*pow(1.0-u,-2.0/3.0) - 0.5; 
  
  if (a<0.0 || b<0.0) return {0.0, 0.0};

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


