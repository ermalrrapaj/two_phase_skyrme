/// \file EOSNSE.cpp
/// \author lroberts
/// \since Sep 15, 2015
///
/// \brief
///
///

#include <algorithm> 
#include <omp.h>
#include "EquationsOfState/EOSNSE.hpp"
#include "Util/OneDimensionalRoot.hpp"
#include "Util/OneDimensionalMinimization.hpp"
#include "Util/MultiDimensionalRoot.hpp" 

EOSData EOSNSE::FromNAndT(const EOSData& eosIn) {
  return EOSData(); 
}

NSEProperties EOSNSE::GetExteriorDensities(const EOSData& eosIn, 
    const EOSData& extGuess) {
  double T = extGuess.T();
  double T0 = eosIn.T();
  auto nse_funcs = [&T,eosIn,this](std::vector<double> xx)->std::vector<double> {
    std::vector<double> yy(2, 0.0);
    EOSData eosOut = mpEos->FromNAndT(EOSData::InputFromTNnNp(T, 
        exp(xx[0]), exp(xx[1])));
    auto nucleiProps = GetNucleiScalars(eosOut, eosIn.Np()); 
     
    yy[0] = (1.0 - nucleiProps.uNuc)*eosOut.Nn() + nucleiProps.nn; 
    yy[1] = (1.0 - nucleiProps.uNuc)*eosOut.Np() + nucleiProps.np; 
    yy[0] = yy[0]/eosIn.Nn() - 1.0;
    yy[1] = yy[1]/eosIn.Np() - 1.0;
    return yy;
  };
  
  MultiDimensionalRoot rootFinder = MultiDimensionalRoot(1.e-6, 1000);
  std::vector<double> in = {log(extGuess.Nn()), log(extGuess.Np())};
  std::vector<double> res;
  
  // Search in temperature if the guess temperature is different from 
  // the desired temperature
  if (T0 < T) {
    double fac = 0.99;
    for (; T>T0; T *=fac)
      try { 
        in = rootFinder(nse_funcs, in, 2);
      } catch(MultiDRootException& e) {
        T /= fac;
        fac = 1.0 - (1.0 - fac)*0.8;
        if (1.0 - fac < 1.e-5) throw e; 
      }
  }

  if (T0 > T) {
    double fac = 1.01;
    for (; T<T0; T *=fac) {
      try {
        in = rootFinder(nse_funcs, in, 2);
      } catch(MultiDRootException& e) {
        T /= fac;
        fac = 1.0 + (fac - 1.0)*0.8;
        if (fac - 1.0 < 1.e-5) throw e; 
      }
    }
  }

  // Now find the solution at the actual temperature
  T = T0;
  res = rootFinder(nse_funcs, in, 2);

  // Find the nuclear properties and return the thermodynamic state
  EOSData eosOut = mpEos->FromNAndT(EOSData::InputFromTNnNp(T, 
        exp(res[0]), exp(res[1])));
  auto nucleiProps = GetNucleiScalars(eosOut, eosIn.Np());
   
  return GetStateNSEprop(NSEProperties(
          nucleiProps.nn + eosOut.Nn()*(1.0 - nucleiProps.uNuc), 
          nucleiProps.np + eosOut.Np()*(1.0 - nucleiProps.uNuc), 
          T, eosOut, nucleiProps.uNuc));
}

std::vector<double> EOSNSE::GetExteriorDensities(const EOSData& eosIn) {
  
  double T = eosIn.T();
  double nQ = pow(Constants::NeutronMassInFm * T / (2.0 * Constants::Pi), 1.5);
  auto nse_funcs = [&](std::vector<double> xx)->std::vector<double> {
    std::vector<double> yy(2, 0.0);
    EOSData eosOut = mpEos->FromNAndT(EOSData::InputFromTNnNp(T, 
        exp(xx[0]), exp(xx[1])));
    auto nucleiProps = GetNucleiScalars(eosOut, eosIn.Np()); 
     
    yy[0] = (1.0 - nucleiProps.uNuc)*eosOut.Nn() + nucleiProps.nn; 
    yy[1] = (1.0 - nucleiProps.uNuc)*eosOut.Np() + nucleiProps.np; 
    yy[0] = yy[0]/eosIn.Nn() - 1.0;
    yy[1] = yy[1]/eosIn.Np() - 1.0;
    return yy;
  };

  MultiDimensionalRoot rootFinder = MultiDimensionalRoot(1.e-9, 1000);
  std::vector<double> pars = {log(eosIn.Nn()*1.e0), log(eosIn.Np()*1.e0)};
  double Tsafe = pow(eosIn.Nb(), 1.0/3.0)*50.0;
  if (T < Tsafe) {
    for (T = Tsafe; T >= eosIn.T(); T *= 0.999) {
      try {
        pars = rootFinder(nse_funcs, pars, 2);
      } catch (MultiDRootException& e) {
        if (e.GetError()<1.e-1) {
          pars = e.GetX();
        } else {
          throw e;
        }  
      }
    }
  }

  T = eosIn.T(); 
  pars = rootFinder(nse_funcs, pars, 2);
  std::vector<double> dens = {exp(pars[0]), exp(pars[1])};
  EOSData eosOut = mpEos->FromNAndT(EOSData::InputFromTNnNp(T, 
      exp(pars[0]), exp(pars[1])));
  double mun = eosOut.Mun(); 
  double mup = eosOut.Mup();
  for (auto& nuc : mNuclei) {
    double v = nuc->GetVolume(eosOut, eosIn.Np());
    double BE = nuc->GetBindingEnergy(eosOut, eosIn.Np(), v);
    double aa = nuc->GetN()*mun + nuc->GetZ()*mup + BE - v*eosOut.P(); 
    dens.push_back(nQ * pow((double) nuc->GetA(), 1.5) * exp(aa/T));
  }  
  return dens;
}


std::vector<NSEProperties> EOSNSE::GetExteriorProtonDensity(double ne, 
    double nno, double T) {
  
  OneDimensionalRoot rootFinder1D(1.e-8, 100);
  
  double np_min = 1.e-300;
  double np_max = 1.2*ne;
  
  // Find a good interval 
  auto nse_func = [nno, ne, T, this](double xx)->double{
    double npo = exp(xx);
    EOSData eosOut = mpEos->FromNAndT(EOSData::InputFromTNnNp(T, nno, npo));
    auto nucleiProps = GetNucleiScalars(eosOut, ne); 
    double yy = 1.0 - (nucleiProps.np + (1.0 - nucleiProps.uNuc)*npo)/ne;
    return yy;
  };

  // Try finding the maximum of nse_func
  //double npmax = exp(minimize(nse_func, log(0.8*ne), log(1.e-300), log(ne)));  
  //double npmax2 = exp(minimize(nse_func, log(0.8*ne), log(0.7*ne), log(ne)));  
  //
  //if (nse_func(log(npmax))*nse_func(log(np_low)) < 0.0) np_high = npmax;
  //else if (nse_func(log(npmax))*nse_func(log(np_high)) < 0.0) np_low = npmax;
  
  double two, one, zero; 
  double ntwo, none, nzero; 
  
  double nstart = 1.e-5; 
  ntwo = nstart; none = nstart; nzero = nstart; 
  two = nse_func(log(nstart));
  one = nse_func(log(nstart));
  zero = nse_func(log(nstart));
  std::vector<double> extremum(1, 1.e-300);

  // The 1.1 seems to be critical, increasing it's value makes things not work 
  for (double npt = nstart; npt<=np_max; npt*=1.05) {
    ntwo = none;
    two = one;
    none = nzero;
    one = zero; 
    nzero = std::min(npt, ne);
    zero = nse_func(log(nzero));
    double npmax2 = -1.0;
    
    if (one>=two && one>=zero) {
      OneDimensionalMinimization minimize(1.e-7, 25);
      npmax2 = exp(minimize(nse_func, log(none), log(ntwo), log(nzero), true));  
      extremum.push_back(npmax2);
      //std::cerr << " Should be local maximum " << ntwo << " " << npmax2 << " " 
      //  << nzero; 
    }

    if (one<=two && one<=zero) {
      OneDimensionalMinimization minimize(1.e-7, 25);
      npmax2 = exp(minimize(nse_func, log(none), log(ntwo), log(nzero)));  
      extremum.push_back(npmax2);
      //std::cerr << " Should be local minimum " << ntwo << " " << npmax2 << " " 
      //  << nzero; 
    }

    if (npt >= ne) {
      OneDimensionalMinimization minimize(1.e-7, 25);
      npmax2 = exp(minimize(nse_func, log(none), log(ntwo), log(nzero), true));  
      if (npmax2>ntwo*(1.0+1.e-6) && npmax2<nzero*(1.0-1.e-6)) {
        extremum.push_back(npmax2);
      } else {
        npmax2 = -1.0;
      }
    }
  }
  
  //extremum.push_back(ne); 
  extremum.push_back(np_max);
  std::sort(extremum.begin(), extremum.end(), 
    [](double a, double b) {
      return b>a;
  });

  std::vector<double> npsol;
  std::vector<NSEProperties> output;
  for (int i=0; i<extremum.size()-1; i++) {
    double np_low = extremum[i];
    double np_high = extremum[i+1];
    if (nse_func(log(np_high))*nse_func(log(np_low))>0.0) continue;
    try { 
      double npo = exp(rootFinder1D(nse_func, log(np_low), log(np_high)));
      // Make sure we don't have two solutions that are the same
      //for (auto npt : npsol) {
      //  if (fabs(npo/npt - 1.0)<1.e-2) continue;
      //}
      npsol.push_back(npo);
      EOSData eosOut = mpEos->FromNAndT(EOSData::InputFromTNnNp(T, nno, npo));
      auto nucleiProps = GetNucleiScalars(eosOut, ne); 
      output.push_back(NSEProperties(
          nucleiProps.nn + eosOut.Nn()*(1.0 - nucleiProps.uNuc), 
          nucleiProps.np + eosOut.Np()*(1.0 - nucleiProps.uNuc), 
          T, eosOut, nucleiProps.uNuc));
    } catch (...) { std::cout << "Did not find root " << np_low << " " 
        << np_high << std::endl;}
  } 
  
  return output;
}

NSEProperties EOSNSE::GetTotalDensities(const EOSData& eosOut, double ne) {
  NucleiProperties nucleiProps = GetNucleiScalars(eosOut, ne);
  return NSEProperties(nucleiProps.nn + eosOut.Nn()*(1.0 - nucleiProps.uNuc), 
    nucleiProps.np + eosOut.Np()*(1.0 - nucleiProps.uNuc), eosOut.T(), eosOut, 
    nucleiProps.uNuc);
}

// find the total electron number density from the exterior
// proton and neutron number densities 
NSEProperties EOSNSE::GetTotalDensities(const EOSData& eosIn) {

  EOSData eosOut = mpEos->FromNAndT(eosIn);
  auto nse_func = [&eosOut, this](double xx)->double{
    auto nucleiProps = GetNucleiScalars(eosOut, exp(xx)); 
    double yy = (nucleiProps.np + (1.0 - nucleiProps.uNuc)*eosOut.Np())/exp(xx) 
        - 1.0;
    //std::cout << yy << " " << exp(xx) << std::endl;
    return yy;
  };

  OneDimensionalRoot rootFinder1D(1.e-7, 1000);
  double ne_low = log(1.e-5*eosIn.Np());
  double ne_high = log(std::max(eosIn.Np(), 0.1));
  double ne = exp(rootFinder1D(nse_func, ne_low, ne_high));
  auto nucleiProps = GetNucleiScalars(eosOut, ne); 
  return NSEProperties(nucleiProps.nn + eosOut.Nn()*(1.0 - nucleiProps.uNuc), 
    nucleiProps.np + eosOut.Np()*(1.0 - nucleiProps.uNuc), eosIn.T(), eosOut, 
    nucleiProps.uNuc);
}

EOSData EOSNSE::GetState(const NSEProperties& Prop){
  NSEProperties f = GetStateNSEprop(Prop);
  return EOSData::Output(f.T, f.nnTot, f.npTot, f.mun, f.mup, f.P, f.S, f.E);
}	

NSEProperties EOSNSE::GetStateNSEprop(const NSEProperties& Prop){
  
  auto nucTherm = GetNucleiScalars<true>(Prop.eosExterior, Prop.npTot); 
  
  double nn     = nucTherm.nn; 
  double np     = nucTherm.np; 
  double Ftot   = nucTherm.F;
  double Stot   = nucTherm.S;
  double uNuc   = nucTherm.uNuc;
  double vNuc   = nucTherm.vNuc;
  double Ptot   = nucTherm.P;
  double muntot = nucTherm.mun;
  double muptot = nucTherm.mup;

  double u0 = 1.0 - uNuc;
  double nbo = Prop.eosExterior.Nb();
  double T = Prop.eosExterior.T();
  nn     += u0*Prop.eosExterior.Nn();
  np     += u0*Prop.eosExterior.Np();
  // I think these should have a u0 out front
	Ftot   += u0*nbo*(Prop.eosExterior.E() - Prop.eosExterior.S()*T);
	Stot   += u0*nbo*Prop.eosExterior.S();
	muntot += Prop.eosExterior.Mun(); 
  muptot += Prop.eosExterior.Mup();
  Ptot   += Prop.eosExterior.P(); 
   
  Ftot/=(nn+np+1.e-40);
  Stot/=(nn+np+1.e-40);
  
  double Etot = Ftot + T*Stot;  
  double avgEc = nucTherm.avgEc;
  double avgBe = nucTherm.avgBe;
  double avgPv = nucTherm.avgP;
  return NSEProperties(nn, np, T, Prop.eosExterior, uNuc, avgEc, avgBe, avgPv, 
      Etot, Stot, Ftot, Ptot, muntot, muptot);
}

// Calculate the number densities and other properties of the nuclear 
// ensemble given a set of exterior conditions and total electron density
template <bool getEosContributions>
EOSNSE::NucleiProperties EOSNSE::GetNucleiScalars(const EOSData& eosOut, 
    double ne, double uog) {
  double mun = eosOut.Mun(); 
  double mup = eosOut.Mup(); 
  double T   = eosOut.T(); 
  double nQ = pow(Constants::NeutronMassInFm * T / (2.0 * Constants::Pi), 1.5);
  
  double nn = 0.0; 
  double np = 0.0;
  double uNuc  = 0.0;
  double vNuc  = 0.0;
  const double Pout = eosOut.P();
  #pragma omp parallel for reduction(+:nn,np,uNuc,vNuc)
  for (int i=0; i<mNuclei.size(); ++i) {
    double v = mNuclei[i]->GetVolume(eosOut, ne);
    double ni = mNuclei[i]->GetDensity(eosOut, ne, uog, v);
    
    nn += mNuclei[i]->GetN()*ni;
    np += mNuclei[i]->GetZ()*ni;
    uNuc += v*ni;
    vNuc += v; 
  }
  NucleiProperties out; 
  out.nn = nn;   
  out.np = np;   
  out.uNuc = uNuc;   
  out.vNuc = vNuc; 
    
  if (getEosContributions) { 
    double uo = 1.0 - uNuc; 
    double Ftot = 0.0, Stot = 0.0, Ptot = 0.0, muntot = 0.0, muptot = 0.0; 
    double avgEc = 0.0, avgBe = 0.0, avgP = 0.0, niTot = 0.0;
    // Currently have to have this second loop because nuclei contributions 
    // to thermodynamic quantities care about u0. May eventually fix how 
    // nucleus class returns EoS quantities
    #pragma omp parallel for default(shared) schedule(static) \
      reduction(+:Ftot, Stot, Ptot, muntot, muptot, avgEc, avgBe, niTot) 
    for (int i=0; i<mNuclei.size(); ++i) {
      double v = mNuclei[i]->GetVolume(eosOut, ne);
      double BE = mNuclei[i]->GetBindingEnergy(eosOut, ne, v);
      double ni = mNuclei[i]->GetDensity(eosOut, ne, 0.0, v);
      if (ni > 1.e-200) { // Ignore contributions from nuclei with small abundances
        Ftot   += ni*mNuclei[i]->FreeEnergy(eosOut, ne, ni);
        Stot   += ni*mNuclei[i]->Entropy(eosOut, ne, ni);
	      Ptot   += ni*mNuclei[i]->NucleusPressure(eosOut, ne, uo);
	      muntot += mNuclei[i]->Nucleusmun(eosOut, ne, uo, ni);
	      muptot += mNuclei[i]->Nucleusmup(eosOut, ne, uo, ni);
        avgEc  += mNuclei[i]->GetCoulombEnergy(eosOut, ne); 
        avgBe  += BE;
        niTot  += ni;
      }
    }
    avgEc /= niTot + 1.e-80;
    avgBe /= niTot + 1.e-80;
    avgP = Ptot/(niTot + 1.e-80);
    out.F = Ftot;
    out.S = Stot; 
    out.P = Ptot;
    out.mun = muntot;
    out.mup = muptot;
    out.avgEc = avgEc;
    out.avgBe = avgBe;
    out.avgP = avgP;
  } 
  return out;
}

