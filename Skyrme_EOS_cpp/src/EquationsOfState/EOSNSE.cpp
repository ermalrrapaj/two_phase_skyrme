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

EOSNSE::EOSNSE(const EOSNSE& other) : 
    EOSNSE(std::vector<std::unique_ptr<NucleusBase>>(), *other.mpEos) {
  for (auto& nuc : other.mNuclei) {
    mNuclei.push_back(nuc->MakeUniquePtr());
  } 
}

EOSData EOSNSE::FromNAndT(const EOSData& eosIn) {
  return EOSData(); 
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

NSEProperties EOSNSE::GetExteriorNeutronDensity(double ne, double npo, 
    double T) {
  
  OneDimensionalRoot rootFinder1D(1.e-8, 100);
  
  double nn_low = 1.e-10;
  double nn_high = 2.0;
  
  // Find a good interval 
  auto nse_func = [npo, ne, T, this](double xx)->double{
    double nno = exp(xx);
    EOSData eosOut = mpEos->FromNAndT(EOSData::InputFromTNnNp(T, nno, npo));
    auto nucleiProps = GetNucleiScalars(eosOut, ne); 
    double yy = 1.0 - (nucleiProps.np + (1.0 - nucleiProps.uNuc)*npo)/ne;
    std::cerr << nno << " " << " " << nucleiProps.np << " " << eosOut.Mun() 
        << " " << eosOut.Mup() << " " << yy << std::endl;
    return yy;
  };

  //while(fabs(nse_func(log(nn_high)))>1.e10) nn_high *= 0.9;
  
  //while (nse_func(log(nn_high))*nse_func(log(nn_low))>0.0) {
  //  nn_low = nn_high;
  //  nn_high *= 10.0;
  //}
  while (nse_func(log(nn_high))*nse_func(log(nn_low))>0.0) {
    nn_high *=0.99;
  }
  
  double nno = exp(rootFinder1D(nse_func, log(nn_low), log(nn_high)));
 
  // Set up output  
  EOSData eosOut = mpEos->FromNAndT(EOSData::InputFromTNnNp(T, nno, npo));
  auto nucleiProps = GetNucleiScalars(eosOut, ne); 
  return NSEProperties(nucleiProps.nn + eosOut.Nn()*(1.0 - nucleiProps.uNuc), 
      nucleiProps.np + eosOut.Np()*(1.0 - nucleiProps.uNuc), 
      T, eosOut, nucleiProps.uNuc);
}

std::vector<NSEProperties> EOSNSE::GetExteriorProtonDensity(double ne, 
    double nno, double T) {
  
  OneDimensionalRoot rootFinder1D(1.e-8, 100);
  OneDimensionalMinimization minimize(1.e-12, 1000);
  
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
  for (double npt = nstart; npt<=np_max; npt*=1.1) {
    ntwo = none;
    two = one;
    none = nzero;
    one = zero; 
    nzero = std::min(npt, ne);
    zero = nse_func(log(nzero));
    double npmax2 = -1.0;
    
    if (one>=two && one>=zero) {
      npmax2 = exp(minimize(nse_func, log(none), log(ntwo), log(nzero), true));  
      extremum.push_back(npmax2);
      //std::cerr << " Should be local maximum " << ntwo << " " << npmax2 << " " 
      //  << nzero; 
    }

    if (one<=two && one<=zero) {
      npmax2 = exp(minimize(nse_func, log(none), log(ntwo), log(nzero)));  
      extremum.push_back(npmax2);
      //std::cerr << " Should be local minimum " << ntwo << " " << npmax2 << " " 
      //  << nzero; 
    }

    if (npt >= ne) {
      npmax2 = exp(minimize(nse_func, log(none), log(ntwo), log(nzero), true));  
      if (npmax2>ntwo*(1.0+1.e-7) && npmax2<nzero*(1.0-1.e-7)) {
        extremum.push_back(npmax2);
        //std::cerr << " Maybe there is a maxima here " << ntwo << " " << npmax2 << " "; 
      } else {
        npmax2 = -1.0;
      }
    }
  }

  extremum.push_back(np_max);
  std::vector<double> npsol;
  std::vector<NSEProperties> output;
  for (int i=0; i<extremum.size()-1; i++) {
    double np_low = extremum[i];
    double np_high = extremum[i+1];
    if (nse_func(log(np_high))*nse_func(log(np_low))>0.0) continue;
    try { 
      double npo = exp(rootFinder1D(nse_func, log(np_low), log(np_high)));
      npsol.push_back(npo);
      EOSData eosOut = mpEos->FromNAndT(EOSData::InputFromTNnNp(T, nno, npo));
      auto nucleiProps = GetNucleiScalars(eosOut, ne); 
      output.push_back(NSEProperties(
          nucleiProps.nn + eosOut.Nn()*(1.0 - nucleiProps.uNuc), 
          nucleiProps.np + eosOut.Np()*(1.0 - nucleiProps.uNuc), 
          T, eosOut, nucleiProps.uNuc));
    } catch (...) {}
  } 
  
   
  //np_low = exp(minimize(nse_func, log(0.95*ne), log(0.9*ne), log(ne), true));  
  //np_high = ne; 
  //while(fabs(nse_func(log(np_high)))>1.e10) np_high *= 0.9;
  //
  //while (nse_func(log(np_high))*nse_func(log(np_low))>0.0) {
  //  np_low = np_high;
  //  np_high *= 10.0;
  //}
  //while (nse_func(log(np_high))*nse_func(log(np_low))>0.0) {
  //  np_high *=0.99;
  //}
  
  // Find solution in interval
  //if (nse_func(log(np_high))*nse_func(log(np_low))>0.0) std::cerr 
  //    << "Expect bad interval." << std::endl;
  //std::cerr << np_low << " " << np_high << " " << nse_func(log(0.9*ne)) << " " 
  //<< nse_func(log(np_low)) 
  //    << " " <<nse_func(log(np_high)) << std::endl;
  //double npo = exp(rootFinder1D(nse_func, log(np_low), log(np_high)));
  //std::cerr << np_low << " " << npo << " " << np_high << std::endl;
 
  // Set up output  
  //double npo = npsol.back();
  //EOSData eosOut = mpEos->FromNAndT(EOSData::InputFromTNnNp(T, nno, npo));
  //auto nucleiProps = GetNucleiScalars(eosOut, ne); 
  //return NSEProperties(nucleiProps[0] + eosOut.Nn()*(1.0 - nucleiProps[2]), 
  //    nucleiProps[1] + eosOut.Np()*(1.0 - nucleiProps[2]), 
  //    T, eosOut, nucleiProps[2]);
  return output;
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


void EOSNSE::SeTNSEdata (const std::vector<NSEProperties> & NSEDat) {
	NSEprop = NSEDat;
}
std::vector<NSEProperties> EOSNSE::GetNSEdata() {
  return NSEprop;
}

// Calculate the number densities and other properties of the nuclear 
// ensemble given a set of exterior conditions and total electron density
template <bool getEosContributions>
EOSNSE::NucleiProperties EOSNSE::GetNucleiScalars(const EOSData& eosOut, 
    double ne) {
  double mun = eosOut.Mun(); 
  double mup = eosOut.Mup(); 
  double T   = eosOut.T(); 
  double nQ = pow(Constants::NeutronMassInFm * T / (2.0 * Constants::Pi), 1.5);
  
  double nn = 0.0; 
  double np = 0.0;
  double uNuc  = 0.0;
  double vNuc  = 0.0;
  
  #pragma omp parallel for 
  #pragma omp default(shared) 
  #pragma omp schedule(static) 
  #pragma omp reduce(+:nn,np,uNuc,vNuc) 
  for (int i=0; i<mNuclei.size(); ++i) {
    double v = mNuclei[i]->GetVolume(eosOut, ne);
    double BE = mNuclei[i]->GetBindingEnergy(eosOut, ne, v);
    double aa = mNuclei[i]->GetN()*mun + mNuclei[i]->GetZ()*mup 
        + BE - v*eosOut.P(); 
    double ni = std::min(nQ * pow((double) mNuclei[i]->GetA(), 1.5) 
        * exp(aa/T), 1.e200);
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
    double u0 = 1.0 - uNuc; 
    double Ftot = 0.0, Stot = 0.0, Ptot = 0.0, muntot = 0.0, muptot = 0.0; 
    double avgEc = 0.0, avgBe = 0.0, avgP = 0.0, niTot = 0.0;
    // Currently have to have this second loop because nuclei contributions 
    // to thermodynamic quantities care about u0. May eventually fix how 
    // nucleus class returns EoS quantities
    #pragma omp parallel for 
    #pragma omp default(shared) 
    #pragma omp schedule(static) 
    #pragma omp reduce(+:Ftot, Stot, Ptot, muntot, muptot, avgEc, avgBe, niTot) 
    for (int i=0; i<mNuclei.size(); ++i) {
      double v = mNuclei[i]->GetVolume(eosOut, ne);
      double BE = mNuclei[i]->GetBindingEnergy(eosOut, ne, v);
      double aa = mNuclei[i]->GetN()*mun + mNuclei[i]->GetZ()*mup 
          + BE - v*eosOut.P(); 
      double ni = std::min(nQ * pow((double) mNuclei[i]->GetA(), 1.5) 
          * exp(aa/T), 1.e200);
      niTot  += ni;
      Ftot   += ni*mNuclei[i]->FreeEnergy(eosOut, ne, ni);
      Stot   += ni*mNuclei[i]->Entropy(eosOut, ne, ni);
	    Ptot   += ni*mNuclei[i]->NucleusPressure(eosOut, ne, u0);
	    muntot += mNuclei[i]->Nucleusmun(eosOut, ne, u0, ni);
	    muptot += mNuclei[i]->Nucleusmup(eosOut, ne, u0, ni);
      avgEc  += mNuclei[i]->CoulombEnergy(v, ne, eosOut.Np()); 
      avgBe  += BE;
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

