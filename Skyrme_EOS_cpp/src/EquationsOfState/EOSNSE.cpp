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
     
    yy[0] = (1.0 - nucleiProps[2])*eosOut.Nn() + nucleiProps[0]; 
    yy[1] = (1.0 - nucleiProps[2])*eosOut.Np() + nucleiProps[1]; 
    yy[0] = yy[0]/eosIn.Nn() - 1.0;
    yy[1] = yy[1]/eosIn.Np() - 1.0;
    //std::cout << T*197.3 << " " << yy[0] << " " <<yy[1] << " " << 
    //nucleiProps[2] << " " << eosOut.Nb() << std::endl;
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
    double yy = 1.0 - (nucleiProps[1] + (1.0 - nucleiProps[2])*npo)/ne;
    std::cerr << nno << " " << " " << nucleiProps[1] << " " << eosOut.Mun() 
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
  return NSEProperties(nucleiProps[0] + eosOut.Nn()*(1.0 - nucleiProps[2]), 
      nucleiProps[1] + eosOut.Np()*(1.0 - nucleiProps[2]), 
      T, eosOut, nucleiProps[2]);
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
    double yy = 1.0 - (nucleiProps[1] + (1.0 - nucleiProps[2])*npo)/ne;
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
      //std::cerr << "Expect bad interval." << std::endl;
    try { 
      double npo = exp(rootFinder1D(nse_func, log(np_low), log(np_high)));
      npsol.push_back(npo);
      EOSData eosOut = mpEos->FromNAndT(EOSData::InputFromTNnNp(T, nno, npo));
      auto nucleiProps = GetNucleiScalars(eosOut, ne); 
      output.push_back(NSEProperties(
      nucleiProps[0] + eosOut.Nn()*(1.0 - nucleiProps[2]), 
      nucleiProps[1] + eosOut.Np()*(1.0 - nucleiProps[2]), 
      T, eosOut, nucleiProps[2]));
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

EOSData EOSNSE::GetState(const NSEProperties& Prop){
	double T = Prop.eosExterior.T();
	double ne = Prop.npTot;
	double np0 = Prop.eosExterior.Np();
	double nn0 = Prop.eosExterior.Nn();
	double mun = Prop.eosExterior.Mun(); 
	double mup = Prop.eosExterior.Mup();  
	double muntot = mun, muptot	= muptot;
	double nQ = pow(Constants::NeutronMassInFm * T / (2.0 * Constants::Pi), 1.5);  
	double nn = 0.0; 
	double np = 0.0;
	double uNuc  = 0.0;
	double vNuc  = 0.0;
	double Ftot = (nn0+np0)*(Prop.eosExterior.E() - Prop.eosExterior.S()*T);
	double Stot = (nn0+np0)*Prop.eosExterior.S();
	double Ptot = Prop.eosExterior.P();
  
  // This is a good candidate for OpenMP parallelization
  #pragma omp parallel for 
  #pragma omp default(shared) 
  #pragma omp schedule(static) 
  #pragma omp reduce(+:nn,np,uNuc,vNuc,Ftot,Stot) 
  for (int i=0; i<mNuclei.size(); ++i) {
    double v = mNuclei[i]->GetVolume(Prop.eosExterior, ne);
    double BE = mNuclei[i]->GetBindingEnergy(Prop.eosExterior, ne, v);
    double aa = mNuclei[i]->GetN()*mun + mNuclei[i]->GetZ()*mup 
        + BE - v*Prop.eosExterior.P(); 
    double ni = std::min(nQ * pow((double) mNuclei[i]->GetA(), 1.5) 
        * exp(aa/T), 1.e200);
    nn += mNuclei[i]->GetN()*ni;
    np += mNuclei[i]->GetZ()*ni;
    Ftot += ni*mNuclei[i]->FreeEnergy(Prop.eosExterior,ne,ni);
    Stot += ni*mNuclei[i]->Entropy(Prop.eosExterior,ne,ni);
    uNuc += v*ni;
    vNuc += v;
  }
  
  double u0 = 1 - uNuc;
  nn += u0*nn0;
  np += u0*np0;
  
  // This is another good candidate for OpenMP parallelization
  #pragma omp parallel for 
  #pragma omp default(shared) 
  #pragma omp schedule(static) 
  #pragma omp reduce(+:Ptot,muntot,muptot)
  for (int i=0; i<mNuclei.size(); ++i){ 
	  double v = mNuclei[i]->GetVolume(Prop.eosExterior, ne);
	  double BE = mNuclei[i]->GetBindingEnergy(Prop.eosExterior, ne, v);
	  double aa = mNuclei[i]->GetN()*mun + mNuclei[i]->GetZ()*mup 
	  		+ BE - v*Prop.eosExterior.P(); 
	  double ni = std::min(nQ * pow((double) mNuclei[i]->GetA(), 1.5) 
	  		* exp(aa/T), 1.e200);
	  Ptot += ni*mNuclei[i]->NucleusPressure(Prop.eosExterior,ne,u0);
	  muntot += mNuclei[i]->Nucleusmun(Prop.eosExterior,ne,u0,ni);
	  muptot += mNuclei[i]->Nucleusmup(Prop.eosExterior,ne,u0,ni);
  }
  Ftot/=(nn+np);
  Stot/=(nn+np);
  double Etot = Ftot + T*Stot;
  EOSData totalData = EOSData::Output(T, nn, np, mun, mup, Ptot, Stot, Etot);
  return totalData;
}	


// find the total electron number density from the exterior
// proton and neutron number densities 
NSEProperties EOSNSE::GetTotalDensities(const EOSData& eosIn) {

  EOSData eosOut = mpEos->FromNAndT(eosIn);
  auto nse_func = [&eosOut, this](double xx)->double{
    auto nucleiProps = GetNucleiScalars(eosOut, exp(xx)); 
    double yy = (nucleiProps[1] + (1.0 - nucleiProps[2])*eosOut.Np())/exp(xx) - 1.0;
    //std::cout << yy << " " << exp(xx) << std::endl;
    return yy;
  };

  OneDimensionalRoot rootFinder1D(1.e-7, 1000);
  double ne_low = log(1.e-5*eosIn.Np());
  double ne_high = log(std::max(eosIn.Np(), 0.1));
  double ne = exp(rootFinder1D(nse_func, ne_low, ne_high));
  auto nucleiProps = GetNucleiScalars(eosOut, ne); 
  return NSEProperties(nucleiProps[0] + eosOut.Nn()*(1.0 - nucleiProps[2]), 
    nucleiProps[1] + eosOut.Np()*(1.0 - nucleiProps[2]), eosIn.T(), eosOut, 
    nucleiProps[2]);
}

std::array<double, 4> EOSNSE::GetNucleiScalars(const EOSData& eosOut, double ne) {
  double mun = eosOut.Mun(); 
  double mup = eosOut.Mup(); 
  double T   = eosOut.T(); 
  double nQ = pow(Constants::NeutronMassInFm * T / (2.0 * Constants::Pi), 1.5);
  
  double nn = 0.0; 
  double np = 0.0;
  double uNuc  = 0.0;
  double vNuc  = 0.0;
  // This is a good candidate for OpenMP parallelization
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
  return {nn, np, uNuc, vNuc};
}

