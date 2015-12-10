/// \file EOSNSE.cpp
/// \author lroberts
/// \since Sep 15, 2015
///
/// \brief
///
///

#include <algorithm> 
#include "EquationsOfState/EOSNSE.hpp"
#include "Util/OneDimensionalRoot.hpp"
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

NSEProperties EOSNSE::GetExteriorProtonDensity(double ne, double nno, 
    double T) {
  
  auto nse_func = [nno, ne, T, this](double xx)->double{
    double npo = exp(xx);
    EOSData eosOut = mpEos->FromNAndT(EOSData::InputFromTNnNp(T, nno, npo));
    auto nucleiProps = GetNucleiScalars(eosOut, ne); 
    double yy = (nucleiProps[1] + (1.0 - nucleiProps[2])*npo)/ne - 1.0;
    return yy;
  };

  // Find a good interval 
  double np_low = 1.e-300;
  double np_high = ne;
  while(fabs(nse_func(log(np_high)))>1.e10) np_high *= 0.1;
  if (nse_func(log(np_high))*nse_func(log(np_low))>0.0) {
    while (nse_func(log(np_high))*nse_func(log(np_low))>0.0) {
      np_low = np_high;
      np_high *= 10.0;
    }
  }
   
  // Find solution in interval
  OneDimensionalRoot rootFinder1D(1.e-7, 100);
  double npo = exp(rootFinder1D(nse_func, log(np_low), log(np_high)));
 
  // Set up output  
  EOSData eosOut = mpEos->FromNAndT(EOSData::InputFromTNnNp(T, nno, npo));
  auto nucleiProps = GetNucleiScalars(eosOut, ne); 
  return NSEProperties(nucleiProps[0] + eosOut.Nn()*(1.0 - nucleiProps[2]), 
    nucleiProps[1] + eosOut.Np()*(1.0 - nucleiProps[2]), T, eosOut, nucleiProps[2]);
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
  for (auto& nuc : mNuclei) {
    double v = nuc->GetVolume(eosOut, ne);
    double BE = nuc->GetBindingEnergy(eosOut, ne, v);
    double aa = nuc->GetN()*mun + nuc->GetZ()*mup + BE - v*eosOut.P(); 
    double ni = std::min(nQ * pow((double) nuc->GetA(), 1.5) * exp(aa/T), 1.e200);
    nn += nuc->GetN()*ni;
    np += nuc->GetZ()*ni;
    uNuc += v*ni;
    vNuc += v; 
  }
  return {nn, np, uNuc, vNuc};
}

