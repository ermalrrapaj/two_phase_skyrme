/// \file EOSNSE.cpp
/// \author lroberts
/// \since Sep 15, 2015
///
/// \brief
///
///

#include "EquationsOfState/EOSNSE.hpp"
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
    //std::cout << exp(xx[0]) << " " << exp(xx[1]) << " ";
    std::vector<double> yy(2, 0.0);
    double uNuclei = 0.0;
    EOSData eosOut = mpEos->FromNAndT(EOSData::InputFromTNnNp(T, 
        exp(xx[0]), exp(xx[1])));
    double mun = eosOut.Mun(); 
    double mup = eosOut.Mup();
     
    for (auto& nuc : mNuclei) {
      double v = nuc->GetVolume(eosOut, eosIn.Np());
      double BE = nuc->GetBindingEnergy(eosOut, eosIn.Np(), v);
      double aa = nuc->GetN()*mun + nuc->GetZ()*mup + BE - v*eosOut.P(); 
      double nn = nQ * pow((double) nuc->GetA(), 1.5) * exp(aa/T);
      yy[0] += nuc->GetN()*nn;
      yy[1] += nuc->GetZ()*nn;
      uNuclei += v*nn; 
      //std::cout << aa << " " << nn << " " << BE*197.3/nuc->GetA() << " ";
    }
    //std::cout << yy[0] << " " << yy[1] << " ";
    yy[0] += (1.0 - uNuclei)*eosOut.Nn(); 
    yy[1] += (1.0 - uNuclei)*eosOut.Np(); 
    yy[0] = yy[0]/eosIn.Nn() - 1.0;
    yy[1] = yy[1]/eosIn.Np() - 1.0;
    std::cout << T*197.3 << " " << yy[0] << " " <<yy[1] << " " << 
    uNuclei << " " << eosOut.Nb() << std::endl;
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
