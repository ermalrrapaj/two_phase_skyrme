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
  double T = eosIn.T();
  double nQ = pow(Constants::NeutronMassInFm * T / (2.0 * Constants::Pi), 1.5);
  auto nse_funcs = [&](std::vector<double> xx)->std::vector<double> {
    std::vector<double> yy(2, 0.0);
    double uNuclei = 0.0;
    EOSData eosOut = mpEos->FromNAndT(EOSData::InputFromTNnNp(T, 
        exp(xx[0]), exp(xx[1])));
    double mun = eosOut.Mun(); 
    double mup = eosOut.Mup();
     
    for (auto& nuc : mNuclei) {
      double v = nuc->GetVolume(eosOut, eosIn.Np());
      double aa = nuc->GetN()*mun + nuc->GetZ()*mup
           - nuc->GetBindingEnergy(eosOut, eosIn.Np(), v); 
      double nn = nQ * pow(nuc->GetA(), 1.5) * exp(aa/T);
      yy[0] += nuc->GetN()*nn;
      yy[1] += nuc->GetZ()*nn;
      uNuclei += v*nn; 
    }
    
    yy[0] += (1.0 - uNuclei)*eosOut.Nn(); 
    yy[1] += (1.0 - uNuclei)*eosOut.Np(); 
    yy[0] = yy[0]/eosIn.Nn() - 1.0;
    yy[1] = yy[1]/eosIn.Np() - 1.0;
    return yy;
  };

  MultiDimensionalRoot rootFinder = MultiDimensionalRoot(1.e-12, 100);
  std::vector<double> pars = rootFinder(nse_funcs, 
      {log(1.e-80), log(1.e-80)}, 2); 

  return EOSData(); 
}

