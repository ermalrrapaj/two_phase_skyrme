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

  double T  = eosIn.T();
  double nb = eosIn.Nb(); 
  double nQ = pow(Constants::NeutronMassInFm * T / (2.0 * Constants::Pi), 1.5);
  double lam = 1.0;
  double lamV = 1.0;
  auto nse_funcs = [&](std::vector<double> xx)->std::vector<double> {
    std::vector<double> yy(2, 0.0);
    double uNuclei = 0.0;
    EOSData eosOut;
    double p;
    try {
    eosOut = mpEos->FromNAndT(EOSData::InputFromTNnNp(T, 
        exp(xx[0]), exp(xx[1])));
    p=eosOut.P();
    } catch (std::logic_error& e) { 
      std::cout << e.what() << " " << exp(xx[0]) << " " << exp(xx[1]) << std::endl; 
      throw e;
    }
    double mun = eosOut.Mun(); 
    double mup = eosOut.Mup();
    
    for (auto& nuc : mNuclei) {
      double v = nuc->GetVolume(eosOut, eosIn.Np());
      double BE = nuc->GetBindingEnergy(eosOut, eosIn.Np(), v);
      double aa = nuc->GetN()*mun + nuc->GetZ()*mup + (lam*BE - lamV*v*eosOut.P()); 
      double nn = nQ * pow((double) nuc->GetA(), 1.5) * exp(aa/T);
      //std::cout << nn << " " << v  << " " << eosOut.P() << " " << BE << " " << std::endl;
      yy[0] += nuc->GetN()*nn;
      yy[1] += nuc->GetZ()*nn;
      uNuclei += lamV*v*nn; 
      //std::cout << aa << " " << nn << " " << BE*197.3/nuc->GetA() << " ";
    }
    //std::cout << yy[0] << " " << yy[1] << " ";
    yy[0] += (1.0 - uNuclei)*eosOut.Nn(); 
    yy[1] += (1.0 - uNuclei)*eosOut.Np(); 
    yy[0] = yy[0]/((1.0-eosIn.Ye())*nb) - 1.0;
    yy[1] = yy[1]/(eosIn.Ye()*nb) - 1.0;
    ///std::cout << T*197.3 << " " << yy[0] << " " <<yy[1] << " " << 
    //uNuclei << " " << eosOut.Nb() << std::endl;
    return yy;
  };

  // Try some tricks to find the NSE solution
  MultiDimensionalRoot rootFinder = MultiDimensionalRoot(1.e-10, 2000);
  std::vector<double> pars = {log(eosIn.Nn()*1.e0), log(eosIn.Np()*1.e0)};
  double delta = 0.5; 
  int nfail = 0;
  lam = -1.e5;
  lamV = 0.0;
  do {
    lam += delta;
    if (lam>1.0) lam = 1.0;
    try {  
      pars = rootFinder(nse_funcs, pars, 2);
    } catch(...) { 
      lam -= delta;
      delta = 0.8*delta;
      ++nfail;
      if (delta<=1.e-6) 
          throw std::overflow_error("Too many tries in BE iter. " 
              + std::to_string(lam)
              + " " + std::to_string(delta));
      continue;
    }
    delta *= 1.01; 
  } while (lam< 1.0);
  
  lamV = 1.e-14;  
  delta = 100.0;
  do {
    lamV *= delta;
    if (lamV>1.0) lamV = 1.0;
    try {  
      pars = rootFinder(nse_funcs, pars, 2);
    } catch(...) { 
      lamV /= delta;
      delta = 0.8*(delta - 1.0) + 1.0;
      ++nfail;
      if (nfail>100)
          throw std::overflow_error("Too many tries in v iter. " 
              + std::to_string(lam)
              + " " + std::to_string(delta));
    }
  } while (lamV< 1.0);

  //for (lam = -1.e2; lam <= 1.0; lam+=0.01) {
  //  //std::cout << lam << std::endl;
  //  pars = rootFinder(nse_funcs, pars, 2);
  //}
  //lam = 1.0;
  //for (lamV = 1.e-10; lamV <= 1.0; lamV*=1.0001) {
  //  //std::cout << lamV << std::endl;
  //  pars = rootFinder(nse_funcs, pars, 2);
  //}
  //lamV = 1.0;
    
  //double Tsafe = pow(eosIn.Nb(), 1.0/3.0)*5.0;
  //double Nbsafe = pow(eosIn.T()/5.0, 3.0);
  //lamV = 1.0;
  //lam = 1.0;
  ////nb = Nbsafe;
  //if (T < Tsafe) {
  //  for (T = Tsafe; T >= eosIn.T(); T *= 0.999) {
  //    try {
  //      pars = rootFinder(nse_funcs, pars, 2);
  //    } catch (MultiDRootException& e) {
  //      if (e.GetError()<1.e-1) {
  //        pars = e.GetX();
  //      } else {
  //        std::cout << "temp " << nb << " " << T*197.3 <<std::endl;
  //        throw e;
  //      }  
  //    }
  //  }
  //}
  //}
  //T = eosIn.T();
  //if (nb >= Nbsafe) {
  //  for (nb = Nbsafe; nb <= eosIn.Nb(); nb *= 1.01) {
  //    try {
  //      pars = rootFinder(nse_funcs, pars, 2);
  //    } catch (MultiDRootException& e) {
  //      if (e.GetError()<1.e-3) {
  //        pars = e.GetX();
  //      } else {
  //        std::cout << "dens " << nb << " " << T*197.3 <<std::endl;
  //        throw e;
  //      } 
  //    } catch(std::logic_error& e) {
  //      std::cout << "Logic error : " << e.what()<< std::endl;
  //      throw e; 
  //    } catch(...) {
  //      std::cout << "Other "<< nb << " " << T*197.3 <<std::endl;
  //      throw; 
  //    }
  //  }
  //}

   
  nb = eosIn.Nb();
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
    double aa = nuc->GetN()*mun + nuc->GetZ()*mup + lam*BE - lamV*v*eosOut.P(); 
    dens.push_back(nQ * pow((double) nuc->GetA(), 1.5) * exp(aa/T));
  }  
  return dens;
}
