/// \file EOSNSE_constructor.cpp
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

EOSNSE::EOSNSE(const std::vector<std::unique_ptr<NucleusBase>>& nuclei,
    const EOSBase& eos, bool buildGuessArray, double npomin, double npomax, 
    double nniMin, double nniMax, double T0) 
    : mTMin(0.1/Constants::HBCFmMeV), mpEos(eos.MakeUniquePtr()) { 
  for (auto& nuc : nuclei) 
    mNuclei.push_back(nuc->MakeUniquePtr());
  if (buildGuessArray) {
    for (double npo = npomin; npo<=npomax; npo *= 1.5) {
      nseGuesses.push_back(GetValidPoints(npo, T0, nniMin, nniMax));  
    }
  }
}

EOSNSE::EOSNSE(const EOSNSE& other) : 
    EOSNSE(std::vector<std::unique_ptr<NucleusBase>>(), *other.mpEos) {
  for (auto& nuc : other.mNuclei) {
    mNuclei.push_back(nuc->MakeUniquePtr());
  }
  nseGuesses = other.nseGuesses;
  mTMin = other.mTMin; 
}

std::vector<NSEProperties> EOSNSE::GetAllPoints(double npo, double T, 
    double nniMin, double nniMax) {
  std::vector<NSEProperties> allPts;
  double nn0 = nniMin;
  double delta = 0.01;
  double deltaMin = 1.e-2; 
  double deltaMax = 5.e-2;
  double nno = nn0;
  while (nn0<=nniMax) {
    try {
      auto nsev = GetExteriorProtonDensity(npo, nn0, T);
      allPts.insert(allPts.end(), nsev.begin(), nsev.end());
      if (nsev.size()>0) {
          NSEProperties nse = nsev.back();
          double npt = nse.npTot;
          double nnt = nse.nnTot;
          double npe = nse.eosExterior.Np();
          double nne = nse.eosExterior.Nn();

          if (fabs(nnt/nno - 1.0) > 0.1 && delta>deltaMin) {
            nn0 /= 1.0 + delta; 
            delta *= 0.5;
          } else { 
            if (fabs(nnt/nno - 1.0)<0.01 && delta<deltaMax) delta /= 0.5;
            nno = nnt;
          }
      }
    } catch(...) {}
    nn0 *= 1.0 + delta; 
  }
  return allPts;
}

bool EOSNSE::KeepPointBasedOnF(NSEProperties pt) {
  // Determines if the given point is on a branch we should keep
  // At low density, we need to keep all three branches while at high density 
  // we only want to keep one
  // The idea is to check if the free energy of the homogeneous phase is 
  // lower than the free energy of the input point and check that the properties 
  // of the homogeneous phase given NSE abundances that are effectively zero 
  
  auto bulk = mpEos->FromNAndT(
      EOSData::InputFromTNnNp(pt.T, pt.nnTot, pt.npTot)); 
  double fbulk = bulk.E() - pt.T*bulk.S();

  // Need to offset fbulk slightly to account for 
  // finite numerical precision 
  if (fbulk < 0.0) fbulk *= 1.0 + 1.e-4;
  else fbulk *= 1.0 - 1.e-4;

  // Get the properties of the nuclei for the given external conditions  
  auto nProp = GetTotalDensities(bulk, pt.npTot);
  ///todo Come up with some better criteria for choosing solutions 

  // This has a lower free energy than bulk matter would have, so keep it
  if (fbulk > pt.F) return true;
  
  // The bulk solution is invalid, since it would predict substantial NSE
  // abundances.  Therefore we need to keep this solution  
  if ((fabs((nProp.nnTot + nProp.npTot)/bulk.Nb()-1.0)>1.e-6 
          && bulk.Nb() < 0.04)) return true;

  // This is basically a heuristic that prevents us from cutting off some low 
  // density valid solutions
  if (bulk.Nn() < 1.e-5) return true;
  
  return false;
}

std::vector<NSEProperties> EOSNSE::GetValidPoints(double npo, double T, 
    double nniMin, double nniMax) {
  // Get all possible solutions 
  auto pts = GetAllPoints(npo, T, nniMin, nniMax); 
  
  // Keep the ones that satisfy our free energy criteria
  std::vector<NSEProperties> gdPts;
  for (auto& pt:pts) 
    if (KeepPointBasedOnF(pt)) gdPts.push_back(pt);
  
  // Sort the points by exterior neutron density
  std::sort(gdPts.begin(), gdPts.end(), 
      [](NSEProperties a, NSEProperties b) {
      return b.nnTot > a.nnTot;
  }); 

  // Tack on bulk solution at high density
  for (double nnb = gdPts.back().nnTot*1.01; nnb <= nniMax; nnb *=1.01) {
    gdPts.push_back(NSEProperties(mpEos->FromNAndT(
        EOSData::InputFromTNnNp(T, nnb, npo))));
  }

  // Tack on bulk solution at low density  
  for (double nnb = gdPts[0].nnTot*0.99; nnb >= nniMin; nnb *= 0.99) {
    gdPts.push_back(NSEProperties(mpEos->FromNAndT(
        EOSData::InputFromTNnNp(T, nnb, npo))));
  }
  
  return gdPts;
}


