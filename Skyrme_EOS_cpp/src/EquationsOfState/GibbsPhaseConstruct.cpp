/// \file GibbsPhaseConstruct.cpp
/// \authorr lroberts
/// \since Aug 18, 2015
///
/// \brief
///
///

#include <math.h>
#include <iostream> 
#include <algorithm> 
#include <iterator> 

#include "EquationsOfState/GibbsPhaseConstruct.hpp" 
#include "Util/Constants.hpp"
#include "Util/MultiDimensionalRoot.hpp"
#include "Util/OneDimensionalRoot.hpp"

GibbsPhaseConstruct::GibbsPhaseConstruct(const EOSBase& eos, 
    bool findPhaseBound, double TMeVMin, double TMeVMax) : 
    mTMult(1.25),
    mpEos(eos.MakeUniquePtr()), 
    mVerbose(false),
    mTMin(TMeVMin/Constants::HBCFmMeV),
    mTMax(TMeVMax/Constants::HBCFmMeV),
    mPhaseBounds(std::vector<std::vector<std::pair<EOSData, EOSData>>>()) {
  if (findPhaseBound) FindPhaseBoundary();
}

void GibbsPhaseConstruct::FindPhaseBoundary() {
  
  // Do a first pass for finding phase boundaries  
  double TFail = -1.0;
  double TMax = 50.0/Constants::HBCFmMeV;
  if (mTMax < 0.1/Constants::HBCFmMeV) 
    throw std::range_error("Minimum temperature must be > 0.1 MeV"); 
  for (double lT = log10(std::max(mTMin,0.1/Constants::HBCFmMeV)); 
      lT <= log10(TMax) && lT <= log10(mTMax); lT +=  log10(mTMult)) {
    std::cerr << pow(10.0, lT)*Constants::HBCFmMeV << " ";
    auto phaseBound = FindFixedTPhaseBoundary(pow(10.0, lT));
    std::cerr << phaseBound.size() << std::endl;
    if (phaseBound.size()>10) {
      mPhaseBounds.push_back(phaseBound);
      mTCrit = 1.01 * pow(10.0, lT);
    } else {
      TFail = pow(10.0, lT);
      break;
    }
  }

  // Find low temperature phase boundaries
  if (mTMin < 0.1) {
    for (double lT = log10(0.1/Constants::HBCFmMeV/mTMult); lT >= log10(mTMin); 
        lT -= log10(mTMult)){
      std::cerr << pow(10.0, lT)*Constants::HBCFmMeV << " ";
      auto phaseBound = FindFixedTPhaseBoundary(pow(10.0,lT), mPhaseBounds[0]); 
      std::cerr << phaseBound.size() << std::endl;
      mPhaseBounds.push_back(phaseBound);
    } 
  }
  
  // Now bisect to get close to the critical temperature
  if (TFail > 0.0) {
    double Tlo = mTCrit;
    double Thi = TFail;
    double Tmid = 0.0;
    for (int i=0; i<8; i++) {
      Tmid = 0.5*(Tlo + Thi); 
      std::cerr << Tmid*Constants::HBCFmMeV << " ";
      auto phaseBound = FindFixedTPhaseBoundary(Tmid);
      std::cerr << phaseBound.size() << std::endl;
      if (phaseBound.size()>10) { 
        mPhaseBounds.push_back(phaseBound);
        Tlo = Tmid; 
      } else { 
        Thi = Tmid;
      }
    }
    mTCrit = Tmid;
    std::cerr << "Tcritical : " << mTCrit*Constants::HBCFmMeV 
        << std::endl;
  } else {
    mTCrit = mTMax;
  }
   
  // Sort the phase boundaries in terms of temperature 
  std::sort (mPhaseBounds.begin(), mPhaseBounds.end(), 
      [](std::vector<std::pair<EOSData, EOSData>> a, 
      std::vector<std::pair<EOSData, EOSData>> b) { 
      return ((a[0].first).T() < (b[0].first).T());});
}

EOSData GibbsPhaseConstruct::FromNAndT(const EOSData& eosIn) {
  double T = eosIn.T(); 
  double nb = eosIn.Nb();
  double np = eosIn.Np();

  // Check that we are within the temperature bounds. 
  if (T < mTMin) {
    std::cerr << T << " " << mTMin << std::endl;
    throw std::logic_error("Trying to call the EoS below the minimum allowed T.");
  }
     
  if (T > mTCrit) {
    if (mVerbose) std::cout << "# We are above the critical temperature \n";
    return mpEos->FromNAndT(eosIn); 
  } 
    
  // Check if a close by temperature phase boundary has been calculated 
  std::vector<std::pair<EOSData, EOSData>> *phaseBound = NULL;
  phaseBound = &mPhaseBounds.back();
  for (int i=1; i<mPhaseBounds.size(); ++i) {
    if (T <= (mPhaseBounds[i])[0].first.T()) {
      phaseBound = &mPhaseBounds[i-1];
      break;
    }
  }
  
  // Didn't find a phase boundary that was close enough, bite the bullet and 
  // calculate a new one 
  if (!phaseBound) {
    if (mVerbose) std::cerr << "Calculating new phase boundary." << std::endl;
    mPhaseBounds.push_back(FindFixedTPhaseBoundary(T)); 
    phaseBound = &mPhaseBounds.back(); 
  } 
  
  if (phaseBound->size()<1) {
    if (mVerbose) std::cout << "# There are not two phases for this temperature \n";
    return mpEos->FromNAndT(eosIn); 
  }
   
  // Check to see if we are in or out of the mixed phase region
  std::vector<double> bracketNn;
  for (int i = 1; i < phaseBound->size(); ++i) {
    if ( ((*phaseBound)[i].first.Np() <= np) && ((*phaseBound)[i-1].first.Np() > np) 
        || ((*phaseBound)[i-1].first.Np() <= np) && ((*phaseBound)[i].first.Np() > np)) {
      bracketNn.push_back((*phaseBound)[i].first.Nn());
      bracketNn.push_back((*phaseBound)[i-1].first.Nn());
    }  
    if ( ((*phaseBound)[i].second.Np() <= np) && ((*phaseBound)[i-1].second.Np() > np) 
        || ((*phaseBound)[i-1].second.Np() <= np) && ((*phaseBound)[i].second.Np() > np)) {
      bracketNn.push_back((*phaseBound)[i].second.Nn());
      bracketNn.push_back((*phaseBound)[i-1].second.Nn());
    }  
  }

  if ( ((*phaseBound)[0].first.Np() <= np) && ((*phaseBound)[0].second.Np() > np) 
      || ((*phaseBound)[0].second.Np() <= np) && ((*phaseBound)[0].first.Np() > np)) {     
    bracketNn.push_back((*phaseBound)[0].first.Nn());
    bracketNn.push_back((*phaseBound)[0].second.Nn());
  }  

  if ( ((*phaseBound).back().first.Np() <= np) && ((*phaseBound).back().second.Np() > np) 
      || ((*phaseBound).back().second.Np() <= np) && ((*phaseBound).back().first.Np() > np)) {     
    bracketNn.push_back((*phaseBound).back().first.Nn());
    bracketNn.push_back((*phaseBound).back().second.Nn());
  }  

  if (bracketNn.size() < 1) { 
    if (mVerbose) std::cout << "# Next point found no bracketing proton fractions \n";
    return mpEos->FromNAndT(eosIn); 
  }

  std::sort(bracketNn.begin(), bracketNn.end()); 
  
  if (eosIn.Nn() < bracketNn[0] || eosIn.Nn() > bracketNn.back()) {
    if (mVerbose) std::cout << "# Next point outside of neutron brackets " 
        << bracketNn[0] << " " << bracketNn.back() << std::endl;
    return mpEos->FromNAndT(eosIn); 
  }
  
  // Ok, we are inside the phase boundary so find a guess for the up and down 
  // densities and fraction in each phase  
  std::vector<double> u, npArr;
  u.reserve(phaseBound->size()); 
  npArr.reserve(phaseBound->size()); 
  for (auto &p : *phaseBound) {
    double uC = (nb - p.first.Nb())/(p.second.Nb() - p.first.Nb());
    double npC = (1.0 - uC)*p.first.Np() + uC*p.second.Np(); 
    u.push_back(uC);
    npArr.push_back((1.0 - uC)*p.first.Np() + uC*p.second.Np()); 
  }

  // Find the bracketing points
  // The phaseBound array has been sorted in order of increasing mun, which 
  // means it is sorted in decreasing order of np 
  int up = 0;
  int lo = npArr.size()-1;
  while (up+1 < lo) {
    int mid = (lo + up)/2; 
    if (np > npArr[mid]) {
      lo = mid;
    } else { 
      up = mid;
    } 
  }
  
  double del = (np - npArr[lo]) / (npArr[up] - npArr[lo]);
  double nnLo = del * (*phaseBound)[up].first.Nn() 
      + (1.0 - del) * (*phaseBound)[lo].first.Nn(); 
  double npLo = del * (*phaseBound)[up].first.Np() 
      + (1.0 - del) * (*phaseBound)[lo].first.Np(); 
  double nnHi = del * (*phaseBound)[up].second.Nn() 
      + (1.0 - del) * (*phaseBound)[lo].second.Nn(); 
  double npHi = del * (*phaseBound)[up].second.Np() 
      + (1.0 - del) * (*phaseBound)[lo].second.Np(); 
  double uu = del * u[up] + (1.0 - del) * u[lo]; 
  
  if (mVerbose) std::cout << "# Next point is mixed phase \n";
  
  return GetState(eosIn, 
      EOSData::InputFromTNnNp(T, nnLo, npLo), 
      EOSData::InputFromTNnNp(T, nnHi, npHi), 
      uu);
}

std::vector<std::pair<EOSData, EOSData>>
GibbsPhaseConstruct::FindPhaseRange(double T, bool doMun, 
    double muStart, double muEnd, double deltaMu, double NLoG, double NHiG) const {
  
  std::vector<std::pair<EOSData, EOSData>> phasePoints;
   
  auto comp_func = [&](double a, double b) -> bool {
      if (muStart < muEnd) {return (a<b);} else {return (a>b);}};
   
  auto iter_func = [&](double dmu, bool adapt=false) -> void { 
    double mu = muStart;
    double dmu_init = dmu; 
    while (comp_func(mu, muEnd)) {
       
      std::pair<EOSData, EOSData> outDat;
      try { 
        outDat = FindPhasePoint(T, mu, NLoG, NHiG, doMun); 
      } catch (std::exception& e) {
        if (adapt && dmu>5.e-2 * dmu_init) {
          dmu = dmu*0.8;
        } else { 
          mu += dmu;
        } 
        continue;
      } catch (...) {
        if (adapt && dmu>5.e-2 * dmu_init) {
          dmu = dmu*0.8;
        } else { 
          mu += dmu;
        } 
        continue;
      }
      if (doMun && outDat.first.Np()*1.e-3 < outDat.first.Nn()) 
        phasePoints.push_back(outDat);
      if (!doMun && outDat.first.Nn()*1.e-3 < outDat.first.Np()) 
        phasePoints.push_back(outDat);
      
      if (doMun) { 
        NLoG = outDat.first.Np()*0.95;
        NHiG = outDat.second.Np()*1.05;
      } else { 
        NLoG = outDat.first.Nn()*0.95;
        NHiG = outDat.second.Nn()*1.05;
      }
      if (dmu < dmu_init) dmu = dmu*1.2;
      mu += dmu;
    }
    return;
  };
  
  // Do the basic iteration 
  iter_func(deltaMu); 

  // ok, that was a total failure 
  if (phasePoints.size() == 0) return phasePoints; 
  
  // Work forward with smaller steps from the last point 
  double muStartOld = muStart; 
  double muEndOld = muEnd; 
  if (doMun) {
    muStart = phasePoints.back().first.Mun();
  } else {
    muStart = phasePoints.back().first.Mup();
  }
  muEnd = muStart + 5.0 * deltaMu;
  iter_func(deltaMu*0.01); 
  muStart = muStartOld;
  muEnd = muEndOld;
   
  // Now work back towards muStart from the first succesful point 
  muEnd = muStart; 
  if (doMun) {
    muStart = phasePoints[0].first.Mun(); 
    NLoG = phasePoints[0].first.Np(); 
    NHiG = phasePoints[0].second.Np(); 
  } else { 
    muStart = phasePoints[0].first.Mup();
    NLoG = phasePoints[0].first.Nn(); 
    NHiG = phasePoints[0].second.Nn(); 
  } 
  if ((muStart-muEnd)/(muEnd + 1.e-4) < -0.02
      || (muStart-muEnd)/(muEnd + 1.e-4) > 0.02) 
      iter_func(-deltaMu*3.e-1, true);
  
  return phasePoints; 

} 

// Find a phase boundary for one temperature using boundary from another
// temperature
std::vector<std::pair<EOSData, EOSData>> 
GibbsPhaseConstruct::FindFixedTPhaseBoundary(double T, 
    std::vector<std::pair<EOSData, EOSData>> oldPhase) const {
  if (oldPhase.size() < 1) throw std::length_error("oldPhase contains no points.");
  double Told = oldPhase[0].first.T();
  std::vector<std::pair<EOSData, EOSData>> newPhase;
  for (auto& old : oldPhase) {
    std::pair<EOSData, EOSData> cur; 
    double Tc = T;
    for (bool doMun : {true, false} ) {
      cur = old;   
      try {
        for (double Tc = old.first.T(); Tc > T; Tc *= 0.9) {
          cur = FindPhasePoint(Tc, cur.first.Mun(), cur.first.Np(), 
              cur.second.Np(), doMun);
        }
        cur = FindPhasePoint(T, cur.first.Mun(), cur.first.Np(), 
            cur.second.Np(), doMun);
        newPhase.push_back(cur); 
      } catch(...) {}
      
      cur = old; 
      try {
        for (double Tc = old.first.T(); Tc > T; Tc *= 0.9) {
          cur = FindPhasePoint(Tc, cur.first.Mup(), cur.first.Nn(), 
              cur.second.Nn(), doMun);
        }
        cur = FindPhasePoint(T, cur.first.Mup(), cur.first.Nn(), 
            cur.second.Nn(), doMun);
        newPhase.push_back(cur); 
      } catch(...) {}
    }
  }
  
  std::sort (newPhase.begin(), newPhase.end(), 
      [](std::pair<EOSData, EOSData> a, std::pair<EOSData, EOSData> b) { 
      return (fabs(0.5-(a.second).Ye()) > fabs(0.5-(b.second).Ye()));});
  
  double YeMax = newPhase.back().second.Ye(); 
  YeMax = std::min(YeMax, 1.0-YeMax);
  std::cout << YeMax << std::endl; 
  auto selfBound = SelfBoundPoints(T, YeMax, 1.0 - YeMax);
  newPhase.insert(newPhase.end(), selfBound.begin(), selfBound.end());
   
  std::sort (newPhase.begin(), newPhase.end(), 
      [](std::pair<EOSData, EOSData> a, std::pair<EOSData, EOSData> b) { 
      return ((a.second).Mun() < (b.second).Mun());});
  
  return newPhase; 
}

std::vector<std::pair<EOSData, EOSData>>
GibbsPhaseConstruct::FindFixedTPhaseBoundary(double T, double NLoG, 
    double NHiG, double deltaMu) const {
  
  std::vector<std::pair<EOSData, EOSData>> phaseBound; 
  
  // Scan over chemical potentials 
  deltaMu = deltaMu * std::max(T, 2.0/197.3);
  double murange = 20.0 * std::max(T, 2.0/197.3);  
  auto phaseRange = FindPhaseRange(T, true, -murange, murange, deltaMu, NLoG, NHiG);
  phaseBound.insert(phaseBound.end(), phaseRange.begin(), phaseRange.end()); 
  
  phaseRange = FindPhaseRange(T, false, -murange, murange, deltaMu, NLoG, NHiG);
  phaseBound.insert(phaseBound.end(), phaseRange.begin(), phaseRange.end()); 
  
  // Sort the phase boundary data by the neutron chemical potential
  std::sort (phaseBound.begin(), phaseBound.end(), 
      [](std::pair<EOSData, EOSData> a, std::pair<EOSData, EOSData> b) { 
      return ((a.first).Mun() < (b.first).Mun());});
  
  // Now check for and remove (likely) bad points
  int i = 1; 
  while (i < phaseBound.size()) {
    double dnb = phaseBound[i-1].second.Nb()/phaseBound[i].second.Nb();
    if (dnb<0.80 || dnb>1.2) {
      phaseBound.erase(phaseBound.begin()+i);
    } else {
      i++;
    }
  }  

  return phaseBound; 
}

EOSData GibbsPhaseConstruct::GetState(const EOSData& eosIn, 
    const EOSData& lo, const EOSData& hi, double ug) const { 
  double T = eosIn.T();
  bool linear = true;
  double PScale = 1.e-20; 
  double MunScale = 1.e-20; 
  double MupScale = 1.e-20; 
  auto root_f = [this, T, &eosIn, &linear, &PScale, &MunScale, &MupScale]
      (std::vector<double> xx) -> std::vector<double> {
    
    EOSData eLo = mpEos->FromNAndT(
        EOSData::InputFromTNnNp(T, exp(xx[0]), exp(xx[1]))); 
    EOSData eHi = mpEos->FromNAndT(
        EOSData::InputFromTNnNp(T, exp(xx[2]) + exp(xx[0]), 
        exp(xx[3])+exp(xx[1])));
    if (linear) {
      return { (eHi.P() - eLo.P()),
           ((1.0 - xx[4])*eLo.Np() + xx[4]*eHi.Np())/eosIn.Np() - 1.0,
           ((1.0 - xx[4])*eLo.Nn() + xx[4]*eHi.Nn())/eosIn.Nn() - 1.0,
           (eHi.Mun() - eLo.Mun()),
           (eHi.Mup() - eLo.Mup())};
    } else {
      return { (eHi.P() - eLo.P()) / (PScale),
           ((1.0 - xx[4])*eLo.Np() + xx[4]*eHi.Np())/eosIn.Np() - 1.0,
           ((1.0 - xx[4])*eLo.Nn() + xx[4]*eHi.Nn())/eosIn.Nn() - 1.0,
           (eHi.Mun() - eLo.Mun()) / (MunScale),
           (eHi.Mup() - eLo.Mup()) / (MupScale)};
    }
  };
  
  MultiDimensionalRoot rootFinder = MultiDimensionalRoot(1.e-14, 200);
 
  std::vector<double> pars; 
  try {  
    pars = rootFinder(root_f, {log(lo.Nn()), log(lo.Np()), log(hi.Nn()-lo.Nn()), 
      log(hi.Np()-lo.Np()), ug}, 5);
  } catch (MultiDRootException& e) {
    pars = e.GetX();
  } catch (...) {
    std::cerr << "Gibbs linear solve failed." << std::endl;
    return mpEos->FromNAndT(
        EOSData::InputFromTNnNp(T, eosIn.Nn(), eosIn.Np())); 
  }
  
  EOSData eLo = mpEos->FromNAndT(
      EOSData::InputFromTNnNp(T, exp(pars[0]), exp(pars[1]))); 
  EOSData eHi = mpEos->FromNAndT(
      EOSData::InputFromTNnNp(T, exp(pars[2]) + exp(pars[0]), 
      exp(pars[3])+exp(pars[1])));
  rootFinder = MultiDimensionalRoot(1.e-10, 200);
  linear = false;
  
  try { 
    PScale = std::max(1.e0*fabs(eHi.P()), 1.e-6);
  } catch(...) {
    try { 
      PScale = std::max(1.e0*fabs(eLo.P()), 1.e-6);
    } catch(...) {
      if (mVerbose) std::cerr << "Pressure problems encountered, giving up." 
          << std::endl;
      return mpEos->FromNAndT(
          EOSData::InputFromTNnNp(T, eosIn.Nn(), eosIn.Np())); 
    }
  }

  MunScale = fabs(eHi.Mup());
  MupScale = fabs(eHi.Mun());
  try {
    pars = rootFinder(root_f, pars, 5);
    eLo = mpEos->FromNAndT(
        EOSData::InputFromTNnNp(T, exp(pars[0]), exp(pars[1]))); 
    eHi = mpEos->FromNAndT(
        EOSData::InputFromTNnNp(T, exp(pars[2]) + exp(pars[0]), 
        exp(pars[3])+exp(pars[1])));

  } catch (std::exception& e) {
    std::cerr << e.what(); 
    std::cerr << " in second try." << std::endl;
    //return mpEos->FromNAndT(
    //    EOSData::InputFromTNnNp(T, eosIn.Nn(), eosIn.Np())); 
  }

  if ((exp(pars[2]) < 1.e-1*exp(pars[0]) ) 
      && (exp(pars[3]) < 1.e-1*exp(pars[1]))) {
    if (mVerbose) std::cerr << "Failed to converge to separate points." 
        << std::endl;
    return mpEos->FromNAndT(
        EOSData::InputFromTNnNp(T, eosIn.Nn(), eosIn.Np())); 
  }
  
  if (pars[4]<0.0 || pars[4]>1.0) {
    return mpEos->FromNAndT(
        EOSData::InputFromTNnNp(T, eosIn.Nn(), eosIn.Np())); 
  }

  EOSData out = EOSData::Output(T, eosIn.Np(), eosIn.Nn(), 
      eHi.Mun(), eHi.Mup(), eHi.P(), 
      ((1.0 - pars[4])*eLo.S()*eLo.Nb() + pars[4]*eHi.S()*eHi.Nb())/eosIn.Nb(),
      ((1.0 - pars[4])*eLo.E()*eLo.Nb() + pars[4]*eHi.E()*eHi.Nb())/eosIn.Nb());
  out.SetPhases({eLo, eHi});
  return out;
} 

std::vector<std::pair<EOSData, EOSData>> 
GibbsPhaseConstruct::SelfBoundPoints(double T, 
    double yeMin, double yeMax) const {
  std::vector<EOSData> selfBound;
  // Calculate self bound points over a large range of Yes 
  for (double Ye = yeMin; Ye <= yeMax; Ye += 0.01) {
    try {
      selfBound.push_back(FindSelfBoundPoint(T, Ye)); 
    } catch(...) {} 
  }

  // Sort by Ye
  std::sort (selfBound.begin(), selfBound.end(), 
      [](EOSData a, EOSData b) { return (a.Ye() < b.Ye()); });
  std::vector<std::pair<EOSData, EOSData>> selfBoundPair;
  for (auto& pt : selfBound) {
    auto low = mpEos->FromMuAndT(
      EOSData::InputFromTMunMup(pt.T(), pt.Mun(), pt.Mup()));
    //auto lowTemp = EOSData::Output(pt.T(), 1.e-250, 1.e-250, 
    //    pt.Mun(), pt.Mup(), 0.0);
    if (low.size()>0) {
      if (low[0].Nb()<1.e-1) 
        selfBoundPair.push_back(std::pair<EOSData, EOSData>(low[0], pt)); 
    }
  }

  // Find matching points, ignoring pressure  
  return selfBoundPair;
} 

EOSData GibbsPhaseConstruct::FindSelfBoundPoint(double T, double Ye) const {
  // Get set of bracketing densities (How do we know the minimum is in there) 
  // Maybe try brute force search 
  double nbMin = 1.e5; 
  double nbMax = 0.2;
  double pMin = 1.e10;
  double pMax = mpEos->FromNAndT(EOSData::InputFromTNbYe(T, nbMax, Ye)).P();
  for (double lnb = -4.0; lnb< log10(nbMax); lnb += 0.1) { 
    double pc = 
        (mpEos->FromNAndT(EOSData::InputFromTNbYe(T, pow(10.0, lnb), Ye))).P();
    if (pc < pMin) {
      nbMin = pow(10.0, lnb);
      pMin = pc; 
    }
  }

  // Check that our two limits have different signs on the pressure 
  if (pMax*pMin>0.0) throw EOSException("No zero pressure.");  
  
  // Search for the zero, throw an error if it is not found
  auto PZero = [&T, &Ye, this](double nb) {
    return mpEos->FromNAndT(EOSData::InputFromTNbYe(T, nb, Ye)).P();
  };
  
  OneDimensionalRoot rootFinder(1.e-10, 150); 
  try {
    double nbZero = rootFinder(PZero, nbMin, nbMax);
    return mpEos->FromNAndT(EOSData::InputFromTNbYe(T, nbZero, Ye));
  } catch (...) { 
    throw EOSException("Zero pressure solver failed.");
  }
}
 
std::pair<EOSData, EOSData> GibbsPhaseConstruct::FindPhasePoint(double T, 
    double mu, double NLoG, double NHiG, bool doMun) const {
  
  if (NLoG>NHiG) {
    double tmp = NHiG;
    NHiG = NLoG;
    NLoG = tmp;  
  }

  // Declare pointers so we can switch between neutron and proton chemical
  // potentials 
  EOSData (EOSBase::*eosCall)(const EOSData&) const;
  double (EOSData::*getPotential)(void) const;
  EOSData (*eosDat)(double, double, double);
  if (doMun) {
    eosDat = EOSData::InputFromTNpMun;
    eosCall = &EOSBase::FromNpMunAndT;
    getPotential = &EOSData::Mup;
  } else {
    eosDat = EOSData::InputFromTNnMup;
    eosCall = &EOSBase::FromNnMupAndT;
    getPotential = &EOSData::Mun;
  }  
  
  // Using these guesses, search for the phase boundary points
  bool initial;
  double muScale = 1.0;
  double PScale = 1.0;
  auto root_f = [this, &mu, &T, initial, &eosCall, &eosDat, &doMun, 
      &getPotential, &muScale, &PScale]
      (std::vector<double> xx) -> std::vector<double> {
    EOSData eLo = (mpEos.get()->*eosCall)(eosDat(T, exp(xx[0]), mu)); 
    EOSData eHi = (mpEos.get()->*eosCall)(eosDat(T, exp(xx[0])+exp(xx[1]), mu));
    if (initial) {
      return {eHi.P() - eLo.P(), (eHi.*getPotential)() - (eLo.*getPotential)()}; 
    } else {  
      return {(eHi.P() - eLo.P())/(eHi.P() + PScale), 
          ((eHi.*getPotential)() - (eLo.*getPotential)()) 
          / muScale}; 
    }
  };
  
  // First get close with linear equations  
  MultiDimensionalRoot rootFinder(1.e-5*NHiG, 100);
  initial = true;
  std::vector<double> logN = rootFinder(root_f, 
      {log(NLoG), log(NHiG-NLoG)}, 2); 

  // Now get to high precision with non-linear scaled version of equations
  initial = false;
  auto eosLo = (mpEos.get()->*eosCall)(eosDat(T, exp(logN[0]), mu)); 
  auto eosHi = (mpEos.get()->*eosCall)(eosDat(T, exp(logN[0])+exp(logN[1]), mu));
  muScale = std::max(fabs(eosLo.Mun()), fabs(eosHi.Mup()));
  PScale = 1.e3*fabs(eosHi.P());
  rootFinder = MultiDimensionalRoot(1.e-13, 200);
  logN = rootFinder(root_f, logN, 2);
  
  eosLo = (mpEos.get()->*eosCall)(eosDat(T, exp(logN[0]), mu)); 
  eosHi = (mpEos.get()->*eosCall)(eosDat(T, exp(logN[0])+exp(logN[1]), mu));
  
  if (fabs(eosHi.Np()/eosLo.Np()-1.0) < 1.e-2 &&
      fabs(eosHi.Nn()/eosLo.Nn()-1.0) < 1.e-2) {
    throw std::runtime_error("GibbsPhaseConstruct converged to the same point.");
  }
  
  if (eosHi.Nb() > 1.e2) {
    throw std::runtime_error("GibbsPhaseConstruct converged to a huge density.");
  } 
  
  if (eosHi.Nb() < 1.e-3) {
    throw std::runtime_error("GibbsPhaseConstruct converged to a small density.");
  } 
  
  return std::pair<EOSData, EOSData>(eosLo, eosHi); 
}  

