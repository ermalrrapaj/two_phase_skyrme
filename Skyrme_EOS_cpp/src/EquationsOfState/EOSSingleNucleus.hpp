/// \file EOSSingleNucleus.hpp
/// \author lroberts
/// \since Aug 29, 2015
///
/// \brief
///
///

#ifndef EOS_EOSSINGLENUCLEUS_HPP_
#define EOS_EOSSINGLENUCLEUS_HPP_

#include <memory> 
#include <vector> 
#include <utility> 

#include "EquationsOfState/EOSBase.hpp"
#include "EquationsOfState/EOSData.hpp"
#include "EquationsOfState/GibbsPhaseConstruct.hpp"

///
/// Calculate a two-phase EOS within the single nucleus 
/// approximation.
///
class EOSSingleNucleus : public GibbsPhaseConstruct {
public:
  EOSSingleNucleus(const EOSBase& eos, bool createPhaseBound = true, 
      double TMeVMin = 0.1, double TMeVMax = 2.e2) : 
      GibbsPhaseConstruct(eos, createPhaseBound, TMeVMin, TMeVMax), 
      mA0(60.0),
      mSigma0(1.15/Constants::HBCFmMeV),
      mSs0(45.8/Constants::HBCFmMeV), 
      mKs0(220.0/Constants::HBCFmMeV){};
  
  EOSData FromNAndT(const EOSData& eosIn);

protected:
  std::vector<EOSData> EquilibriumConditions(const EOSData& eosIn, 
      const EOSData& eosLo, const EOSData& eosHi, double llam0 = -5.0); 

  EOSData GetStateFromPhases(double u, const EOSData& low, const EOSData& high);
  double GetMuh(double u, double nn, double T); 
  std::array<double, 2> DSurf(double u);
  std::array<double, 2> Sigma(double xp);
  std::array<double, 3> HFunc(double T, double xp);
   
  double mA0;
  double mSigma0;
  double mSs0;
  double mKs0; 
  std::vector<double> mLastSet;
};

#endif // EOS_EOSSINGLENUCLEUS_HPP_
