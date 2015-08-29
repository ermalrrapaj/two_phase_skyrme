/// \file GibbsPhaseConstruct.hpp
/// \author lroberts
/// \since Aug 18, 2015
///
/// \brief
///
///

#ifndef EOS_GIBBSPHASECONSTRUCT_HPP_
#define EOS_GIBBSPHASECONSTRUCT_HPP_

#include <memory> 
#include <vector> 
#include <utility> 

#include "EOSBase.hpp"
#include "EOSData.hpp"

///
/// Finds the Gibbs phase boundaries for a given EoS.
///
class GibbsPhaseConstruct : public EOSBase {
public:
  /// Initialize with an EoS that has a non-convex region
  GibbsPhaseConstruct(const EOSBase& eos);  
  
  EOSData FromNAndT(const EOSData& eosIn);
  
  /// This function is not implemented and will throw an error if called 
  std::vector<EOSData> FromMuAndT(const EOSData& eosIn) const {
    throw std::logic_error("FromMuAndT has not been implemented.");
    return std::vector<EOSData>();
  } 
  
  /// This function is not implemented and will throw an error if called 
  EOSData FromNpMunAndT(const EOSData& eosIn) const {
    throw std::logic_error("FromNpMunAndT has not been implemented.");
    return EOSData();
  } 
  
  /// This function is not implemented and will throw an error if called 
  EOSData FromNnMupAndT(const EOSData& eosIn) const {
    throw std::logic_error("FromNnMupAndT has not been implemented.");
    return EOSData();
  } 
  
  std::unique_ptr<EOSBase> MakeUniquePtr() const {
    return std::unique_ptr<EOSBase>(new GibbsPhaseConstruct(*mpEos));
  } 
   
  /// Find a phase boundary in the np, nn plane for a fixed temperature 
  std::vector<std::pair<EOSData, EOSData>> FindFixedTPhaseBoundary(double T,
      double NLoG=1.e-20, double NHiG=0.08, double deltaMu=0.03) const;
  
protected:
  /// Solve the three Gibbs equilibrium equations and the constraint equations 
  /// for the neutron and proton density and return the mixed phase state
  EOSData GetState(const EOSData& eosIn, const EOSData& lo, const EOSData& hi, 
      double ug) const;
  
  /// Find a pair of phase points for a fixed temperature and chemical potential 
  std::pair<EOSData, EOSData> FindPhasePoint(double T, double mu, double NLoG, 
      double NHiG, bool doMun=true) const;

  std::vector<std::pair<EOSData, EOSData>> FindPhaseRange(double T, bool doMun, 
      double muStart, double muEnd, double deltaMu, double NLoG, double NHiG) const;
  
  /// Copy of the input bulk EOS 
  std::unique_ptr<EOSBase> mpEos; 
  
  /// Vector of constant temperature phase boundaries that have already been 
  /// calculated 
  std::vector<std::vector<std::pair<EOSData, EOSData>>> mPhaseBounds;
  
  double mTMult;
  double mTCrit; 
}; 
#endif // EOS_GIBBSPHASECONSTRUCT_HPP_
