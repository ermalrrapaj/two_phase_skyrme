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

#include "EquationsOfState/EOSBase.hpp"
#include "EquationsOfState/EOSData.hpp"
#include "Util/Constants.hpp" 

#include <boost/archive/text_iarchive.hpp> 
#include <boost/archive/text_oarchive.hpp> 

#include <boost/serialization/base_object.hpp> 
#include <boost/serialization/utility.hpp> 
#include <boost/serialization/list.hpp> 
#include <boost/serialization/assume_abstract.hpp> 

///
/// Finds the Gibbs phase boundaries for a given EoS.
///
class GibbsPhaseConstruct : public EOSBase {
public:
  /// Initialize with an EoS that has a non-convex region
  GibbsPhaseConstruct(const EOSBase& eos, bool findPhaseBound=true);  
  
  virtual EOSData FromNAndT(const EOSData& eosIn);
  
  /// Get the critical temperature for this EoS.
  double GetCriticalT() const {return mTCrit;}
  
  double GetMinimumT() const {return mTMin;}
  double GetMaximumT() const {return 200.0/Constants::HBCFmMeV;}
   
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
  
  /// Allows for easy serialization of the class so the phase boundary can 
  /// easily be read from file.  
  friend class boost::serialization::access; 
  template<class Archive> 
  void serialize(Archive & ar, const unsigned int /* File Version */) {
    ar & mTMin & mTMult & mTCrit & mVerbose & mPhaseBounds;
  }
    
protected:
  /// Solve the three Gibbs equilibrium equations and the constraint equations 
  /// for the neutron and proton density and return the mixed phase state
  EOSData GetState(const EOSData& eosIn, const EOSData& lo, const EOSData& hi, 
      double ug) const;
  
  /// Find a pair of phase points for a fixed temperature and chemical potential 
  std::pair<EOSData, EOSData> FindPhasePoint(double T, double mu, double NLoG, 
      double NHiG, bool doMun=true) const;
  
  /// Try to find the phase boundary over a range of chemical potentials.  This 
  /// routine is used by FindFixedTPHaseBoundary. 
  std::vector<std::pair<EOSData, EOSData>> FindPhaseRange(double T, bool doMun, 
      double muStart, double muEnd, double deltaMu, double NLoG, double NHiG) const;
  
  /// Calculate the entire phase boundary and store in mPhaseBounds 
  void FindPhaseBoundary();
  
  /// Copy of the input bulk EOS 
  std::unique_ptr<EOSBase> mpEos; 
  
  /// Vector of constant temperature phase boundaries that have already been 
  /// calculated 
  std::vector<std::vector<std::pair<EOSData, EOSData>>> mPhaseBounds;
  
  double mTMin;  ///< Minimum temperature calculated for the EoS 
  double mTMult; ///< How close together are the phase boundaries
  double mTCrit; ///< The critical temperature 
  bool mVerbose; ///< How verbose should I be?
}; 
#endif // EOS_GIBBSPHASECONSTRUCT_HPP_
