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
class GibbsPhaseConstruct {
public:
  /// Initialize with an EoS that has a non-convex region
  GibbsPhaseConstruct(const EOSBase& eos);  
  
  /// Find the phase boundary in the np, nn plane for a fixed temperature 
  std::vector<std::pair<EOSData, EOSData>> FindFixedTPhaseBoundary(double T);
  
  /// Find a pair of phase points for a fixed temperature and chemical potential 
  std::vector<EOSData> FindPhasePoint(double T, double mu, double NLoG, 
      double NHiG, bool doMun=true);

protected:
  std::vector<std::pair<EOSData, EOSData>> FindPhaseRange(double T, bool doMun, 
      double muStart, double muEnd, double deltaMu, double NLoG, double NHiG);
  
  /// Copy of the input bulk EOS 
  std::unique_ptr<EOSBase> mpEos; 

}; 
#endif // EOS_GIBBSPHASECONSTRUCT_HPP_
