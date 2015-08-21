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

#include "EOSBase.hpp"
#include "EOSData.hpp"

///
/// Finds the Gibbs phase boundaries for a given EoS.
///
class GibbsPhaseConstruct {
public:
 
  GibbsPhaseConstruct(const EOSBase& eos);  
  
  std::vector<EOSData> FindPhasePoint(double T, double mu, double NLoG, 
      double NHiG, bool doMun=true);

protected: 
  std::unique_ptr<EOSBase> mpEos; 

}; 
#endif // EOS_GIBBSPHASECONSTRUCT_HPP_
