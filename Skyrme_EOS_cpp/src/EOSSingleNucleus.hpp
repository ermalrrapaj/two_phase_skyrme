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

#include "EOSBase.hpp"
#include "EOSData.hpp"
#include "GibbsPhaseConstruct.hpp"

///
/// Calculate a two-phase EOS within the single nucleus 
/// approximation.
///
class EOSSingleNucleus : public GibbsPhaseConstruct {
public:
  EOSSingleNucleus(const EOSBase& eos) : GibbsPhaseConstruct(eos) {};
  
  EOSData FromNAndT(const EOSData& eosIn);

protected:
  std::vector<double> DSurf(double u);
};

#endif // EOS_EOSSINGLENUCLEUS_HPP_
