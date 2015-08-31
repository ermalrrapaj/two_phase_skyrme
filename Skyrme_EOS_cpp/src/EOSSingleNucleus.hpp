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
  EOSSingleNucleus(const EOSBase& eos, bool createPhaseBound=true
  ) : 
      GibbsPhaseConstruct(eos, createPhaseBound), mA0(56.0) {};
  
  EOSData FromNAndT(const EOSData& eosIn);

protected:
  std::vector<EOSData> EquilibriumConditions(const EOSData& eosIn, 
      const EOSData& eosLo, const EOSData& eosHi); 

  std::vector<double> DSurf(double u);
  
  double GetNQ(double u, double nn, double T); 
   
  double mA0;

};

#endif // EOS_EOSSINGLENUCLEUS_HPP_
