/// \file EOSSingleNucleus.cpp
/// \authorr lroberts
/// \since Aug 30, 2015
///
/// \brief
///
///

#include <math.h>
#include <iostream> 
#include <algorithm> 
#include <iterator> 

#include "EOSSingleNucleus.hpp"
#include "Constants.hpp"
#include "MultiDimensionalRoot.hpp"
#include "OneDimensionalRoot.hpp" 

EOSData EOSSingleNucleus::FromNAndT(const EOSData& eosIn) {
  
  EOSData gibbsState = GibbsPhaseConstruct::FromNAndT(eosIn);
  
  return gibbsState;
}

