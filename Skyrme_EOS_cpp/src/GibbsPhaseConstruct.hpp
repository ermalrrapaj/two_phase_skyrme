/// \file GibbsPhaseConstruct.hpp
/// \authorr lroberts
/// \since Aug 18, 2015
///
/// \brief
///
///

#ifndef EOS_GIBBSPHASECONSTRUCT_HPP_
#define EOS_GIBBSPHASECONSTRUCT_HPP_

#include <memory> 

#include "EOSBase.hpp"

class GibbsPhaseConstruct {
public: 
  GibbsPhaseConstruct(const EOSBase& eos);  
protected: 
  std::unique_ptr<EOSBase> mpEos; 
}; 
#endif // EOS_GIBBSPHASECONSTRUCT_HPP_
