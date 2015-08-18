/// \file EOSBase.hpp
/// \author lroberts
/// \since Aug 18, 2015
///
/// \brief
///
///

#ifndef EOS_EOSBASE_HPP_
#define EOS_EOSBASE_HPP_

#include "EOSData.hpp" 

class EOSBase {
public:
  EOSBase(){}
  //EOSData FromMuAndT(const EOSData& eosIn) const =0; 
  virtual EOSData FromNAndT(const EOSData& eosIn) const =0; 
  //EOSData FromMixedAndT(const EOSData& eosIn) const =0; 
protected:
private: 
};

#endif // EOS_EOSBASE_HPP_

