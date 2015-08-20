/// \file EOSBase.hpp
/// \author lroberts
/// \since Aug 18, 2015
///
/// \brief
///
///

#ifndef EOS_EOSBASE_HPP_
#define EOS_EOSBASE_HPP_

#include <memory>
#include <vector> 
 
#include "EOSData.hpp" 

class EOSBase {
public:
  virtual std::vector<EOSData> FromMuAndT(const EOSData& eosIn) const =0; 
  virtual EOSData FromNAndT(const EOSData& eosIn) const =0; 
  virtual EOSData FromNpMunAndT(const EOSData& eosIn) const =0; 
  virtual EOSData FromNnMupAndT(const EOSData& eosIn) const =0; 
  virtual std::unique_ptr<EOSBase> MakeUniquePtr() const =0; 
protected:
private: 
};

#endif // EOS_EOSBASE_HPP_

