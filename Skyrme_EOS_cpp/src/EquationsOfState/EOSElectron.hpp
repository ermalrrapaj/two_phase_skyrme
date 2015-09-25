/// \file EOSElectron.hpp
/// \author lroberts
/// \since Sep 23, 2015
///
/// \brief
///
///

#ifndef EOS_EOSELECTRON_HPP_
#define EOS_EOSELECTRON_HPP_

#include <memory>
#include <vector> 
#include <iostream> 
 
#include "EquationsOfState/EOSData.hpp" 
#include "EquationsOfState/EOSBase.hpp" 
#include "Util/OneDimensionalRoot.hpp"


class EOSElectron : public EOSBase { 

public:
   
  EOSData FromNAndT(const EOSData& eosIn);
  
  double GetMinimumT() const {return 0.001/Constants::HBCFmMeV;}
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
    return std::unique_ptr<EOSBase>(new EOSElectron());
  } 

protected: 
};

#endif // EOS_EOSELECTRON_HPP_
