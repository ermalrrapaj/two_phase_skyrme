/// \file EOSNSE.hpp
/// \author lroberts
/// \since Sep 15, 2015
///
/// \brief
///
///

#ifndef EOS_EOSNSE_HPP_
#define EOS_EOSNSE_HPP_

#include "EquationsOfState/EOSBase.hpp" 
#include "EquationsOfState/NucleusBase.hpp" 
#include "Util/Constants.hpp"

class EOSNSE : public EOSBase {
public:
  EOSNSE(std::vector<std::unique_ptr<NucleusBase>> nuclei,
      const EOSBase& eos) : 
      mTMin(0.1/Constants::HBCFmMeV), 
      mpEos(eos.MakeUniquePtr()) { 
    for (auto & nuc : nuclei) 
      mNuclei.push_back(nuc->MakeUniquePtr());
  }
  
  EOSNSE(const EOSNSE& other); 
   
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
  
  EOSData FromNAndT(const EOSData& eosIn);

  double GetMinimumT() const {return mTMin;}
  double GetMaximumT() const {return 200.0/Constants::HBCFmMeV;}

  std::unique_ptr<EOSBase> MakeUniquePtr() const {
    return std::unique_ptr<EOSBase>(new EOSNSE(*this));
  }  

private: 
  std::vector<std::unique_ptr<NucleusBase>> mNuclei; 
  double mTMin;
  std::shared_ptr<EOSBase> mpEos; 

};

#endif // EOS_EOSSKYRME_HPP_
