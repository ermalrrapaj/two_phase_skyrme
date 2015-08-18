/// \file EOSSkyrme.hpp
/// \author lroberts
/// \since Aug 18, 2015
///
/// \brief
///
///

#ifndef EOS_EOSSKYRME_HPP_
#define EOS_EOSSKYRME_HPP_

#include "EOSBase.hpp" 
//#include "EOSData.hpp"

class EOSSkyrme : public EOSBase {
public:
  EOSSkyrme();
  //EOSData FromMuAndT(const EOSData& eosIn) const; 
  EOSData FromNAndT(const EOSData& eosIn) const; 
  //EOSData FromNpMunAndT(const EOSData& eosIn) const;
   
protected:
  double mA, mB, mC, mD, mF, mG, mDelta; 
private: 
};

#endif // EOS_EOSSKYRME_HPP_

