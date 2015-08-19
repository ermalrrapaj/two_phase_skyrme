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
  EOSData FromNpMunAndT(const EOSData& eosIn) const;
  std::unique_ptr<EOSBase> MakeUniquePtr() const {
    return std::unique_ptr<EOSBase>(new EOSSkyrme(*this));
  }   
protected:

  EOSData BaseEOSCall(const double T, const double nn, const double np) const;
  double mA;
  double mB;
  double mC;
  double mD;
  double mF;
  double mG;
  double mDelta;
   
private: 
};

#endif // EOS_EOSSKYRME_HPP_

