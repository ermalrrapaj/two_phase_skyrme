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

///
/// Implements a uniform matter, Skyrme type EOS for with arbitrary parameters.
/// 
class EOSSkyrme : public EOSBase {
public:
  EOSSkyrme();
  EOSSkyrme(double A,  double B,  double C,  double D,  double F,
            double G,  double delta) : mA(A), mB(B), mC(C), mD(D),
            mF(F), mG(G), mDelta(delta) {};
  static EOSSkyrme FromErmalSkyrme(std::vector<double>& param);
  static EOSSkyrme FromSaturation(std::vector<double>& param);
  
  std::vector<EOSData> FromMuAndT(const EOSData& eosIn) const; 
  EOSData FromNAndT(const EOSData& eosIn); 
  EOSData FromNpMunAndT(const EOSData& eosIn) const;
  EOSData FromNnMupAndT(const EOSData& eosIn) const;
  
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

