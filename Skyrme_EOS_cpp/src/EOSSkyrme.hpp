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
  
  EOSSkyrme(double A,  double B,  double C,  double D,  double F,
            double G,  double delta) : mA(A), mB(B), mC(C), mD(D),
            mF(F), mG(G), mDelta(delta) {};
  
  static EOSSkyrme FromErmalSkyrme(const double a, const double b,
      const double t0, const double x0, const double t3, 
      const double x3, const double alpha);
  
  static EOSSkyrme FromErmalSkyrme(std::vector<double>& param);
  
  static EOSSkyrme FromEosObs(const double ns, const double BE, const double K,
      const double Sv, const double mstarom, const double L, const double Ks);
  
  
  std::vector<EOSData> FromMuAndT(const EOSData& eosIn) const; 
  EOSData FromNAndT(const EOSData& eosIn) const; 
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

