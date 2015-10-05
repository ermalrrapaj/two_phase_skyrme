/// \file EOSSkyrme.hpp
/// \author lroberts
/// \since Aug 18, 2015
///
/// \brief
///
///

#ifndef EOS_EOSSKYRME_HPP_
#define EOS_EOSSKYRME_HPP_

#include "EquationsOfState/EOSBase.hpp" 
#include "Util/Constants.hpp"

///
/// Implements a uniform matter, Skyrme type EOS for with arbitrary parameters.
/// 
class EOSSkyrme : public EOSBase {
public:
  EOSSkyrme();
  EOSSkyrme(double A,  double B,  double C,  double D,  double F,
            double G,  double delta) : mA(A), mB(B), mC(C), mD(D),
            mF(F), mG(G), mDelta(delta) {};
  static EOSSkyrme FreeGas();
  static EOSSkyrme FromErmalSkyrme(const std::array<const double, 7>& param);
  static EOSSkyrme FromSaturation(const std::array<const double, 7>& param);
  
  std::vector<EOSData> FromMuAndT(const EOSData& eosIn) const; 
  EOSData FromNAndT(const EOSData& eosIn); 
  EOSData FromNpMunAndT(const EOSData& eosIn) const;
  EOSData FromNnMupAndT(const EOSData& eosIn) const;
  
  std::unique_ptr<EOSBase> MakeUniquePtr() const {
    return std::unique_ptr<EOSBase>(new EOSSkyrme(*this));
  }  
  
  double GetMinimumT() const { return 0.001/Constants::HBCFmMeV; }  
  double GetMaximumT() const { return 200.0/Constants::HBCFmMeV; }
    
  friend class boost::serialization::access; 
  template<class Archive> 
  void serialize(Archive & ar, const unsigned int /* File Version */) {
    ar & mA & mB & mC & mD & mF & mG & mDelta;
  }
   
protected:

  EOSData BaseEOSCall(const double T, const double nn, const double np);
  void BaseEOSCallD1( EOSData& eosout );
  double mA;
  double mB;
  double mC;
  double mD;
  double mF;
  double mG;
  double mDelta;
   
private: 
  double G[2]; //{Gn, Gp};
  double mom[2]; //{momsn, momsp};
  double dm[2][3];//{{dmndnn, dmndnp, dmndT}, {dmpdnn, dmpdnp, dmpdT}};
  double eta[2];//{etan, etap}; 
  double deta[2][3];//{{detandnn, detandnp, detandT}, {detapdnn, detapdnp, detapdT}};
  double tau[2];//{taun, taup};
  double dtau[2][3];//{{dtaundnn, dtaundnp, dtaundT}, {dtaupdnn, dtaupdnp, dtaupdT}};
  double U[2];//{Un, Up};
  double dU[2][3];//{{dUndnn, dUndnp, dUndT}, {dUpdnn, dUpdnp, dUpdT}};
};

#endif // EOS_EOSSKYRME_HPP_

