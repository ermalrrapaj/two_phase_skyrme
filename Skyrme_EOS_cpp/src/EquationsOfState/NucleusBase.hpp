/// \file NucleusBase.hpp
/// \author lroberts
/// \since Sep 14, 2015
///
/// \brief
///
///

#ifndef EOS_NUCLEUSBASE_HPP_
#define EOS_NUCLEUSBASE_HPP_

#include <memory>
#include <vector> 
#include <iostream> 

#include "EquationsOfState/EOSData.hpp" 
#include "EquationsOfState/EOSBase.hpp"
#include "Util/Constants.hpp" 

class NucleusBase {
public:
  NucleusBase(int Z, int A) : mZ(Z), mA(A), mN(A-Z) {}
  virtual double GetBindingEnergy(const EOSData& eosIn, double ne) const =0; 
  //virtual double GetVolume(const EOSData& eosIn) const =0; 
  
protected:
  int mZ, mN, mA;

}; 

class FreeNucleus : public NucleusBase {
public:
  FreeNucleus(int Z, int A, double BE, const std::vector<double>& TGrid, 
      const std::vector<double>& PFGrid) : NucleusBase(Z, A), mBE(BE), 
      mTg(TGrid), mPFg(PFGrid) {}

  /// \todo Still need to include partition function in binding energy
  double GetBindingEnergy(const EOSData& /*eosIn*/, double /*ne*/) const { 
    return mBE;
  }
   
protected: 
  double mBE;
  std::vector<double> mTg, mPFg;

}; 

class LDNucleus : public NucleusBase { 
public:
  LDNucleus(int Z, int A, const EOSBase& eos) : NucleusBase(Z, A), 
      mpEos(eos.MakeUniquePtr()),
      mSigma0(1.15/Constants::HBCFmMeV),
      mSs0(45.8/Constants::HBCFmMeV) {} 
  double GetBindingEnergy(const EOSData& eosIn, double ne) const;

protected:
  double SurfacePressure(double v) const;
  double SurfaceEnergy(double v) const;
  double CoulombPressure(double v, double npo, double ne) const;
  double CoulombEnergy(double v, double npo, double ne) const;
  std::unique_ptr<EOSBase> mpEos;
  double mSs0; 
  double mSigma0; 

};

#endif // EOS_NUCLEUSBASE_HPP_
