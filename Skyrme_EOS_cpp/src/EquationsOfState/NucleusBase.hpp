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
  virtual double GetBindingEnergy(const EOSData& eosIn, double ne, 
      double v) const =0; 
  virtual double GetVolume(const EOSData& eosIn, double ne) const =0; 
  virtual std::unique_ptr<NucleusBase> MakeUniquePtr() const =0;
   
  double GetN() const {return (double) mN;} 
  double GetZ() const {return (double) mZ;} 
  double GetA() const {return (double) mA;} 
    
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
  
  double GetBindingEnergy(const EOSData& eosIn, double ne, double /*v*/) const { 
    return GetBindingEnergy(eosIn, ne);
  }

  double GetVolume(const EOSData& /*eosIn*/, double /*ne*/) const {return 0.0;}
  
  std::unique_ptr<NucleusBase> MakeUniquePtr() const { 
    return std::unique_ptr<NucleusBase>(new FreeNucleus(*this));
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
  double GetVolume(const EOSData& eosIn, double ne) const;
  double GetBindingEnergy(const EOSData& eosIn, double ne) const;
  double GetBindingEnergy(const EOSData& eosIn, double ne, double v) const;

  std::unique_ptr<NucleusBase> MakeUniquePtr() const { 
    return std::unique_ptr<NucleusBase>(new LDNucleus(GetZ(), 
      GetA(), *mpEos));
  }

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
