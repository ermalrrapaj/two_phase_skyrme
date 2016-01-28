/// \file LDNucleus.hpp
/// \author lroberts
/// \since Sep 14, 2015
///
/// \brief
///
///

#ifndef EOS_LDNUCLEUS_HPP_
#define EOS_LDNUCLEUS_HPP_

#include <memory>
#include <vector> 
#include <iostream> 

#include "EquationsOfState/EOSData.hpp" 
#include "EquationsOfState/EOSBase.hpp"
#include "EquationsOfState/NucleusBase.hpp"
#include "Util/Constants.hpp" 

class LDNucleus : public NucleusBase { 
public:
  
  LDNucleus(int Z, int A, const EOSBase& eos);
  
  StaticNucleus GetStaticNucleus() const;

  virtual double GetVolume(const EOSData& eosIn, double ne) const;
  double GetBindingEnergy(const EOSData& eosIn, double ne) const;
  double GetBindingEnergy(const EOSData& eosIn, double ne, double v) const;
  double FreeEnergy(const EOSData& eosIn, double ne, double ni) const;
  double Entropy (const EOSData& eosIn, double ne, double ni) const;
  double NucleusPressure (const EOSData& eosIn, double ne, double uo) const;
  double Nucleusmup (const EOSData& eosIn, double ne, double uo, double ni) const;
  double Nucleusmun (const EOSData& eosIn, double ne, double uo, double ni) const;
 
  std::unique_ptr<NucleusBase> MakeUniquePtr() const { 
    return std::unique_ptr<NucleusBase>(new LDNucleus(GetZ(), 
      GetA(), *mpEos));
  }

protected:
  
  double SurfacePressure(double v) const;
  double SurfaceEnergy(double v) const;
  double CoulombPressure(double v, double npo, double ne) const;
  double CoulombEnergy(double v, double npo, double ne) const;
  
  virtual EOSData GetBulk(double T, double v) const { 
      return mpEos->FromNAndT(EOSData::InputFromTNnNp(T, 
      (double) NucleusBase::mN/v, (double) NucleusBase::mZ/v));
  }
    
  std::unique_ptr<EOSBase> mpEos;
  
  double mSs0; 
  double mSigma0; 
  double mV0;
    
  std::vector<std::pair<double,double>> PvsV;  

};

/// Should behave exactly like a liquid drop nucleus, except for with a constant 
/// volume
class StaticLDNucleus : public LDNucleus { 
public:  
  StaticLDNucleus(int Z, int A, const EOSBase& eos) : LDNucleus(Z, A, eos) {
    mBulk = LDNucleus::GetBulk(1.e-3/Constants::HBCFmMeV, mV0);
  }
  
  virtual double GetVolume(const EOSData& eosIn, double ne) const { return mV0;}
  
  std::unique_ptr<NucleusBase> MakeUniquePtr() const { 
    return std::unique_ptr<NucleusBase>(new StaticLDNucleus(GetZ(), GetA(), 
        *mpEos));
  }

protected:
  EOSData GetBulk(double T, double v) const { 
      return mBulk;
  }
  
  EOSData mBulk;

};

#endif // EOS_LDNUCLEUS_HPP_
