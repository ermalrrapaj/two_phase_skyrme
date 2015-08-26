/// \file EOSBase.hpp
/// \author lroberts
/// \since Aug 18, 2015
///
/// \brief
///
///

#ifndef EOS_EOSBASE_HPP_
#define EOS_EOSBASE_HPP_

#include <memory>
#include <vector> 
 
#include "EOSData.hpp" 

/// Abstract base class for equations of state -
/// Provides a number of virtual functions for finding the state of the material
/// given the temperature and various combinations of the chemical potentials
/// and densities of neutrons and protons.
/// Additionally, any inheriting class is forced to include a method for 
/// creating copies via a unique pointer to simplify polymorphism.
class EOSBase {
public:
  /// Get the state of the material from the temperature and the neutron and 
  /// proton chemical potentials.  Since there may be multiple points in the 
  /// density temperature space which have a given set of chemical potentials
  /// the a vector of EOSData is returned.
  virtual std::vector<EOSData> FromMuAndT(const EOSData& eosIn) const =0; 

  /// Get the state of the material from the temperature, neutron density,
  /// and proton density.  This routine will also work if the input EOSData 
  /// has been set using the baryon number density and the electron fraction.
  virtual EOSData FromNAndT(const EOSData& eosIn) =0; 
  
  /// Get the state of the material from the temperature, neutron chemical 
  /// potential, and proton density.  This is guaranteed to be single 
  /// valued.
  virtual EOSData FromNpMunAndT(const EOSData& eosIn) const =0; 
  
  /// Get the state of the material from the temperature, proton chemical 
  /// potential, and neutron density.  This is guaranteed to be single 
  /// valued.
  virtual EOSData FromNnMupAndT(const EOSData& eosIn) const =0;
  
  /// Return a unique pointer to a copy of the current object. 
  virtual std::unique_ptr<EOSBase> MakeUniquePtr() const =0; 

protected:
private: 
};

#endif // EOS_EOSBASE_HPP_

