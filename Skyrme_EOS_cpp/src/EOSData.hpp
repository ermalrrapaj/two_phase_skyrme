/// \file EOSData.hpp
/// \authorr lroberts
/// \since Aug 18, 2015
///
/// \brief
///
///

#ifndef EOS_EOSDATA_HPP_
#define EOS_EOSDATA_HPP_

#include <utility>
#include <limits>
#include <vector>

#include <boost/archive/text_iarchive.hpp> 
#include <boost/archive/text_oarchive.hpp> 

#include <boost/serialization/base_object.hpp> 
#include <boost/serialization/utility.hpp> 
#include <boost/serialization/list.hpp> 
#include <boost/serialization/vector.hpp> 
#include <boost/serialization/assume_abstract.hpp> 

///
/// Class for storing state data passed to or returned from equations of state.
/// Performs various error checking to make sure that values being asked for are 
/// have been set.  If the variable asked for is not set, a std::logic_error is 
/// thrown.  This class provides a uniform, extensible interface to 
/// EoSs such that new variables can easily be added.  All inputs should assume
/// hbar = c = 1 and be expressed in terms of fm.   
///
class EOSData {
public:
  
  EOSData();
  
  /// Create EOSData object from temperature and densities.
  static EOSData InputFromTNnNp(const double T, const double nn, const double np);
  
  /// Create EOSData object from temperature, baryon density, and Ye.
  static EOSData InputFromTNbYe(const double T, const double nb, const double ye);
  
  /// Create EOSData object from temperature and chemical potentials.
  static EOSData InputFromTMunMup(const double T, const double mun, const double mup);
  
  /// Create EOSData object from temperature, proton number density, and neutron 
  /// chemical potential.
  static EOSData InputFromTNpMun(const double T, const double np, const double mun);
  
  /// Create EOSData object from temperature, neutron number density, and proton 
  /// chemical potential.
  static EOSData InputFromTNnMup(const double T, const double nn, const double mup);
  
  /// Create EOSData object with most fields initialized.  Intended for output 
  /// from EOSBase calls.
  static EOSData Output(const double T, const double nn, const double np, 
      const double mun, const double mup,  
      const double pp = std::numeric_limits<double>::quiet_NaN(), 
      const double ss = std::numeric_limits<double>::quiet_NaN(), 
      const double ee = std::numeric_limits<double>::quiet_NaN());

  void SetPhases(std::vector<EOSData> phases) { mPhases = phases;}
  
  void SetT(double T) {mT=std::pair<double, bool>(T,true);} 
  
  /// 
  friend class boost::serialization::access; 
  template<class Archive> 
  void serialize(Archive & ar, const unsigned int /* File Version */) {
    ar & mT & mNp & mNn & mP & mMun & mMup & mE & mS & mPhases;
  }
   
  /// Get a vector of the subphases of this point
  std::vector<EOSData> Phases() const {return mPhases;}  
  double T() const; ///< Return the temperature in [1/fm]
  double Ye() const; ///< Return the electron fraction 
  double Nb() const; ///< Return the baryon density [1/fm^3] 
  double Nn() const; ///< Return the neutron density [1/fm^3] 
  double Np() const; ///< Return the proton density [1/fm^3] 
  double P() const;  ///< Return the pressure [1/fm^4] 
  double Mun() const; ///< Return the neutron chemical potential [1/fm] 
  double Mup() const; ///< Return the proton chemical potential [1/fm]
  double E() const; ///< Return the energy per baryon [1/fm] 
  double S() const; ///< Return the entropy per baryon 

protected:
  std::vector<EOSData> mPhases;
  std::pair<double, bool> mT, mNp, mNn;
  std::pair<double, bool> mP, mMun, mMup, mE, mS;

};
#endif // EOS_EOSDATA_HPP_

