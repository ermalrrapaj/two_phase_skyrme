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
#include <string> 
#include <map> 

#include <boost/archive/text_iarchive.hpp> 
#include <boost/archive/text_oarchive.hpp> 

#include <boost/serialization/base_object.hpp> 
#include <boost/serialization/utility.hpp> 
#include <boost/serialization/list.hpp> 
#include <boost/serialization/vector.hpp> 
#include <boost/serialization/map.hpp> 
#include <boost/serialization/assume_abstract.hpp> 

#define STDNAN std::numeric_limits<double>::quiet_NaN()
#define UNINIT std::pair<bool, double>(false,0.0) 
class EOSDatum {
public:  
  EOSDatum() : name(" "), val(UNINIT), 
      dNp(UNINIT), dNn(UNINIT), dT(UNINIT), dNpdT(UNINIT), 
      dNndT(UNINIT), dNndNp(UNINIT) {} 
  
  EOSDatum(std::string name) : EOSDatum() {name = name;}
  
  EOSDatum(double valin, std::string namein="Un-named") : EOSDatum() {
    name = namein; 
    val = std::pair<bool, double>(true, valin); 
  } 
  
  double Get() const {
    if (val.first) return val.second;
    else throw std::logic_error(name + " not set."); 
  }

  double GetDT() const {
    if (dT.first) return dT.second;
    else throw std::logic_error(name + " dT not set."); 
  }

  double GetDNp() const {
    if (dNp.first) return dNp.second;
    else throw std::logic_error(name + " dNp not set."); 
  }

  double GetDNn() const {
    if (dNn.first) return dNn.second;
    else throw std::logic_error(name + " dNn not set."); 
  }

  void Set(const double& valin) {
    val.first = true;
    val.second = valin;
  }

  void SetDT(const double& valin) {
    dT.first = true;
    dT.second = valin;
  }

  void SetDNn(const double& valin) {
    dNn.first = true;
    dNn.second = valin;
  }

  void SetDNp(const double& valin) {
    dNp.first = true;
    dNp.second = valin;
  }

  friend class boost::serialization::access; 
  template<class Archive> 
  void serialize(Archive & ar, const unsigned int /* File Version */) {
    ar & name & val & dNp & dNn & dT & dNpdT & dNndT & dNndNp;
  }
   
protected: 
  std::string name;
  std::pair<bool, double> val, dNp, dNn, dT, dNpdT, dNndT, dNndNp;
};

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
  
  /// Create EOSData object from entropy per baryon and densities.
  static EOSData InputFromSNnNp(const double S, const double nn, const double np);
  
  /// Create EOSData object from entropy per baryon, baryon density, and Ye.
  static EOSData InputFromSNbYe(const double S, const double nb, const double ye);
  
  /// Create EOSData object from temperature and chemical potentials.
  static EOSData InputFromTMunMup(const double T, const double mun, const double mup);
  
  /// Create EOSData object from temperature, proton number density, and neutron 
  /// chemical potential.
  static EOSData InputFromTNpMun(const double T, const double np, const double mun);
  
  /// Create EOSData object from temperature, neutron number density, and proton 
  /// chemical potential.
  static EOSData InputFromTNnMup(const double T, const double nn, const double mun);
  
  /// Create EOSData object with most fields initialized.  Intended for output 
  /// from EOSBase calls.
  static EOSData Output(const double T, const double nn, const double np, 
      const double mun, const double mup,  
      const double pp = STDNAN, const double ss = STDNAN, const double ee = STDNAN,
      const double dpdnn=STDNAN, const double dpdnp=STDNAN, const double dpdt=STDNAN,
      const double dmundnn=STDNAN, const double dmundnp=STDNAN, const double dmundt=STDNAN,
      const double dmupdnn=STDNAN, const double dmupdnp=STDNAN, const double dmupdt=STDNAN,
      const double dsdnn=STDNAN, const double dsdnp=STDNAN, const double dsdt=STDNAN,
      const double dednn=STDNAN, const double dednp=STDNAN, const double dedt=STDNAN,
      const double mue = STDNAN);

  void SetPhases(std::vector<EOSData> phases) { mPhases = phases;} 
  
  /// 
  friend class boost::serialization::access; 
  template<class Archive> 
  void serialize(Archive & ar, const unsigned int /* File Version */) {
    ar & mT & mNp & mNn & mP & mMun & mMup & mMue & mE & mS & mPhases;
  }
   
  /// Get a vector of the subphases of this point
  std::vector<EOSData> Phases() const {return mPhases;}  
  double T()   const {return mT.Get();} ///< Return the temperature in [1/fm]
  double Ye()  const {return mNp.Get()/(mNp.Get() + mNn.Get() + 1.e-50);}; ///< Return the electron fraction 
  double Nb()  const {return mNp.Get() + mNn.Get();} ///< Return the baryon density [1/fm^3] 
  double Nn()  const {return mNn.Get();} ///< Return the neutron density [1/fm^3] 
  double Np()  const {return mNp.Get();} ///< Return the proton density [1/fm^3] 
  double P()   const {return mP.Get();}  ///< Return the pressure [1/fm^4] 
  double Mun() const {return mMun.Get();} ///< Return the neutron chemical potential [1/fm] 
  double Mup() const {return mMup.Get();} ///< Return the proton chemical potential [1/fm]
  double Mue() const {return mMue.Get();} ///< Return the proton chemical potential [1/fm]
  double E()   const {return mE.Get();} ///< Return the energy per baryon [1/fm] 
  double S()   const {return mS.Get();} ///< Return the entropy per baryon 
  double dPdNn()   const {return mP.GetDNn();} ///< Return the derivative of P with repsect to Nn [1/fm]
  double dPdNp()   const {return mP.GetDNp();} ///< Return the derivative of P with repsect to Np [1/fm]
  double dPdT()   const {return mP.GetDT();} ///< Return the derivative of P with repsect to T [1/fm^3]
  double dMundNn()   const {return mMun.GetDNn();} ///< Return the derivative of Mun with repsect to Nn [fm^2]
  double dMundNp()   const {return mMun.GetDNp();} ///< Return the derivative of Mun with repsect to Np [fm^2]
  double dMundT()   const {return mMun.GetDT();} ///< Return the derivative of Mun with repsect to T [1]
  double dMupdNn()   const {return mMup.GetDNn();} ///< Return the derivative of Mup with repsect to Nn [fm^2]
  double dMupdNp()   const {return mMup.GetDNp();} ///< Return the derivative of Mup with repsect to Np [fm^2]
  double dMupdT()   const {return mMup.GetDT();} ///< Return the derivative of Mup with repsect to T [1]
  double dMuedNn()   const {return mMue.GetDNn();} ///< Return the derivative of Mup with repsect to Nn [fm^2]
  double dMuedNp()   const {return mMue.GetDNp();} ///< Return the derivative of Mup with repsect to Np [fm^2]
  double dMuedT()   const {return mMue.GetDT();} ///< Return the derivative of Mup with repsect to T [1]
  double dSdNn () const {return mS.GetDNn();} /// ///< Return the derivative of S with repsect to Nn [1]
  double dSdNp () const {return mS.GetDNp();} /// ///< Return the derivative of S with repsect to Np [1]
  double dSdT () const {return mS.GetDT();} /// ///< Return the derivative of S with repsect to T [1/fm^2]
  double dEdNn () const {return mE.GetDNn();} /// ///< Return the derivative of E with repsect to Nn [1]
  double dEdNp () const {return mE.GetDNp();} /// ///< Return the derivative of E with repsect to Np [1]
  double dEdT () const {return mE.GetDT();} /// ///< Return the derivative of E with repsect to T [1/fm^2]
  
  void Set(const std::string name, const double val) { 
      try {mVars[name]->Set(val);} catch(...){std::cerr << name << " causes an error.\n";}}
  double Get(const std::string name) {return mVars[name]->Get();}
  
  void SetDNp(const std::string name, const double val) { 
      try {mVars[name]->SetDNp(val);} catch(...){std::cerr << name << " causes an error.\n";}}
  double GetDNp(const std::string name) {return mVars[name]->GetDNp();}
   
  void SetDNn(const std::string name, const double val) { 
      try {mVars[name]->SetDNn(val);} catch(...){std::cerr << name << " causes an error.\n";}}
  double GetDNn(const std::string name) {return mVars[name]->GetDNn();}
   
  void SetDT(const std::string name, const double val) {
      try {mVars[name]->SetDT(val);} catch(...){std::cerr << name << " causes an error.\n";}}
  double GetDT(const std::string name) {return mVars[name]->GetDT();}
   
protected:
  EOSDatum mT, mNp, mNn;
  EOSDatum mP, mMun, mMup, mMue, mE, mS;
  std::vector<EOSData> mPhases;
  std::map<std::string, EOSDatum*> mVars;

};
#endif // EOS_EOSDATA_HPP_

