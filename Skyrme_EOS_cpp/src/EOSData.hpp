/// \file EOSData.hpp
/// \authorr lroberts
/// \since Aug 18, 2015
///
/// \brief
///
///

#ifndef EOS_EOSBASE_HPP_
#define EOS_EOSBASE_HPP_

#include <utility>
#include <limits>

class EOSData {

public:
  
  EOSData();
  
  static EOSData InputFromTNnNp(const double T, const double nn, const double np);
  static EOSData InputFromTNbYe(const double T, const double nb, const double ye);
  static EOSData InputFromTMunMup(const double T, const double mun, const double mup);
  static EOSData InputFromTMunNp(const double T, const double mun, const double np);
  static EOSData Output(const double T, const double nn, const double np, 
      const double mun, const double mup,  
      const double pp = std::numeric_limits<double>::quiet_NaN(), 
      const double ss = std::numeric_limits<double>::quiet_NaN(), 
      const double ee = std::numeric_limits<double>::quiet_NaN());
  
  double T();  
  double Ye(); 
  double Nb(); 
  double Nn();  
  double Np();  
  double P(); 
  double Mun();  
  double Mup();  
  double E();  
  double S();  

  std::pair<double, bool> mT, mNp, mNn;
  std::pair<double, bool> mP, mMun, mMup, mE, mS;

};
#endif // EOS_EOSBASE_HPP_

