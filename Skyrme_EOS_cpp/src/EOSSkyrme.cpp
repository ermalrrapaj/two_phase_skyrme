/// \file EOSSkyrme.cpp
/// \author lroberts
/// \since Aug 18, 2015
///
/// \brief
///
///

#include "EOSSkyrme.hpp" 

static double HBC = 197.3; 

extern "C" {
  double ifermi12_(double* scale_density);
  double zfermi32_(double* eta);
}

EOSSkyrme::EOSSkyrme() :
    mA(-286.1/HBC), 
    mB(-107.1/HBC), 
    mC(968.0/HBC),
    mD(0.0), 
    mF(0.0), 
    mG(0.0),
    mDelta(2.002) {}

EOSData EOSSkyrme::FromNAndT(const EOSData& eosIn) const {
  double T = eosIn.T(); 
  double nn = eosIn.Nn();
  double np = eosIn.Np(); 
  
  double momsp = 1.0 + mF*(nn + np) + mG*(nn - np);
  double momsn = 1.0 + mF*(nn + np) - mG*(nn - np);
  
  //double etap = 
  //double etan = 

  //double taup = 0.0;
  //double taun = 0.0;

}
