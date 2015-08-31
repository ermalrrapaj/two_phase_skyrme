/// \file EOSSingleNucleus.cpp
/// \authorr lroberts
/// \since Aug 30, 2015
///
/// \brief
///
///

#include <math.h>
#include <iostream> 
#include <algorithm> 
#include <iterator> 

#include "EOSSingleNucleus.hpp"
#include "Constants.hpp"
#include "MultiDimensionalRoot.hpp"
#include "OneDimensionalRoot.hpp" 

EOSData EOSSingleNucleus::FromNAndT(const EOSData& eosIn) {
  
  EOSData gibbsState = GibbsPhaseConstruct::FromNAndT(eosIn);
  
  return gibbsState;
}

std::vector<double> EOSSingleNucleus::DSurf(double u) {
  
  double a  = 1.0 - 1.5*pow(u, 1.0/3.0) + 0.5*u + 1.e-50;
  double uda = -0.5*(u, 1.0/3.0) + 0.5*u;
  double b = std::max(1.0 - 1.5*pow(1.0-u, 1.0/3.0) + 0.5*(1.0-u) + 1.e-50,0.0);
  double mudb = 0.5*pow(1.0-u, 1.0/3.0) - 0.5*(1.0-u);
  double denom = u*u + (1.0-u)*(1.0-u) + 0.6*u*u*pow(1.0-u, 2.0); 
  double ddeno = 2.0*u - 2.0*(1.0 - u) + 1.2*u*pow(1.0-u, 2) - 1.2*u*u*(1.0-u);  
  
  double D  = (1.0-u)*((1.0-u)*pow(a, 1.0/3.0) + u*pow(b, 1.0/3.0))/denom;
  
  double Dp = (1.0-2.0*u)*((1.0-u)*pow(a, 1.0/3.0) + u*pow(b, 1.0/3.0))/denom; 
  Dp +=(pow(1.0-u, 2)*pow(a, -2.0/3.0)/3.0*uda 
      + u*u*pow(b, -2.0/3.0)/3.0*mudb)/denom;
  Dp += u*(1.0-u)*(-pow(a, 1.0/3.0) + pow(b, 1.0/3.0))/denom;
  Dp -= u*(1.0-u)*((1.0-u)*pow(a, 1.0/3.0) 
      + u*pow(b, 1.0/3.0))/(denom*denom)*ddeno;
  return {D, Dp}; 
}

