/// \file EOSTestSuite.hpp
/// \author lroberts
/// \since Sep 30, 2015
///
/// \brief
///
///

#include "EquationsOfState/EOSTestSuite.hpp" 
#include "Util/OneDimensionalRoot.hpp"
#include "EquationsOfState/EOSData.hpp" 
#include "EquationsOfState/EOSBase.hpp" 

#define SMALL 1.e-50
#define COMPARE(x,y) fabs((x - y)/(x+SMALL)) 
#define COMPARES(x,y,z) fabs((x - y)/(z+SMALL)) 
 
int EOSTestSuite::CheckThermodynamicConsistency(double T, double nn, double np) const { 
  
  EOSData state = mpEos->FromNAndT(EOSData::InputFromTNnNp(T, nn, np));
  
  // Calculate derivatives in terms of baryon density and electron abundance
  double dmudt = 0.0; 
  double dmudn = 0.0;
  double muscale = 0.0; 
  try { 
    dmudt += state.dMupdT();
    dmudn += state.Ye()*state.dMupdNp() + (1.0 - state.Ye())*state.dMupdNn();
    muscale += state.Mup();
  } catch(...) {
  }
  try { 
    dmudt -= state.dMundT();
    dmudn -= state.Ye()*state.dMundNp() + (1.0 - state.Ye())*state.dMundNn();
    muscale += state.Mun();
  } catch(...) {
  }
  try { 
    dmudt += state.dMuedT();
    dmudn += state.Ye()*state.dMuedNp() + (1.0 - state.Ye())*state.dMuedNn();
    muscale += state.Mue();
  } catch(...) {
  }
  
  double dsdy = state.Nb()*state.dSdNp() - state.Nb()*state.dSdNn();
  double dsdn = state.Ye()*state.dSdNp() + (1.0 - state.Ye())*state.dSdNn();
  
  double dpdt = state.dPdT(); 
  double dpdy = state.Nb()*state.dPdNp() - state.Nb()*state.dPdNn();
  
  // Compare derivatives for consistency 
  double f1 = (dpdt + pow(state.Nb(),2)*dsdn)/(state.P()/state.T() + 1.e-50); 
  double f2 = (dpdy - pow(state.Nb(),2)*dmudn)/(state.P()/state.Ye() + 1.e-50); 
  double f3 = (dmudt + dsdy)/(muscale/state.T() + 1.e-50);  
  if (mVerbose) std::cout << "Thermodynamic Consistency " << f1 <<  " " << f2 
      << " " << f3 << std::endl;
  if (fabs(f1) + fabs(f2) + fabs(f3) < mTol) return 0;
  return 1; 
}

int EOSTestSuite::CheckAnalyticDerivatives(double T, double nn, double np) const { 
  double delta = 1.e-6; 
  EOSData state = mpEos->FromNAndT(EOSData::InputFromTNnNp(T, nn, np));
  double errtot = 0.0; 
  double nerr = 0.0;  
  auto CompFunc = [this, &errtot, &nerr](double nderiv, double deriv, 
      double scale, std::string name) {
    double err = (nderiv - deriv)/scale; 
    errtot += err;
    nerr++; 
    if (mVerbose) std::cout << name << " " << deriv/scale << " " << nderiv/scale 
        << " " << err << std::endl;  
  };

  // When the entropy is very low there seem to be problems with accuracy  
  double sscale = state.S() > 1.0 ? state.S() : 1.0; 
   
  // Check T derivatives  
  {  
    EOSData sp = mpEos->FromNAndT(EOSData::InputFromTNnNp(T*(1.0 + delta), nn, np));
    EOSData sm = mpEos->FromNAndT(EOSData::InputFromTNnNp(T*(1.0 - delta), nn, np));
    CompFunc((sp.P() - sm.P())/(sp.T() - sm.T()), state.dPdT(), state.P()/state.T(),"dPdT"); 
    CompFunc((sp.E() - sm.E())/(sp.T() - sm.T()), state.dEdT(), state.E()/state.T(),"dEdT"); 
    CompFunc((sp.S() - sm.S())/(sp.T() - sm.T()), state.dSdT(), sscale/state.T(),"dSdT"); 
    try { 
      CompFunc((sp.Mun() - sm.Mun())/(sp.T() - sm.T()), state.dMundT(), state.Mun()/state.T(),"dMundT"); 
    } catch(...) {} 
    try { 
      CompFunc((sp.Mup() - sm.Mup())/(sp.T() - sm.T()), state.dMupdT(), state.Mup()/state.T(),"dMupdT"); 
    } catch(...) {} 
    try { 
      CompFunc((sp.Mue() - sm.Mue())/(sp.T() - sm.T()), state.dMuedT(), state.Mue()/state.T(),"dMuedT"); 
    } catch(...) {} 
  } 
  
  // Check Nn derivatives  
  {  
    EOSData sp = mpEos->FromNAndT(EOSData::InputFromTNnNp(T, nn*(1.0 + delta), np));
    EOSData sm = mpEos->FromNAndT(EOSData::InputFromTNnNp(T, nn*(1.0 - delta), np));
    CompFunc((sp.P() - sm.P())/(sp.Nn() - sm.Nn()), state.dPdNn(), state.P()/state.Nn(),"dPdNn"); 
    CompFunc((sp.E() - sm.E())/(sp.Nn() - sm.Nn()), state.dEdNn(), state.E()/state.Nn(),"dEdNn"); 
    CompFunc((sp.S() - sm.S())/(sp.Nn() - sm.Nn()), state.dSdNn(), sscale/state.Nn(),"dSdNn"); 
    try { 
      CompFunc((sp.Mun() - sm.Mun())/(sp.Nn() - sm.Nn()), state.dMundNn(), state.Mun()/state.Nn(),"dMundNn"); 
    } catch(...) {} 
    try { 
      CompFunc((sp.Mup() - sm.Mup())/(sp.Nn() - sm.Nn()), state.dMupdNn(), state.Mup()/state.Nn(),"dMupdNn"); 
    } catch(...) {} 
    try { 
      CompFunc((sp.Mue() - sm.Mue())/(sp.Nn() - sm.Nn()), state.dMuedNn(), state.Mue()/state.Nn(),"dMuedNn"); 
    } catch(...) {} 
  } 
  
  // Check Np derivatives  
  {  
    EOSData sp = mpEos->FromNAndT(EOSData::InputFromTNnNp(T, nn, np*(1.0 + delta)));
    EOSData sm = mpEos->FromNAndT(EOSData::InputFromTNnNp(T, nn, np*(1.0 - delta)));
    CompFunc((sp.P() - sm.P())/(sp.Np() - sm.Np()), state.dPdNp(), state.P()/state.Np(),"dPdNp"); 
    CompFunc((sp.E() - sm.E())/(sp.Np() - sm.Np()), state.dEdNp(), state.E()/state.Np(),"dEdNp"); 
    CompFunc((sp.S() - sm.S())/(sp.Np() - sm.Np()), state.dSdNp(), sscale/state.Np(),"dSdNp"); 
    try { 
      CompFunc((sp.Mun() - sm.Mun())/(sp.Np() - sm.Np()), state.dMundNp(), state.Mun()/state.Np(),"dMundNp"); 
    } catch(...) {} 
    try { 
      CompFunc((sp.Mup() - sm.Mup())/(sp.Np() - sm.Np()), state.dMupdNp(), state.Mup()/state.Np(),"dMupdNp"); 
    } catch(...) {} 
    try { 
      CompFunc((sp.Mue() - sm.Mue())/(sp.Np() - sm.Np()), state.dMuedNp(), state.Mue()/state.Np(),"dMuedNp"); 
    } catch(...) {} 
  } 
  
  if (mVerbose) std::cout << errtot/nerr << " " << mTol << std::endl;
  if (errtot/nerr < mTol) return 0;
  return 1; 
}

int EOSTestSuite::CompressionTest(double Ye, double S) const {
  
  double nold, nnew, Told; 
   
  auto root_f = [this, &Ye, &nold, &nnew, &Told](double xx) -> double {
    EOSData eosc = mpEos->FromNAndT(EOSData::InputFromTNbYe(Told, nold, Ye));  
    EOSData eosn = mpEos->FromNAndT(EOSData::InputFromTNbYe(xx, nnew, Ye));  
    return (eosn.E() - eosc.E() 
        - 0.5*(eosn.P()/pow(eosn.Nb(), 2) + eosc.P()/pow(eosc.Nb(),2)) 
        * (eosn.Nb() - eosc.Nb()))/eosc.E();
  }; 
  
  double delta = 0.005; 
  nold = 1.e-4; 
  Told = mpEos->FromNAndS(EOSData::InputFromSNbYe(S, nold, Ye)).T();  
  
  OneDimensionalRoot rootFinder(1.e-10, 100);  
  for (double ln = log10(nold) + delta; ln < -0.5; ln += delta) {
    nnew = pow(10.0, ln); 
    Told = rootFinder(root_f,0.99*Told, 1.01*Told);
    if(mVerbose) {
      EOSData eosc = mpEos->FromNAndT(EOSData::InputFromTNbYe(Told, nnew, Ye));  
      std::cout << nnew << " " << Told << " " << eosc.S() << std::endl; 
    }
    nold = nnew; 
  }
  EOSData eosc = mpEos->FromNAndT(EOSData::InputFromTNbYe(Told, nnew, Ye));  

  if (fabs(eosc.S() - S) < 5.e-3) return 0;
  return 1; 
}

