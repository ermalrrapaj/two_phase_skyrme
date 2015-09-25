/// \file EOSData.cpp
/// \authorr lroberts
/// \since Aug 18, 2015
///
/// \brief
///
///

#include <utility>
#include <limits>
#include <exception>
#include <stdexcept>
#include <vector> 
#include <string> 

#include "EquationsOfState/EOSData.hpp" 

EOSData::EOSData() : 
    mT  (EOSDatum("Temperature")), 
    mNp (EOSDatum("Proton density")),
    mNn (EOSDatum("Neutron density")), 
    mMun(EOSDatum("Neutron chemical potential")),
    mMup(EOSDatum("Proton chemical potential")), 
    mMue(EOSDatum("Electron chemical potential")), 
    mP  (EOSDatum("Pressure")), 
    mE  (EOSDatum("Energy per baryon")),
    mS  (EOSDatum("Entropy per baryon")) {
  mVars["T"] = &mT;
  mVars["Np"] = &mNp;
  mVars["Nn"] = &mNn;
  mVars["Mup"] = &mMup;
  mVars["Mun"] = &mMun;
  mVars["Mue"] = &mMue;
  mVars["P"] = &mP;
  mVars["E"] = &mE;
  mVars["S"] = &mS;
}
  
EOSData EOSData::InputFromTNnNp(const double T, const double nn, 
    const double np) {
  EOSData out;
  out.mT.Set(T); 
  out.mNn.Set(nn); 
  out.mNp.Set(np);
  return out; 
}

EOSData EOSData::InputFromSNnNp(const double S, const double nn, 
    const double np) {
  EOSData out;
  out.mS.Set(S); 
  out.mNn.Set(nn); 
  out.mNp.Set(np); 
  return out; 
}

EOSData EOSData::InputFromTNbYe(const double T, const double nb, 
    const double ye) {
  EOSData out;
  out.mT.Set(T); 
  out.mNn.Set(nb*(1.0-ye));
  out.mNp.Set(nb*ye); 
  return out; 
} 

EOSData EOSData::InputFromSNbYe(const double S, const double nb, 
    const double ye) {
  EOSData out;
  out.mS.Set(S); 
  out.mNn.Set(nb*(1.0-ye));
  out.mNp.Set(nb*ye);
  return out; 
} 

EOSData EOSData::InputFromTMunMup(const double T, const double mun, 
    const double mup) {
  EOSData out;
  out.mT.Set(T); 
  out.mMun.Set(mun); 
  out.mMup.Set(mup); 
  return out; 
}
 
EOSData EOSData::InputFromTNpMun(const double T, const double np, 
    const double mun) {
  EOSData out;
  out.mT.Set(T); 
  out.mMun.Set(mun); 
  out.mNp.Set(np); 
  return out; 
}

EOSData EOSData::InputFromTNnMup(const double T, const double nn, 
    const double mup) {
  EOSData out;
  out.mT.Set(T); 
  out.mMup.Set(mup); 
  out.mNn.Set(nn); 
  return out; 
}

EOSData EOSData::Output(const double T, const double nn, const double np, 
    const double mun, const double mup, const double pp, 
    const double ss, const double ee) {
  EOSData out;
   
  out.mT.Set(T); 
  out.mMup.Set(mup); 
  out.mMun.Set(mun); 
  out.mNn.Set(nn); 
  out.mNp.Set(np); 
  
  if (pp==pp) out.mP.Set(pp);
  if (ss==ss) out.mS.Set(ss);
  if (ee==ee) out.mE.Set(ee);
  
  return out; 
}


