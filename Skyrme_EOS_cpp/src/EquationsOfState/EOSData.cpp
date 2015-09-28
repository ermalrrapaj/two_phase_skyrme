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
    mS  (EOSDatum("Entropy per baryon")),
    mdPdNn (EOSDatum("Derivative of Pressure with respect to Neutron Density")),
    mdPdNp (EOSDatum("Derivative of Pressure with respect to Protron Density")),
    mdPdT (EOSDatum("Derivative of Pressure with respect to Temperature")),
    mdMundNn (EOSDatum("Derivative of Neutron chemical potential with respect to Neutron Density")),
    mdMundNp (EOSDatum("Derivative of Neutron chemical potential with respect to Protron Density")),
    mdMundT (EOSDatum("Derivative of Neutron chemical potential with respect to Temperature")),
    mdMupdNn (EOSDatum("Derivative of Protron chemical potential with respect to Neutron Density")),
    mdMupdNp (EOSDatum("Derivative of Protron chemical potential with respect to Protron Density")),
    mdMupdT (EOSDatum("Derivative of Protron chemical potential with respect to Temperature")), {
  mVars["T"] = &mT;
  mVars["Np"] = &mNp;
  mVars["Nn"] = &mNn;
  mVars["Mup"] = &mMup;
  mVars["Mun"] = &mMun;
  mVars["Mue"] = &mMue;
  mVars["P"] = &mP;
  mVars["E"] = &mE;
  mVars["S"] = &mS;
  mVars["dPdNn"] = &mdPdNn;
  mVars["dPdNp"] = &mdPdNp;
  mVars["dPdT"] = &mdPdT;
  mVars["dMundNn"] = &mdMundNn;
  mVars["dMundNp"] = &mdMundNp;
  mVars["dMundT"] = &mdMundT;
  mVars["dMupdNn"] = &mdMupdNn;
  mVars["dMupdNp"] = &mdMupdNp;
  mVars["dMupdT"] = &mdMupdT;
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
    const double mun, const double mup, 
    const double dpdnn, const double dpdnp, const double dpdt,
    const double dmundnn, const double dmundnp, const double dmundt,
    const double dmupdnn, const double dmupdnp, const double dmupdt,    
    const double pp, const double ss, const double ee) {
  EOSData out;
   
  out.mT.Set(T); 
  out.mMup.Set(mup); 
  out.mMun.Set(mun); 
  out.mNn.Set(nn); 
  out.mNp.Set(np);
  out.mdPdNn.Set(dpdnn);
  out.mdPdNp.Set(dpdnp);
  out.mdPdT.Set(dpdt);
  out.mdMundNn.Set(dmundnn);
  out.mdMundNp.Set(dmundnp);
  out.mdMundT.Set(dmundt);
  out.mdMupdNn.Set(dmupdnn);
  out.mdMupdNp.Set(dmupdnp);
  out.mdMupdT.Set(dmupdt); 
  
  if (pp==pp) out.mP.Set(pp);
  if (ss==ss) out.mS.Set(ss);
  if (ee==ee) out.mE.Set(ee);
  
  return out; 
}


