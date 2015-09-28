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
    mdMupdT (EOSDatum("Derivative of Protron chemical potential with respect to Temperature")),
    mdSdNn (EOSDatum("Derivative of Entropy per baryon with respect to Neutron Density")),
    mdSdNp (EOSDatum("Derivative of Entropy per baryon with respect to Protron Density")),
    mdSdT (EOSDatum("Derivative of Entropy per baryon with respect to Temperature")) {
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
  mVars["dSdNn"] = &mdSdNn;
  mVars["dSdNp"] = &mdSdNp;
  mVars["dSdT"] = &mdSdT;
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
    const double pp, const double ss, const double ee,
    const double dpdnn, const double dpdnp, const double dpdt,
    const double dmundnn, const double dmundnp, const double dmundt,
    const double dmupdnn, const double dmupdnp, const double dmupdt,
    const double dsdnn, const double dsdnp, const double dsdt) {
  EOSData out;
   
  out.mT.Set(T); 
  out.mMup.Set(mup); 
  out.mMun.Set(mun); 
  out.mNn.Set(nn); 
  out.mNp.Set(np);

  if (!(dpdnn!=dpdnn)) out.mdPdNn.Set(dpdnn);
  if (!(dpdnp!=dpdnp)) out.mdPdNp.Set(dpdnp);
  if (!(dpdt!=dpdt)) out.mdPdT.Set(dpdt);
  if (!(dmundnn!=dmundnn)) out.mdMundNn.Set(dmundnn);
  if (!(dmundnp!=dmundnp)) out.mdMundNp.Set(dmundnp);
  if (!(dmundt!=dmundt)) out.mdMundT.Set(dmundt);
  if (!(dmupdnn!=dmupdnn)) out.mdMupdNn.Set(dmupdnn);
  if (!(dmupdnp!=dmupdnp)) out.mdMupdNp.Set(dmupdnp);
  if (!(dmupdt!=dmupdt)) out.mdMupdT.Set(dmupdt); 
  if (!(dsdnn!=dsdnn)) out.mdSdNn.Set(dsdnn);
  if (!(dsdnp!=dsdnp)) out.mdSdNp.Set(dsdnp);
  if (!(dsdt!=dsdt)) out.mdSdT.Set(dsdt);
  
  if (!(pp!=pp)) out.mP.Set(pp);
  if (!(ss!=ss)) out.mS.Set(ss);
  if (!(ee!=ee)) out.mE.Set(ee);
  
  return out; 
}


