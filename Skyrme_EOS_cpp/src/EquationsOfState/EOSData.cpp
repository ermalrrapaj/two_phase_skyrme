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
    mPhases(std::vector<EOSData>()) { 
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
    const double mun, const double mup, 
    const double pp, const double ss, const double ee,
    const double dpdnn, const double dpdnp, const double dpdt,
    const double dmundnn, const double dmundnp, const double dmundt,
    const double dmupdnn, const double dmupdnp, const double dmupdt,
    const double dsdnn, const double dsdnp, const double dsdt,
    const double dednn, const double dednp, const double dedt,
    const double mue) {
  EOSData out;
   
  out.mT.Set(T); 
  out.mMup.Set(mup); 
  out.mMun.Set(mun); 
  out.mNn.Set(nn); 
  out.mNp.Set(np);
  out.mMue.Set(mue);

  if (!(dpdnn!=dpdnn)) out.mP.SetDNn(dpdnn);
  if (!(dpdnp!=dpdnp)) out.mP.SetDNp(dpdnp);
  if (!(dpdt!=dpdt)) out.mP.SetDT(dpdt);
  if (!(dmundnn!=dmundnn)) out.mMun.SetDNn(dmundnn);
  if (!(dmundnp!=dmundnp)) out.mMun.SetDNp(dmundnp);
  if (!(dmundt!=dmundt)) out.mMun.SetDT(dmundt);
  if (!(dmupdnn!=dmupdnn)) out.mMup.SetDNn(dmupdnn);
  if (!(dmupdnp!=dmupdnp)) out.mMup.SetDNp(dmupdnp);
  if (!(dmupdt!=dmupdt)) out.mMup.SetDT(dmupdt); 
  if (!(dsdnn!=dsdnn)) out.mS.SetDNn(dsdnn);
  if (!(dsdnp!=dsdnp)) out.mS.SetDNp(dsdnp);
  if (!(dsdt!=dsdt)) out.mS.SetDT(dsdt);
  if (!(dednn!=dednn)) out.mE.SetDNn(dednn);
  if (!(dednp!=dednp)) out.mE.SetDNp(dednp);
  if (!(dedt!=dedt)) out.mE.SetDT(dedt);
  
  if (pp!=pp) throw std::range_error("Pressure is a NaN");
  out.mP.Set(pp);
  if (ss!=ss) throw std::range_error("Entropy is a NaN");
  out.mS.Set(ss);
  if (ee!=ee) throw std::range_error("Energy is a NaN");
  out.mE.Set(ee);
  
  return out; 
}


