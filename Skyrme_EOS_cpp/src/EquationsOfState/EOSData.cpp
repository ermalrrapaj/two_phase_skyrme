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

#include "EquationsOfState/EOSData.hpp" 

EOSData::EOSData() : 
    mT(std::pair<double, bool>(0.0, false)), 
    mNp(std::pair<double, bool>(0.0, false)),
    mNn(std::pair<double, bool>(0.0, false)), 
    mMun(std::pair<double, bool>(0.0, false)),
    mMup(std::pair<double, bool>(0.0, false)), 
    mP(std::pair<double, bool>(0.0, false)), 
    mE(std::pair<double, bool>(0.0, false)),
    mS(std::pair<double, bool>(0.0, false)) {}
  
EOSData EOSData::InputFromTNnNp(const double T, const double nn, 
    const double np) {
  EOSData out;
  
  out.mT.first = T; 
  out.mT.second = true; 
  
  out.mNn.first = nn; 
  out.mNn.second = true; 

  out.mNp.first = np; 
  out.mNp.second = true; 

  return out; 
}

EOSData EOSData::InputFromSNnNp(const double S, const double nn, 
    const double np) {
  EOSData out;
  
  out.mS.first = S; 
  out.mS.second = true; 
  
  out.mNn.first = nn; 
  out.mNn.second = true; 

  out.mNp.first = np; 
  out.mNp.second = true; 

  return out; 
}

EOSData EOSData::InputFromTNbYe(const double T, const double nb, 
    const double ye) {
  EOSData out;
  
  out.mT.first = T; 
  out.mT.second = true; 
 
  out.mNn.first = nb * (1.0 - ye); 
  out.mNn.second = true; 

  out.mNp.first = nb * ye; 
  out.mNp.second = true; 

  return out; 
} 

EOSData EOSData::InputFromSNbYe(const double S, const double nb, 
    const double ye) {
  EOSData out;
  
  out.mS.first = S; 
  out.mS.second = true; 
 
  out.mNn.first = nb * (1.0 - ye); 
  out.mNn.second = true; 

  out.mNp.first = nb * ye; 
  out.mNp.second = true; 

  return out; 
} 

EOSData EOSData::InputFromTMunMup(const double T, const double mun, 
    const double mup) {
  EOSData out;
  
  out.mT.first = T; 
  out.mT.second = true; 
  
  out.mMun.first = mun; 
  out.mMun.second = true; 

  out.mMup.first = mup; 
  out.mMup.second = true; 

  return out; 
}
 
EOSData EOSData::InputFromTNpMun(const double T, const double np, 
    const double mun) {
  EOSData out;
  
  out.mT.first = T; 
  out.mT.second = true; 
  
  out.mMun.first = mun; 
  out.mMun.second = true; 

  out.mNp.first = np; 
  out.mNp.second = true; 

  return out; 
}

EOSData EOSData::InputFromTNnMup(const double T, const double nn, 
    const double mup) {
  EOSData out;
  
  out.mT.first = T; 
  out.mT.second = true; 
  
  out.mMup.first = mup; 
  out.mMup.second = true; 

  out.mNn.first = nn; 
  out.mNn.second = true; 

  return out; 
}

EOSData EOSData::Output(const double T, const double nn, const double np, 
    const double mun, const double mup, const double pp, 
    const double ss, const double ee) {
  EOSData out;
   
  out.mT.first = T; 
  out.mT.second = true; 
  
  out.mMun.first = mun; 
  out.mMun.second = true; 

  out.mMup.first = mup; 
  out.mMup.second = true; 
  
  out.mNn.first = nn; 
  out.mNn.second = true; 

  out.mNp.first = np; 
  out.mNp.second = true; 
 
  if (pp==pp) { 
    out.mP.first = pp; 
    out.mP.second = true; 
  }

  if (ss==ss) {
    out.mS.first = ss; 
    out.mS.second = true; 
  }

  if (ee==ee) {
    out.mE.first = ee; 
    out.mE.second = true; 
  }
  
  return out; 
}
 
double EOSData::T() const { 
  if (mT.second) { 
    return mT.first; 
  } else {
    throw std::logic_error("Temperature not set");
  }
}

double EOSData::Ye() const { 
  if (mNp.second && mNn.second) {
    return mNp.first/(mNn.first + mNp.first);
  } else {
    throw std::logic_error("Both densities not set");
  }
}

double EOSData::Nb() const { 
  if (mNp.second && mNn.second) {
    return mNn.first + mNp.first; 
  } else {
    throw std::logic_error("Both densities not set");
  }
} 

double EOSData::Nn() const { 
  if (mNn.second) {
    return mNn.first;
  } else { 
    throw std::logic_error("Neutron density not set");
  }
}

double EOSData::Np() const { 
  if (mNp.second) {
    return mNp.first;
  } else {
    throw std::logic_error("Proton density not set");
  }
}

double EOSData::P() const { 
  if (mP.second) { 
    return mP.first; 
  } else { 
    throw std::logic_error("Pressure not set");
  }
} 

double EOSData::Mun() const { 
  if (mMun.second) {
    return mMun.first;
  } else { 
    throw std::logic_error("Mun not set");
  }
} 

double EOSData::Mup() const { 
  if (mMup.second) {
    return mMup.first;
  } else {
    throw std::logic_error("Mup not set");
  }
} 

double EOSData::E() const { 
  if (mE.second) { 
    return mE.first; 
  } else {
    throw std::logic_error("Energy not set");
  }
} 
  
double EOSData::S() const { 
  if (mS.second) { 
    return mS.first; 
  } else {
    throw std::logic_error("Entropy not set");
  }
}


