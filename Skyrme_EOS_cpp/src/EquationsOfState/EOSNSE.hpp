/// \file EOSNSE.hpp
/// \author lroberts
/// \since Sep 15, 2015
///
/// \brief Solves equations of NSE with an excluded volume correction and 
/// calculates the equation of state for this distribution of nuclei.
///
///

#ifndef EOS_EOSNSE_HPP_
#define EOS_EOSNSE_HPP_

#include "EquationsOfState/EOSBase.hpp" 
#include "EquationsOfState/NucleusBase.hpp" 
#include "Util/Constants.hpp"

/// 
/// For a reaction $i \leftrightarrow j + k$, it is easy to show that the principle
/// of detailed balance implies that
/// \begin{equation*}
/// \mu_i = \mu_j + \mu_k
/// \end{equation*}
/// holds when the forward rate balances the backward rate.  When all strong
/// interactions are in equilibrium with their inverses, invoking detailed balance
/// for all of the reactions implies that
/// \begin{equation*}
/// \mu_i = N_i \mu_n + Z_i \mu_p.
/// \end{equation*}
/// where $N_i$ is the number of neutrons and $Z_i$ is the number of protons in a
/// nucleus of species $i$.
/// This of course implies that if the neutron and proton number densities are
/// known, the number densities of all other species in the medium are known.
/// Assuming that the heavy nuclei obey Boltzmann statistics (i.e. $\mu_i = m_i + T
/// \ln \left[n_i/(G_i n_Q)\right]$ gives
/// \begin{equation*}
/// n_i = A_i^{3/2} G_i(T) n_Q \exp\left(\left[Z_i \mu_p + N_i \mu_n
/// - m_i \right]/T\right),
/// \end{equation*}
/// where $G_i$ is the temperature dependent internal partition function of species
/// $i$.  When neutrons and protons also obey Boltzmann statistics, the number
/// density can be expressed as
/// \begin{equation*}
/// n_i = \frac{A_i^{3/2} G_i(T)}{2^{A_i}} n_Q^{1-A_i} n_n^{N_i} n_p^{Z_i}
/// \exp(BE_i/T), 
/// \end{equation*} 
/// where $BE_i$ is the binding energy of species $i$. 
/// 
/// Rather than being given the free proton and neutron densities, we often
/// only know the total neutron and proton densities.  In that case the equations of
/// neutron and proton number conservation must be solved as functions of the free
/// proton and neutron number densities
/// \begin{eqnarray*}
/// &&\sum_i N_i n_i(n_{n,o}, n_{p,o}) = n_n, \\
/// &&\sum_i Z_i n_i(n_{n,o}, n_{p,o}) = n_p.
/// \end{eqnarray*}
/// 

class EOSNSE : public EOSBase {
public:
  EOSNSE(const std::vector<std::unique_ptr<NucleusBase>>& nuclei,
      const EOSBase& eos) : 
      mTMin(0.1/Constants::HBCFmMeV), 
      mpEos(eos.MakeUniquePtr()) { 
    for (auto& nuc : nuclei) 
      mNuclei.push_back(nuc->MakeUniquePtr());
  }
  
  EOSNSE(const EOSNSE& other); 
   
  /// This function is not implemented and will throw an error if called 
  std::vector<EOSData> FromMuAndT(const EOSData& eosIn) const {
    throw std::logic_error("FromMuAndT has not been implemented.");
    return std::vector<EOSData>();
  } 
  
  /// This function is not implemented and will throw an error if called 
  EOSData FromNpMunAndT(const EOSData& eosIn) const {
    throw std::logic_error("FromNpMunAndT has not been implemented.");
    return EOSData();
  } 
  
  /// This function is not implemented and will throw an error if called 
  EOSData FromNnMupAndT(const EOSData& eosIn) const {
    throw std::logic_error("FromNnMupAndT has not been implemented.");
    return EOSData();
  } 
  
  EOSData FromNAndT(const EOSData& eosIn);
  
  std::vector<double> GetExteriorDensities(const EOSData& eosIn);
  std::vector<double> GetTotalDensities(const EOSData& eosIn);
  
  double GetMinimumT() const {return mTMin;}
  double GetMaximumT() const {return 200.0/Constants::HBCFmMeV;}

  std::unique_ptr<EOSBase> MakeUniquePtr() const {
    return std::unique_ptr<EOSBase>(new EOSNSE(*this));
  }  

private: 
  std::vector<std::unique_ptr<NucleusBase>> mNuclei; 
  double mTMin;
  std::shared_ptr<EOSBase> mpEos; 

  std::array<double, 3> GetNucleiScalars(const EOSData& eosOut, double ne);

};

#endif // EOS_EOSSKYRME_HPP_
