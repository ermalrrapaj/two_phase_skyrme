/// \file EOSTable.cpp
/// \author lroberts
/// \since Sep 03, 2015
///
/// \brief
///
///

#include "EOSTable.hpp"

EOSTable::EOSTable(const EOSBase& eos, double TMin, double TMax, double nbMin, 
    double nbMax, double yeMin, double yeMax, 
    std::size_t nT, std::size_t nNb, std::size_t nYe) : 
    mP(NDArray<double, 3>({nT, nNb, nYe})), 
    mS(NDArray<double, 3>({nT, nNb, nYe})), 
    mE(NDArray<double, 3>({nT, nNb, nYe})),
    mT(NDArray<double, 1>({nT})),
    mNb(NDArray<double, 1>({nNb})),
    mYe(NDArray<double, 1>({nYe})),
    mpEos(eos.MakeUniquePtr()) {
  
  double min = log10(TMin);
  double max = log10(TMax);
  for (int i=0; i<nT; ++i) 
    mT[i] = pow(10.0, min + (double) i / (nT-1.0) * (max - min));
  
  min = log10(nbMin);
  max = log10(nbMax);
  for (int i=0; i<nNb; ++i) 
    mNb[i] = pow(10.0, min + (double) i / (nNb-1.0) * (max - min));
  
  min = log10(yeMin);
  max = log10(yeMax);
  for (int i=0; i<nYe; ++i) 
    mYe[i] = pow(10.0, min + (double) i / (nYe-1.0) * (max - min));
  
  BuildTable();
}

void EOSTable::BuildTable() { 
  for (std::size_t i=0; i<mT.size(); ++i) {
    for (std::size_t k=0; k<mYe.size(); ++k) {
      for (std::size_t j=0; j<mNb.size(); ++j) {
        double T = mT[i];
        double nb = mNb[j];
        double ye = mYe[k]; 
        EOSData eosOut = mpEos->FromNAndT(EOSData::InputFromTNbYe(T, nb, ye));
        mP({i, j, k}) = eosOut.P();
        mS({i, j, k}) = eosOut.S();
        mE({i, j, k}) = eosOut.E();
      }
    }
  }
}
  
void EOSTable::WriteToH5(const H5::CommonFG& group) const {
  hsize_t asize[1] = {0}; 
  H5::DataSpace h5DSpace(1, asize); 
  
  mT.WriteToH5(group, "T")
      ->createAttribute("Units: fm^{-1}", H5::PredType::C_S1, h5DSpace);
  
  mNb.WriteToH5(group, "Nb")
      ->createAttribute("Units: fm^{-3}", H5::PredType::NATIVE_CHAR, h5DSpace);

  mYe.WriteToH5(group, "Ye")
      ->createAttribute("Units: None", H5::PredType::NATIVE_CHAR, h5DSpace);

  mP.WriteToH5(group, "P")
      ->createAttribute("Units: fm^{-4}", H5::PredType::NATIVE_CHAR, h5DSpace);

  mS.WriteToH5(group, "S")
      ->createAttribute("Units: kb / baryon", H5::PredType::NATIVE_CHAR, h5DSpace);

  mE.WriteToH5(group, "E")
      ->createAttribute("Units: fm^{-1} / baryon", H5::PredType::NATIVE_CHAR, h5DSpace);

}



