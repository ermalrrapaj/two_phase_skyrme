/// \file EOSTable.cpp
/// \author lroberts
/// \since Sep 03, 2015
///
/// \brief
///
///

#include "EOSTable.hpp"

//using namespace {
//
//template<class T> 
//H5::DataType GetH5DataType() {
//  if (std::is_same<T, int>::value) {
//    return H5::DataType(H5::PredType::NATIVE_INT);
//  }
//  else if (std::is_same<T, float>::value) {
//    return H5::DataType(H5::PredType::NATIVE_FLOAT);
//  }
//  else if (std::is_same<T, double>::value) {
//    return H5::DataType(H5::PredType::NATIVE_DOUBLE);
//  }
//  else if (std::is_same<T, long>::value) {
//    return H5::DataType(H5::PredType::NATIVE_LONG);
//  }
//  else if (std::is_same<T, char>::value) {
//    return H5::DataType(H5::PredType::NATIVE_CHAR);
//  } else {
//    throw std::invalid_argument("Unknown type");
//  }
//}
//
//}

EOSTable::EOSTable(const EOSBase& eos, double TMin, double TMax, double nbMin, 
    double nbMax, double yeMin, double yeMax, int nT, int nNb, int nYe) : 
    mP(ThreeDArray<double>(nT, nNb, nYe)), 
    mS(ThreeDArray<double>(nT, nNb, nYe)), 
    mE(ThreeDArray<double>(nT, nNb, nYe)),
    mpEos(eos.MakeUniquePtr()) {
  
  mT = std::vector<double>(nT);
  double min = log10(TMin);
  double max = log10(TMax);
  for (int i=0; i<nT; ++i) 
    mT[i] = pow(10.0, min + (double) i / (nT-1.0) * (max - min));
  
  mNb = std::vector<double>(nNb); 
  min = log10(nbMin);
  max = log10(nbMax);
  for (int i=0; i<nNb; ++i) 
    mNb[i] = pow(10.0, min + (double) i / (nNb-1.0) * (max - min));
  
  mYe = std::vector<double>(nYe); 
  min = log10(yeMin);
  max = log10(yeMax);
  for (int i=0; i<nYe; ++i) 
    mYe[i] = pow(10.0, min + (double) i / (nYe-1.0) * (max - min));
  
  BuildTable();
}

void EOSTable::BuildTable() { 
  for (int i=0; i<mT.size(); ++i) {
    for (int k=0; k<mYe.size(); ++k) {
      for (int j=0; j<mNb.size(); ++j) {
        double T = mT[i];
        double nb = mNb[j];
        double ye = mYe[k]; 
        EOSData eosOut = mpEos->FromNAndT(EOSData::InputFromTNbYe(T, nb, ye));
        mP(i, j, k) = eosOut.P();
        mS(i, j, k) = eosOut.S();
        mE(i, j, k) = eosOut.E();
      }
    }
  }
}
