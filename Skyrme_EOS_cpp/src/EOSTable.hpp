/// \file EOSTable.hpp
/// \author lroberts
/// \since Sep 03, 2015
///
/// \brief
///
///

#ifndef EOS_EOSTABLE_HPP_
#define EOS_EOSTABLE_HPP_

#include <memory> 
#include <vector> 
#include <utility> 
#include <H5Cpp.h> 

#include "EOSBase.hpp"

template <class T>
class ThreeDArray {
public:
  ThreeDArray() { mData = new T[1]; } 
  ThreeDArray(std::size_t n1, std::size_t n2, std::size_t n3) :
      mN1(n1), mN2(n2), mN3(n3) {
    mData = new T[mN1*mN2*mN3];
  } 
  
  ~ThreeDArray() {delete [] mData;}  
  
  T& operator()(std::size_t idx1, std::size_t idx2, std::size_t idx3) {
    std::size_t idx = idx1 + mN1*(idx2 + mN2*idx3);
    if (idx >= mN1*mN2*mN3) 
      throw std::length_error("Trying to access element in matrix beyond size.");
    return mData[idx];
  }     
  
  T operator()(std::size_t idx1, std::size_t idx2, std::size_t idx3) const {
    std::size_t idx = idx1 + mN1*(idx2 + mN2*idx3);
    if (idx >= mN1*mN2*mN3) 
      throw std::length_error("Trying to access element in matrix beyond size.");
    return mData[idx];
  } 

protected:
  T *mData;
  std::size_t mN1, mN2, mN3;

};

///
/// Create an EOSTable from a specified EOS object
///
class EOSTable {
public:
  EOSTable(const EOSBase& eos) : mpEos(eos.MakeUniquePtr()) {} 
  EOSTable(const EOSBase& eos, double TMin, double TMax, double nbMin, 
      double nbMax, double yeMin, double yeMax, int nT, int nNb, int nYe); 
   
protected:
  void BuildTable(); 
  std::unique_ptr<EOSBase> mpEos;
  std::vector<double> mT, mNb, mYe; 
  ThreeDArray<double> mP, mS, mE;
};

#endif // EOS_EOSTABLE_HPP_
