/// \file NDArray.hpp
/// \author lroberts
/// \since Sep 09, 2015
///
/// \brief
///
///

#ifndef EOS_NDARRAY_HPP_
#define EOS_NDARRAY_HPP_

#include <memory> 
#include <vector> 
#include <utility> 
#include <H5Cpp.h> 

template<class T> 
H5::DataType GetH5DataType() {
  if (std::is_same<T, int>::value) {
    return H5::DataType(H5::PredType::NATIVE_INT);
  }
  else if (std::is_same<T, float>::value) {
    return H5::DataType(H5::PredType::NATIVE_FLOAT);
  }
  else if (std::is_same<T, double>::value) {
    return H5::DataType(H5::PredType::NATIVE_DOUBLE);
  }
  else if (std::is_same<T, long>::value) {
    return H5::DataType(H5::PredType::NATIVE_LONG);
  }
  else if (std::is_same<T, char>::value) {
    return H5::DataType(H5::PredType::NATIVE_CHAR);
  } else {
    throw std::invalid_argument("Unknown type");
  }
}


template <class T, int Nd>
class NDArray {
public:
  NDArray() { mData = new T[1]; } 
  NDArray(std::array<std::size_t, Nd> n) {
    mN = n;
    mNTot = std::accumulate(mN.begin(), mN.end(), 1, 
        std::multiplies<std::size_t>());
    mData = new T[mNTot];
  } 
  
  ~NDArray() {delete [] mData;}  
  
  T& operator()(std::array<std::size_t, Nd> idxs) {
    return mData[GetIdx(idxs)];
  }     
  
  T operator()(std::array<std::size_t, Nd> idxs) const { 
    return mData[GetIdx(idxs)];
  } 
  
  T& operator[](std::size_t idx) {
    if (idx >= mNTot) 
      throw std::length_error("Index out of allocated memory region.");
    return mData[idx];
  }     
  
  T operator[](std::size_t idx) const {
    if (idx >= mNTot) 
      throw std::length_error("Index out of allocated memory region.");
    return mData[idx];
  }    
   
  std::size_t size() const { return mNTot; }
    
  void WriteToH5(const H5::CommonFG& group, const std::string& dsetName) const {
    H5::DataType h5DType = GetH5DataType<T>(); 
    H5::DataSpace h5DSpace(Nd, mN.data()); 
    auto dset = group.createDataSet(dsetName.c_str(), h5DType, h5DSpace);
    dset.write(mData, h5DType); 
  }

protected:
  std::size_t GetIdx(const std::array<std::size_t, Nd>& idxs) const {
  std::size_t idx = 0;
  for (int i=Nd-1; i>0; --i) {
    if (idxs[i] >= mN[i]) 
      throw std::length_error("Trying to access element in matrix beyond size.");
    idx += idxs[i];
    idx *= mN[i-1];
  }
  idx += idxs[0];
   
  if (idx >= mNTot) 
    throw std::length_error("Index out of allocated memory region.");
}

  T *mData;

  std::size_t mNTot; 
  std::array<std::size_t, Nd> mN;

};

#endif // EOS_NDARRAY_HPP_
