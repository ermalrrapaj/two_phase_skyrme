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
#include <string>
#include <algorithm> 
#include <H5Cpp.h> 
#include <iostream> 

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
  }
  else if (std::is_same<T, std::string>::value) {
    return H5::StrType(H5::PredType::C_S1);
  } else {
    throw std::invalid_argument("Unknown type");
  }
}


template <class T, int Nd>
class NDArray {
public:
  NDArray(): mNTot(1) { mData = new T[mNTot];}
   
  NDArray(std::array<std::size_t, Nd> n) {
    mN = n;
    mNTot = std::accumulate(mN.begin(), mN.end(), 1, 
        std::multiplies<std::size_t>());
    mData = new T[mNTot];
  } 
  
  NDArray(const NDArray<T,Nd>& old) {
    mNTot = old.mNTot;
    mN = old.mN;
    mData = new T[mNTot]; 
    for (unsigned int i=0; i<mNTot; ++i) mData[i] = old.mData[i];
  } 
  
  NDArray<T,Nd>& operator=(const NDArray<T,Nd>& old) {
    if (this == &old) return *this; 
    delete [] mData; 
    mNTot = old.mNTot;
    mN = old.mN;
    mData = new T[mNTot]; 
    for (unsigned int i=0; i<mNTot; ++i) mData[i] = old.mData[i];
  }
  
  ~NDArray() {
    delete [] mData;
   }  
  
  static NDArray<T, Nd> ReadFromH5(const H5::DataSet& h5Dset);
  
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
  
  std::size_t dimSize(int dim) const { return mN[dim];} 
    
  std::shared_ptr<H5::DataSet> 
  WriteToH5(const H5::CommonFG& group, const std::string& dsetName) const;

protected:
  T* mData;
  std::size_t mNTot; 
  std::array<std::size_t, Nd> mN;
    
  std::size_t GetIdx(const std::array<std::size_t, Nd>& idxs) const;

};

// Implementations have to be in header for template classes
template <class T, int Nd>
NDArray<T, Nd> NDArray<T,Nd>::ReadFromH5(const H5::DataSet& h5Dset) {
  H5::DataSpace dspace = h5Dset.getSpace(); 
  int ndim = dspace.getSimpleExtentNdims(); 
  if (ndim>Nd) 
      throw std::range_error("Too many dimensions in H5 dataset for NDArray");
  hsize_t dimSize[ndim];
  dspace.getSimpleExtentDims(dimSize); 
  std::array<std::size_t, Nd> dimSizeArr;
  for (int i=0; i<Nd; ++i) dimSizeArr[i] = dimSize[i];
  NDArray<T, Nd> arr(dimSizeArr);
  // Read in data here
  H5::DataType h5DType = GetH5DataType<T>();
  h5Dset.read(arr.mData, h5DType);
  return arr;
}

template <class T, int Nd>
std::size_t NDArray<T, Nd>::GetIdx(const std::array<std::size_t, Nd>& idxs) const {
  std::size_t idx = 0;
  
  // Storage needs to be in this direction to work with HDF5
  for (int i=0; i<Nd-1; ++i) {
    if (idxs[i] >= mN[i]) 
      throw std::length_error("Trying to access element in matrix beyond size.");
    idx += idxs[i];
    idx *= mN[i+1];
  } 
  idx += idxs[Nd-1];
  if (idx >= mNTot) 
    throw std::length_error("Index out of allocated memory region.");
  
  return idx;
}

template <class T, int Nd>
std::shared_ptr<H5::DataSet> NDArray<T, Nd>::WriteToH5(
    const H5::CommonFG& group, const std::string& dsetName) const {
  H5::DataType h5DType = GetH5DataType<T>();
  hsize_t dims[Nd]; 
  for (int i=0; i<Nd; ++i) dims[i] = mN[i]; 
  H5::DataSpace h5DSpace(Nd, dims); 
  std::shared_ptr<H5::DataSet> dset = std::make_shared<H5::DataSet>( 
      group.createDataSet(dsetName.c_str(), h5DType, h5DSpace));
  dset->write(mData, h5DType); 
  return dset;
}

#endif // EOS_NDARRAY_HPP_
