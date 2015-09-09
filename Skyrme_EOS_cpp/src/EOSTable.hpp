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
#include "NDArray.hpp"

///
/// Create an EOSTable from a specified EOS object
///
class EOSTable {
public:
  EOSTable(const EOSBase& eos) : mpEos(eos.MakeUniquePtr()) {} 
  EOSTable(const EOSBase& eos, double TMin, double TMax, double nbMin, 
      double nbMax, double yeMin, double yeMax, 
      std::size_t nT, std::size_t nNb, std::size_t nYe); 
   
protected:
  void BuildTable(); 
  std::unique_ptr<EOSBase> mpEos;
  NDArray<double, 1> mT, mNb, mYe; 
  NDArray<double, 3> mP, mS, mE;

};

#endif // EOS_EOSTABLE_HPP_
