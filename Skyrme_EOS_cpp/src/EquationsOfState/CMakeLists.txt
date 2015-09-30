set(eos_files
  fermidirac_integrals.f90
  EOSData.cpp
  EOSSingleNucleus.cpp
  EOSSkyrme.cpp
  GibbsPhaseConstruct.cpp 
  NucleusBase.cpp
  EOSNSE.cpp
  EOSElectron.cpp
)
add_library(eos ${eos_files})
target_link_libraries(eos util)
