set(eos_files
  EOSData.cpp
  EOSSingleNucleus.cpp
  EOSSkyrme.cpp
  GibbsPhaseConstruct.cpp 
  NucleusBase.cpp
  EOSNSE.cpp
)
add_library(eos ${eos_files})
target_link_libraries(eos util)
