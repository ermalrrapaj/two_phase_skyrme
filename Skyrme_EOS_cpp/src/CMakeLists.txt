set(eos_files
  EOSData.cpp
  EOSSingleNucleus.cpp
  EOSSkyrme.cpp
  EOSTable.cpp
  GibbsPhaseConstruct.cpp 
  funct_fermi1.f90
)
add_library( eos ${eos_files} )

add_executable(skyrme_exe skyrme_exe.cpp) 
target_link_libraries(skyrme_exe eos ${EOSBUILDER_EXTERNAL_LIBS} )

add_executable(gibbs_archive gibbs_archive.cpp) 
target_link_libraries(gibbs_archive eos ${EOSBUILDER_EXTERNAL_LIBS} )

add_executable(gibbs_boundary gibbs_boundary.cpp) 
target_link_libraries(gibbs_boundary eos ${EOSBUILDER_EXTERNAL_LIBS} )

add_executable(gibbs_test gibbs_test.cpp) 
target_link_libraries(gibbs_test eos ${EOSBUILDER_EXTERNAL_LIBS} )
