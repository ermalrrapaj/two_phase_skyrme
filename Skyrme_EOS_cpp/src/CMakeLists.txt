add_subdirectory(EquationsOfState)
add_subdirectory(Util)
add_subdirectory(Output)

add_executable(skyrme_exe skyrme_exe.cpp) 
target_link_libraries(skyrme_exe eos ${EOSBUILDER_EXTERNAL_LIBS} )

add_executable(gibbs_archive gibbs_archive.cpp) 
target_link_libraries(gibbs_archive eos ${EOSBUILDER_EXTERNAL_LIBS} )

add_executable(gibbs_boundary gibbs_boundary.cpp) 
target_link_libraries(gibbs_boundary eos ${EOSBUILDER_EXTERNAL_LIBS} )

add_executable(gibbs_test gibbs_test.cpp) 
target_link_libraries(gibbs_test eos ${EOSBUILDER_EXTERNAL_LIBS} )

add_executable(table_test table_test.cpp) 
target_link_libraries(table_test output ${EOSBUILDER_EXTERNAL_LIBS} )
