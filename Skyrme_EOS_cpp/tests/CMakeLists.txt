add_executable(skyrme_test skyrme_test.cpp) 
add_test(Skyrme skyrme_test)
target_link_libraries(skyrme_test eos ${EOSBUILDER_EXTERNAL_LIBS} )

add_executable(gibbs_test gibbs_test.cpp) 
add_test(Gibbs gibbs_test)
target_link_libraries(gibbs_test eos ${EOSBUILDER_EXTERNAL_LIBS} )

add_executable(table_test table_test.cpp) 
add_test(Table table_test)
target_link_libraries(table_test output ${EOSBUILDER_EXTERNAL_LIBS} )

add_executable(nucleus_test nucleus_test.cpp) 
#add_test(Nucleus nucleus_test)
target_link_libraries(nucleus_test eos ${EOSBUILDER_EXTERNAL_LIBS} )

add_executable(nse_test nse_test.cpp) 
#add_test(NSE nse_test)
target_link_libraries(nse_test eos ${EOSBUILDER_EXTERNAL_LIBS} )

add_executable(gibbs_archive gibbs_archive.cpp) 
target_link_libraries(gibbs_archive eos ${EOSBUILDER_EXTERNAL_LIBS} )

add_executable(gibbs_boundary gibbs_boundary.cpp) 
target_link_libraries(gibbs_boundary eos ${EOSBUILDER_EXTERNAL_LIBS} )

add_executable(pressure_comparison pressure_comparison.cpp) 
target_link_libraries(pressure_comparison eos ${EOSBUILDER_EXTERNAL_LIBS} )
