add_executable(eosdata_test eosdata_test.cpp) 
add_test(EOSData eosdata_test)
target_link_libraries(eosdata_test eos ${EOSBUILDER_EXTERNAL_LIBS} )

add_executable(skyrme_test skyrme_test.cpp) 
add_test(Skyrme skyrme_test)
target_link_libraries(skyrme_test eos ${EOSBUILDER_EXTERNAL_LIBS} )

add_executable(eoselectron_test eoselectron_test.cpp) 
add_test(EOSElectron eoselectron_test)
target_link_libraries(eoselectron_test eos ${EOSBUILDER_EXTERNAL_LIBS} )

add_executable(gibbs_test gibbs_test.cpp) 
add_test(Gibbs gibbs_test)
target_link_libraries(gibbs_test eos ${EOSBUILDER_EXTERNAL_LIBS} )

add_executable(table_test table_test.cpp) 
add_test(Table table_test)
target_link_libraries(table_test output ${EOSBUILDER_EXTERNAL_LIBS} )

add_executable(fermi_test test_fermi_deriv.f90) 
add_test(Fermi fermi_test)
target_link_libraries(fermi_test eos ${EOSBUILDER_EXTERNAL_LIBS} )

add_executable(nucleus_test nucleus_test.cpp) 
#add_test(Nucleus nucleus_test)
target_link_libraries(nucleus_test eos ${EOSBUILDER_EXTERNAL_LIBS} )

add_executable(nse_test nse_test.cpp) 
#add_test(NSE nse_test)
target_link_libraries(nse_test eos ${EOSBUILDER_EXTERNAL_LIBS} )

add_executable(nse_test_new nse_test_new.cpp) 
#add_test(NSE nse_test)
target_link_libraries(nse_test_new eos ${EOSBUILDER_EXTERNAL_LIBS} )

add_executable(nse_guess_test nse_guess_test.cpp) 
#add_test(NSE nse_guess_test)
target_link_libraries(nse_guess_test eos ${EOSBUILDER_EXTERNAL_LIBS} )

add_executable(gibbs_archive gibbs_archive.cpp) 
target_link_libraries(gibbs_archive eos ${EOSBUILDER_EXTERNAL_LIBS} )

add_executable(gibbs_boundary gibbs_boundary.cpp) 
target_link_libraries(gibbs_boundary eos ${EOSBUILDER_EXTERNAL_LIBS} )

add_executable(pressure_comparison pressure_comparison.cpp) 
target_link_libraries(pressure_comparison eos ${EOSBUILDER_EXTERNAL_LIBS} )
