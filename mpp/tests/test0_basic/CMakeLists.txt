# Creates file for TestChainConfig
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/conf)
file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/conf/m++.conf "ConfigVerbose = 5")

add_mpp_test(TestPoint LIB_BASIC)
add_mpp_test(TestPointComparison LIB_BASIC)

add_mpp_test(utility/TestAssertion LIB_BASIC)
add_mpp_test(utility/TestConfig LIB_BASIC)
add_mpp_test(utility/TestChainConfig LIB_BASIC)
add_mpp_test(utility/TestSaveLoad LIB_BASIC)
add_mpp_test(utility/TestConfigFileCreation LIB_BASIC)
add_mpp_test(utility/TestLogging LIB_BASIC)

#add_mpp_test(parallel/TestBuffer LIB_BASIC)

if (BUILD_IA_TESTS)
    add_mpp_test(intervalarithmetic/TestIAInterval LIB_BASIC)
    add_mpp_test(TestIAPoint LIB_BASIC)
    add_mpp_test(utility/TestIASaveLoad LIB_BASIC)
endif ()

add_mpi_test(parallel/TestExchangeBuffer LIB_BASIC)
add_mpi_test(parallel/TestParallelization LIB_BASIC)
add_mpi_test(utility/TestLoggingParallel LIB_BASIC)
add_mpi_test(utility/TestMemoryLogger LIB_BASIC)