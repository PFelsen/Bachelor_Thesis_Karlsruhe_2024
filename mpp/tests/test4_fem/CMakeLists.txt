file(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/mpp/tests/test4_fem/data/vtu)
file(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/tests/test4_fem/data/vtu)

add_mpp_test(discretization/TestLagrangeDiscretization LIB_FEM)
add_mpp_test(discretization/TestSerendipityDiscretization LIB_FEM)
add_mpp_test(discretization/TestRTDiscretization LIB_FEM)
add_mpp_test(discretization/TestRTLagrangeDiscretization LIB_FEM)
add_mpp_test(discretization/TestArgyrisDiscretization LIB_FEM)

add_mpp_test(elements/TestLagrangeElement LIB_FEM)
add_mpp_test(elements/TestArgyrisElement LIB_FEM)
add_mpp_test(elements/TestDivergenceFreeElement LIB_FEM)
add_mpp_test(elements/TestRTElement LIB_FEM)
add_mpp_test(elements/TestVectorvaluedElements LIB_FEM)

add_mpp_test(TestTransfers LIB_FEM)

add_mpi_test(algebra/TestVector LIB_FEM)
add_mpi_test(TestTransfers LIB_FEM)

if (BUILD_IA_TESTS)
    add_mpp_test(elements/TestBernsteinRanges LIB_FEM)
endif ()

# Todo clean this up below - linking and libs are all over the place
add_executable(TestVtuPlotMPI plot/TestVtuPlot.cpp ${PROJECT_MPP_DIR}/submodules/pugixml/src/pugixml.cpp)
target_include_directories(TestVtuPlotMPI PRIVATE ${PROJECT_MPP_DIR}/submodules/pugixml/src)
if (USE_SPACETIME)
    target_link_libraries(TestVtuPlotMPI LIB_FEM LIB_TEST LIB_APP)
else ()
    target_link_libraries(TestVtuPlotMPI LIB_FEM LIB_TEST)
endif ()

file(APPEND ${PROJECT_BINARY_DIR}/mpi_tests.txt
        "${CMAKE_CURRENT_BINARY_DIR}/TestVtuPlotMPI\n")