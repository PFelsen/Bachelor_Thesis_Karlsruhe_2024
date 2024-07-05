add_mpp_test(quadrature/TestQuadrature LIB_DISC)
add_mpp_test(quadrature/TestGaussLobattoQuadrature LIB_DISC)
add_mpp_test(quadrature/TestTensorQuadrature LIB_DISC)

add_mpp_test(dof/TestLagrangeDoF LIB_DISC)
add_mpp_test(dof/TestSerendipityDoF LIB_DISC)
add_mpp_test(dof/TestRTDoF LIB_DISC)
add_mpp_test(dof/TestArgyrisDoF LIB_DISC)

add_mpp_test(shapes/TestLagrangeShape LIB_DISC)
add_mpp_test(shapes/TestGaussLobattoShapes LIB_DISC)
add_mpp_test(shapes/TestArgyrisShape LIB_DISC)
add_mpp_test(shapes/TestRTShape LIB_DISC)

add_mpp_test(matrixgraph/TestRows LIB_DISC)
add_mpi_test(matrixgraph/TestMatrixGraph LIB_DISC)

if (USE_SPACETIME)
    add_mpp_test(shapes/TestSpaceTimeShape LIB_DISC)
endif ()

if (BUILD_IA_TESTS)
    add_mpp_test(quadrature/TestIAQuadrature LIB_DISC)
endif ()