#file(COPY ${PROJECT_SOURCE_DIR}/src/lib7_application/spacetime/conf DESTINATION ${CMAKE_BINARY_DIR}/spacetime/)
#file(COPY ${PROJECT_SOURCE_DIR}/src/lib7_application/spacetime/models DESTINATION ${CMAKE_BINARY_DIR}/spacetime/)
file(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/spacetime/tests/data/vtu)
file(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/spacetime/tests/data/vtk)
file(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/spacetime/tests/data/gp)
file(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/spacetime/tests/data/py)
file(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/spacetime/tests/log)
file(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/spacetime/tests/data/dual)

# Mpi tests
if (SPACE_DIM STREQUAL "2")
    add_mpi_test(TestViscoAcousticMarmousi LIB_APP)
    add_mpi_test(TestNorms LIB_APP)
    add_mpi_test(TestViscoAcousticDegree LIB_APP)
    add_mpi_test(TestHConvergence LIB_APP)
    add_mpi_test(TestMarmousiWithNewImageClass LIB_APP)
    add_mpi_test(TestViscoAcousticDGTMarmousi LIB_APP)
    add_mpi_test(TestConformingInterpolation LIB_APP)
endif ()
#add_mpi_test(TestRiemannInterpolation LIB_APP)
#add_mpi_test(TestNorms LIB_APP)

#add_mpi_test(TestPolynomialSTProblems LIB_APP)
# Unit tests
add_mpp_test(TestPolynomialSTProblems LIB_APP)
if (SPACE_DIM STREQUAL "2")
    add_mpp_test(TestQuadratic LIB_APP)
    add_mpp_test(TestTransfer LIB_APP)
endif ()
#add_mpp_test(TestConformingProjection LIB_APP)
#add_mpp_test(TestMatrixFreeDGT LIB_APP)
#add_mpp_test(TestPreconditioner LIB_APP)
#add_mpp_test(TestMassMatrix LIB_APP)
#add_mpp_test(TestRestriction LIB_APP)

# Test not in pipeline
add_executable(TestRiemannInterpolation TestRiemannInterpolation.cpp)
target_link_libraries(TestRiemannInterpolation LIB_APP LIB_TEST)

