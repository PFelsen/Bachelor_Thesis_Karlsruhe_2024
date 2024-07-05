add_mpp_test(basicalgebra/TestRVector LIB_MATH)
add_mpp_test(basicalgebra/TestCVector LIB_MATH)
add_mpp_test(basicalgebra/TestSymRMatrix LIB_MATH)
add_mpp_test(basicalgebra/TestAntisymRMatrix LIB_MATH)
add_mpp_test(basicalgebra/TestRMatrix LIB_MATH)
add_mpp_test(basicalgebra/TestHermCMatrix LIB_MATH)
add_mpp_test(basicalgebra/TestCMatrix LIB_MATH)
add_mpp_test(basicalgebra/TestCholesky LIB_MATH)
add_mpp_test(basicalgebra/TestInvert LIB_MATH)
add_mpp_test(basicalgebra/TestMatrixExp LIB_MATH)
add_mpp_test(basicalgebra/TestRTensor LIB_MATH)
add_mpp_test(basicalgebra/TestCTensor LIB_MATH)


add_mpp_test(basicsolver/TestBasicSolver LIB_MATH)


if (USE_FFTW)
    add_mpp_test(basicsolver/TestFFT LIB_MATH)
endif ()

add_mpp_test(spectrum/TestSpectrum LIB_MATH)

add_mpp_test(TestVectorField LIB_MATH)
add_mpp_test(TestVectorFieldComponent LIB_MATH)
add_mpp_test(TestSymTensor LIB_MATH)
add_mpp_test(TestTensor LIB_MATH)
add_mpp_test(TestTensorInvert LIB_MATH)
add_mpp_test(TestTensorRow LIB_MATH)
add_mpp_test(TestTensorComponent LIB_MATH)

add_mpp_bench(BenchVectorField LIB_MATH)


if (BUILD_IA_TESTS)
    add_mpp_test(basicalgebra/TestIARVector LIB_MATH)
    add_mpp_test(basicalgebra/TestIACVector LIB_MATH)
    add_mpp_test(basicalgebra/TestIASymRMatrix LIB_MATH)
    add_mpp_test(basicalgebra/TestIAAntisymRMatrix LIB_MATH)
    add_mpp_test(basicalgebra/TestIARMatrix LIB_MATH)
    add_mpp_test(basicalgebra/TestIAHermCMatrix LIB_MATH)
    add_mpp_test(basicalgebra/TestIACMatrix LIB_MATH)
    add_mpp_test(TestIAVectorField LIB_MATH)
    add_mpp_test(TestIASymTensor LIB_MATH)
    add_mpp_test(TestIATensor LIB_MATH)
endif ()