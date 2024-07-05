set (math_src
        basicalgebra/RVector.cpp
        basicalgebra/SymRMatrix.cpp
        basicalgebra/AntisymRMatrix.cpp
        basicalgebra/RMatrix.cpp
        basicalgebra/CVector.cpp
        basicalgebra/HermCMatrix.cpp
        basicalgebra/CMatrix.cpp
        basicalgebra/SparseRMatrix.cpp
        basicalgebra/RTensor.cpp
        basicalgebra/CTensor.cpp
        basicsolver/BasicSolver.cpp
        basicsolver/BasicNonLinearSolver.cpp
        basicsolver/QRDecomposition.cpp
        spectrum/Spectrum.cpp
        DataSet.cpp
        VectorFieldComponent.cpp
        Transformation.cpp
        )

if (USE_SUPERLU)
    set (math_src ${math_src} basicsolver/BasicSparseSolver.cpp)
endif ()

if (USE_CXSC)
    set(math_src ${math_src}
            spectrum/IAEigenvalueBounds.cpp
            spectrum/IASpectrum.cpp
            )
endif ()

if(USE_FFTW)
    set(math_src ${math_src}
            basicsolver/FFT.cpp
            )
endif ()

add_library(LIB_MATH SHARED ${math_src})
if(USE_FFTW)
    if(USE_SUPERLU)
         target_link_libraries(LIB_MATH LIB_BASIC fftw3 superlu)
    else()
         target_link_libraries(LIB_MATH LIB_BASIC fftw3)
    endif()
else()
    if(USE_SUPERLU)
         target_link_libraries(LIB_MATH LIB_BASIC superlu)
    else()
         target_link_libraries(LIB_MATH LIB_BASIC)
    endif()
endif()
