set(fem_src
        algebra/VectorMatrixBase.cpp
        algebra/Vector.cpp
        algebra/VectorAccess.cpp
        discretization/ArgyrisDiscretization.cpp
        discretization/DGDiscretization.cpp
        discretization/MixedEGDiscretization.cpp
        discretization/DivergenceFreeDiscretization.cpp
        discretization/IDiscretization.cpp
        discretization/LagrangeDisplacement.cpp
        discretization/RTDiscretization.cpp
        discretization/RTLagrangeDiscretization.cpp
        elements/ScalarElement.cpp
        elements/VectorFieldElement.cpp
        elements/TensorElement.cpp
        elements/MixedScalarElement.cpp
        elements/MixedVectorFieldElement.cpp
        elements/TaylorHoodElement.cpp
        elements/RTElement.cpp
        elements/RTLagrangeElement.cpp
        elements/ArgyrisElement.cpp
        elements/DivergenceFreeElement.cpp
        elements/DGElement.cpp
        elements/DGVectorFieldElement.cpp
        elements/DGAcousticElement.cpp
        elements/MixedEGElement.cpp
        elements/MixedEGVectorFieldElement.cpp
        elements/CurlElement.cpp
        elements/MultiPartScalarElement.cpp
        elements/MultiPartVectorFieldElement.cpp
        elements/MultiPartScalarWGElement.cpp
        algebra/Operator.cpp
        algebra/Matrix.cpp
        algebra/MatrixAccess.cpp
        algebra/Sparse.cpp
        algebra/VectorProcSet.cpp
        transfer/Rowlist.cpp
        transfer/LinearTransfer.cpp
        transfer/LagrangeTransfer.cpp
        transfer/MatrixTransfer.cpp
        transfer/MixedLagrangeTransfer.cpp
        transfer/Transfers.cpp
        transfer/TransferWrapper.cpp
        plot/Plotting.cpp
        plot/VtuPlot.cpp
)

include_subdirectories(elements/deprecated_element_header)

if(USE_SUPERLU)
    set(fem_src ${fem_src} algebra/ParallelMatrix.cpp)
endif()

if(USE_SPACETIME)
    set (fem_src ${fem_src}
            discretization/SpaceTimeDiscretization.cpp
            discretization/SpaceTimeDiscretization_PGDG.cpp
            discretization/SpaceTimeDiscretization_DGDG.cpp
            discretization/SpaceTimeDiscretization_DGDG_GLGL.cpp
            )
endif()

add_library(LIB_FEM SHARED ${fem_src})

set(fem_deps LIB_DISC)

if (ZLIB_FOUND)
    set(fem_deps ${fem_deps} ZLIB::ZLIB)
endif ()

target_link_libraries(LIB_FEM ${fem_deps})
