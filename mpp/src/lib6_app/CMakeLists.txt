set(app_src
        problems/ProblemBase.cpp
        problems/TransportProblems.cpp
        problems/ParameterImage.cpp
        problems/lodepng.cpp
        problems/AcousticProblems.cpp
        problems/EllipticProblems.cpp
        problems/ReactionProblems.cpp
        problems/GelfandProblem.cpp
        problems/PlateProblem.cpp
        problems/VectorValuedProblems.cpp

        assemble/IEllipticAssemble.cpp
        assemble/LagrangeEllipticAssemble.cpp
        assemble/HybridEllipticAssemble.cpp
        assemble/MixedEllipticAssemble.cpp
        assemble/DGEllipticAssemble.cpp
        assemble/EGEllipticAssemble.cpp
        assemble/ITransportAssemble.cpp
        assemble/DGTransportAssemble.cpp
        assemble/IReactionAssemble.cpp
        assemble/PGReactionAssemble.cpp
        assemble/DGReactionAssemble.cpp
        assemble/IAcousticAssemble.cpp
        assemble/DGAcousticAssemble.cpp
        assemble/IVectorValuedAssemble.cpp
        assemble/LagrangeVectorValuedAssemble.cpp
        assemble/DGVectorValuedAssemble.cpp
        assemble/EGVectorValuedAssemble.cpp
        assemble/GelfandAssemble.cpp
        assemble/PlateAssemble.cpp

        method/SingleExecution.cpp
        method/ConvergenceStudy.cpp
        method/AdaptiveConfig.cpp
        method/AdaptiveStudy.cpp
        method/AdaptiveConvergenceMethod.cpp

        pdesolvers/EllipticPDESolver.cpp
        pdesolvers/TransportPDESolver.cpp
        pdesolvers/ReactionPDESolver.cpp
        pdesolvers/AcousticPDESolver.cpp
        pdesolvers/GelfandPDESolver.cpp
        pdesolvers/PlatePDESolver.cpp
        pdesolvers/MixedMain.cpp
        pdesolvers/VectorValuedMain.cpp

        seismogram/SeismogramData.cpp
)

if (USE_SPACETIME)
    # TODO: This should not be set here! Rather use suboption!
    set_option(SPACE_DIM 2)
    set_option(AFFINE_LINEAR_TRAFO ON)
    set_option(BUILD_SPACETIME_TESTS ON)

    # TODO: Remove spacetime directory and merge code with other directories
    set(app_src
            ${app_src}
            spacetime/problems/ViscoElasticProblems.cpp

            spacetime/assemble/STPGViscoAcousticAssemble.cpp
            spacetime/assemble/STDGMatrixFreeViscoAcousticAssemble.cpp
            spacetime/assemble/STDGViscoAcousticAssemble.cpp
            spacetime/assemble/STDGTransportAssemble.cpp
            spacetime/assemble/DGViscoAcousticDGDG_GLGL.cpp
            spacetime/assemble/STDGElasticity.cpp
            spacetime/assemble/STPGViscoElasticAssemble.cpp
            spacetime/assemble/STGLViscoAcousticAssemble.cpp
            spacetime/assemble/STAssemble.cpp

            spacetime/main/STAcousticMain.cpp
            spacetime/main/STTransportMain.cpp
            spacetime/main/STAcousticPDESolver.cpp
            spacetime/main/STTransportPDESolver.cpp

            spacetime/elements/STDGDGViscoAcousticElement.cpp
            spacetime/elements/STDGDGViscoAcousticFaceElement.cpp
            spacetime/elements/STDGDGTransportElement.cpp
            spacetime/elements/STDGDGTransportFaceElement.cpp
            spacetime/elements/SpaceTimeElasticityElement.cpp
            spacetime/elements/SpaceTimeViscoElasticElement.cpp
            spacetime/elements/SpaceTimeViscoAcousticFaceElement.cpp
            spacetime/elements/ViscoAcousticElements.cpp
            spacetime/elements/STDGDGViscoAcousticElementGLGL.cpp
            spacetime/elements/STDGDGViscoAcousticFaceElementGLGL.cpp

            spacetime/CImg.hpp
            spacetime/SpaceTimePreconditioner.cpp
            spacetime/SpaceTimeTools.cpp
            spacetime/STPathManager.cpp
            spacetime/SpaceTimeMeasurement.cpp
            spacetime/Results.cpp
            spacetime/PrintUtil.cpp
            spacetime/DebuggingTools.cpp
            spacetime/SpaceTimeTransfer.cpp
            spacetime/ErrorEstimator.cpp
            spacetime/STCallbackStrategy.cpp

            spacetime/method/FWISTExperiment.cpp
    )
endif ()

add_library(LIB_APP STATIC ${app_src})
target_link_libraries(LIB_APP LIB_SOLVE)

if (USE_CXSC)
    add_subdirectory(spectrum)
endif ()
