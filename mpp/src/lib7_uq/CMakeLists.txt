set(uq_src
        generators/Random.cpp
        generators/CirculantEmbedding.cpp
        generators/SymmetricCovariance.cpp
        generators/KLExpansionGenerator.cpp
        generators/HybridFluxGenerator.cpp

        problems/StochasticReactionProblem.cpp
        problems/StochasticTransportProblems.cpp
        problems/StochasticAcousticProblems.cpp

        estimators/SingleLevelEstimator.cpp
        estimators/MultiLevelEstimator.cpp
        estimators/EstimatorMap.cpp
        estimators/WelfordAggregate.cpp
        estimators/MultiSampleFEM.cpp

        sgd/StochasticGradientDescent.cpp

        smc/SequentialMonteCarlo.cpp
        smc/model/InverseProblem.cpp
        smc/model/PDEModel.cpp
        smc/model/ObservationOperator.cpp
        smc/generators/ProposalGenerators.cpp
        smc/generators/SamplerSelection.cpp
        smc/generators/KLExpansion.cpp
        smc/tools/Metropolis.cpp
        smc/tools/ParticleSet.cpp
        smc/tools/Particle.cpp
        smc/tools/Resampling.cpp
        smc/tools/QuantityOfInterest.cpp
        smc/tools/OtherTools.cpp
        smc/statistical_tests/ChiSquared.cpp
        smc/statistical_tests/KolmogorovSmirnov.cpp
        smc/statistical_tests/ShapiroWilk.cpp
        smc/statistical_tests/StatisticalTest.cpp
        smc/statistical_tests/TargetDistributions.cpp
        smc/statistical_tests/Binning.cpp
)

if (USE_TASMANIAN)
    set(uq_src ${uq_src}
            generators/SparseGridGenerator.cpp
    )
endif ()

if (USE_SPACETIME)
    set(uq_src ${uq_src}
            problems/StochasticSTAcousticProblems.cpp
            problems/StochasticSTTransportProblems.cpp
    )
endif ()

add_library(LIB_UQ STATIC ${uq_src})
target_link_libraries(LIB_UQ LIB_APP sprng)

if (USE_TASMANIAN)
    target_link_libraries(LIB_UQ Tasmanian_libsparsegrid)
endif ()