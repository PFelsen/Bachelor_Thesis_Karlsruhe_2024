#include "MultiLevelEstimator.hpp"
#include "TestEnvironment.hpp"


using TestParams = std::tuple<PDESolverConfig, MultiSampleFEMConfig,
    SLEstimatorConfig, MLEstimatorConfig>;

class TestMultilevelMonteCarlo : public TestWithParam<TestParams> {
protected:

  PDESolverConfig pdeSolverConfig;

  MultiSampleFEMConfig msFEMConfig;

  SLEstimatorConfig slEstmConfig;

  MLEstimatorConfig mlEstmConfig;

  MultiLevelEstimator mlmc;

  TestMultilevelMonteCarlo() :
      pdeSolverConfig(PDESolverConfig(std::get<0>(GetParam()))),
      msFEMConfig(std::get<1>(GetParam())),
      slEstmConfig(std::get<2>(GetParam())),
      mlEstmConfig(std::get<3>(GetParam())),

      mlmc(MultiLevelEstimator(mlEstmConfig.
          WithSLEstimatorConfig(slEstmConfig.
          WithMultiSampleFEMConfig(msFEMConfig.
          WithPDESolverConfig(pdeSolverConfig))))) {

    mout << mlEstmConfig << endl;
    mout << slEstmConfig << endl;

    mlmc.Method();
    mout << endl;

    mlmc.EstimatorResults();
    mlmc.MultilevelResults();
    mlmc.ExponentResults();
    mlmc.ContinuationResults();
  }

  void TestTotalErrors() const {
    if (mlEstmConfig.epsilon == 0.0) return;
    EXPECT_LE(mlmc.GetEstimatorMap().GetErrors().rmse, mlEstmConfig.epsilon);
    EXPECT_LE(mlmc.GetEstimatorMap().GetErrors().disc, mlEstmConfig.epsilon);
    EXPECT_LE(mlmc.GetEstimatorMap().GetErrors().input, mlEstmConfig.epsilon);
  }

  void TestTotalCost() const {
    if (mlEstmConfig.timeBudget == 0.0) return;
    // Time budget test is not sharp to prevent random pipeline failures
    EXPECT_LE(mlmc.Cost(), mlEstmConfig.timeBudget * 2);
  }

  void AssertExponents() const {
    EXPECT_TRUE(mlmc.GetEstimatorMap().GetExponents().alpha > 0);
    EXPECT_TRUE(mlmc.GetEstimatorMap().GetExponents().beta > 0);
    EXPECT_TRUE(mlmc.GetEstimatorMap().GetExponents().gamma > 0);
  }

  void TearDown() override { PPM->Barrier(0); }
};

INSTANTIATE_TEST_SUITE_P(TestMultilevelMonteCarlo, TestMultilevelMonteCarlo,
                         testing::Combine(
                             Values(
                                 PDESolverConfig().WithModel("LagrangeElliptic").WithDegree(1)
                             ),
                             Values(
//                                 MultiSampleFEMConfig("StochasticLaplace2DTest", "L2"),
                                 MultiSampleFEMConfig("StochasticLaplace2D", "L2")
                             ),
                             Values(
                                 SLEstimatorConfig()
                             ),
                             Values(
                                 MLEstimatorConfig().
                                     WithInitSamples({100, 100, 100}).
                                     WithInitLevel({3, 4, 5}).
                                     WithEpsilon(0.0),

                                 MLEstimatorConfig().
                                     WithInitSamples({16, 8, 4}).
                                     WithInitLevel({3, 4, 5}).
                                     WithEpsilon(0.01),

                                 MLEstimatorConfig().
                                     WithInitSamples({16, 8, 4}).
                                     WithInitLevel({3, 4, 5}).
                                     WithEpsilon(0.001),

                                 MLEstimatorConfig().
                                     WithInitSamples({16, 8, 4}).
                                     WithInitLevel({3, 4, 5}).
                                     WithTimeBudget(10),

                                 MLEstimatorConfig().
                                     WithInitSamples({16, 8, 4}).
                                     WithInitLevel({3, 4, 5}).
                                     WithTimeBudget(20),

                                 MLEstimatorConfig().
                                     WithInitSamples({16, 8, 4}).
                                     WithInitLevel({3, 4, 5}).
                                     WithTimeBudget(30)
                             )
                         )
);

TEST_P(TestMultilevelMonteCarlo, TestMultilevelMonteCarlo) {
  TestTotalErrors();
  TestTotalCost();
  AssertExponents();
}


int main(int argc, char **argv) {
  return MppTest(
      MppTestBuilder(argc, argv).
          WithConfigEntry("GeneratorVerbose", 0).
          WithConfigEntry("PDESolverVerbose", 0).
          WithConfigEntry("AggregateVerbose", 0).
          WithConfigEntry("NewtonVerbose", 0).
          WithConfigEntry("LinearVerbose", 0).
          WithConfigEntry("ConfigVerbose", 0).
          WithConfigEntry("MeshVerbose", 0).
          WithConfigEntry("MainVerbose", 0).
          WithConfigEntry("MLMCVerbose", 1).
          WithConfigEntry("SLEstimatorVerbose", 1).
          WithRandomInitialized().
          WithScreenLogging().
          WithPPM()
  ).RUN_ALL_MPP_TESTS();
}
