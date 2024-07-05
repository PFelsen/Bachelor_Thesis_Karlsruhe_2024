#include "SingleLevelEstimator.hpp"
#include "TestEnvironment.hpp"

struct TestParams {
  std::string problemName;

  std::string quantity;

  std::string model;

  double epsilon;

  int samples;

  int level = 3;

  bool onlyFine = false;
};

std::ostream &operator<<(std::ostream &s, const TestParams &testParams) {
  return s << "Problem Name: " << testParams.problemName << endl
           << "   Model: " << testParams.model << endl
           << "   Quantity: " << testParams.quantity << endl;
}

class TestMonteCarlo : public TestWithParam<TestParams> {
protected:
  int samples;

  double epsilon;

  PDESolverConfig pdeSolverConf;

  MultiSampleFEMConfig msFEMConf;

  SingleLevelEstimator mcSeriell;

  SingleLevelEstimator mcParallel;

  double Tolerance() const {
    if (epsilon == 0) return 2 * max(mcSeriell.TotalError(), mcParallel.TotalError());
    else return 2 * epsilon;
  }

  TestMonteCarlo() :
      samples(GetParam().samples), epsilon(GetParam().epsilon),
      pdeSolverConf(PDESolverConfig().WithDegree(1).WithModel(GetParam().model)),
      msFEMConf(MultiSampleFEMConfig(GetParam().problemName, GetParam().quantity)
                    .WithPDESolverConfig(pdeSolverConf)),
      mcSeriell(SingleLevelEstimator(SLEstimatorConfig()
                                         .WithMultiSampleFEMConfig(msFEMConf)
                                         .WithOnlyFine(GetParam().onlyFine)
                                         .WithInitLevel(GetParam().level)
                                         .WithInitSamples(samples)
                                         .WithEpsilon(epsilon)
                                         .WithParallel(false))),
      mcParallel(SingleLevelEstimator(SLEstimatorConfig()
                                          .WithMultiSampleFEMConfig(msFEMConf)
                                          .WithOnlyFine(GetParam().onlyFine)
                                          .WithInitLevel(GetParam().level)
                                          .WithInitSamples(samples)
                                          .WithEpsilon(epsilon)
                                          .WithParallel(true))) {

    mout << GetParam() << endl;

    mcSeriell.Method();
    mout << endl;
    mcSeriell.EstimatorResults();

    mcParallel.Method();
    mout << endl;
    mcParallel.EstimatorResults();
  }

  void TestTotalErrors() const {
    if (epsilon == 0.0) return;
    EXPECT_LE(mcSeriell.TotalError(), epsilon);
    EXPECT_LE(mcParallel.TotalError(), epsilon);
    EXPECT_LE(mcSeriell.NumericError(), epsilon);
    EXPECT_LE(mcParallel.NumericError(), epsilon);
    EXPECT_LE(mcSeriell.StochasticError(), epsilon);
    EXPECT_LE(mcParallel.StochasticError(), epsilon);
  }

  void CompareValues() const { EXPECT_NEAR(mcParallel.Value(), mcSeriell.Value(), Tolerance()); }

  void CompareData() const {
    EXPECT_NEAR(mcParallel.GetAggregate().Q.GetMean(), mcSeriell.GetAggregate().Q.GetMean(),
                Tolerance());
    EXPECT_NEAR(mcParallel.GetAggregate().Y.GetMean(), mcSeriell.GetAggregate().Y.GetMean(),
                Tolerance());

    EXPECT_NEAR(mcParallel.GetAggregate().Q.GetSVar(), mcSeriell.GetAggregate().Q.GetSVar(),
                Tolerance());
    EXPECT_NEAR(mcParallel.GetAggregate().Y.GetSVar(), mcSeriell.GetAggregate().Y.GetSVar(),
                Tolerance());
  }

  void TearDown() override { PPM->Barrier(0); }
};

INSTANTIATE_TEST_SUITE_P(
    TestMonteCarlo, TestMonteCarlo,
    Values(
        TestParams{"StochasticLaplace2D", "L2", "LagrangeElliptic", 0.0, 100},
        TestParams{"StochasticLaplace2D", "L2", "LagrangeElliptic", 0.05, 40},
        TestParams{"StochasticLaplace2DTest", "L2", "LagrangeElliptic", 0.0, 100},
        TestParams { "StochasticLaplace2DTest", "L2", "LagrangeElliptic", 0.05, 40 }
#ifdef USE_SPACETIME
//    ,
//    TestParams{"StochasticGaussHatAndRicker2D", "L2", "STDGViscoAcoustic", 0.05, 40} TODO:
//    "StochasticGaussHatAndRicker2D" is not a valid protocol!
#endif
        ));

TEST_P(TestMonteCarlo, TestAndCompareEstimators) {
  TestTotalErrors();
  CompareValues();
  CompareData();
}

int main(int argc, char **argv) {
  return MppTest(MppTestBuilder(argc, argv)
#if USE_SPACETIME
                     .WithConfigEntry("Preconditioner", "PointBlockGaussSeidel")
                     .WithConfigEntry("Distribution", "deformed_optimized")
                     .WithConfigEntry("Overlap", "STCellsWithCorners")
#endif
                     .WithConfigEntry("SLEstimatorVerbose", 2)
                     .WithConfigEntry("ParallelPlotting", 1)
                     .WithConfigEntry("PDESolverVerbose", 0)
                     .WithConfigEntry("GeneratorVerbose", 0)
                     .WithConfigEntry("AggregateVerbose", 2)
                     .WithConfigEntry("NewtonVerbose", 0)
                     .WithConfigEntry("LinearVerbose", 0)
                     .WithConfigEntry("ConfigVerbose", 0)
                     .WithConfigEntry("MeshVerbose", 0)
                     .WithConfigEntry("MainVerbose", 0)
                     .WithConfigEntry("MLMCVerbose", 0)
                     .WithConfigEntry("VtuPlot", 0)
                     .WithScreenLogging()
                     .WithRandomInitialized()
                     .WithPPM())
      .RUN_ALL_MPP_TESTS();
}