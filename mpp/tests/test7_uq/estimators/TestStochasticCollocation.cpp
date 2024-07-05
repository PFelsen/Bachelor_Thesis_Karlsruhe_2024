#include "TestStochasticCollocation.hpp"


// Test cases for convergence rate
// depending on dimension and level -> Smoothness of problem has to be known

// Todo: fix

INSTANTIATE_TEST_SUITE_P(TestStochasticCollocation, TestStochasticCollocationWithoutEpsilon,
                         Values(
                             // Dummy PDESolver
                             //  TestParams{"SparseGrid2DGeneratorProblem", "FunctionEvaluation",
                             //             "DummyPDESolver", 2.513723354063905},
                             //  // Todo
                             //  TestParams{"SparseGrid2DHermiteProblem", "FunctionEvaluation",
                             //             "DummyPDESolver", 1.0},

                             // EllipticPDESolver
                             TestParams{"SparseGridLaplace2D", "L2", "LagrangeElliptic",
                                        0.5827913023304272, 3, false}

                             ));

TEST_P(TestStochasticCollocationWithoutEpsilon, TestSeriellAgainstParallel) {
  mout << GetParam() << endl;

  scSeriell->Method();
  mout << endl;
  scSeriell->EstimatorResults();

  EXPECT_NEAR(scSeriell->aggregate.mean.Q, GetParam().refValue, SC_TEST_TOLERANCE);
}

int main(int argc, char **argv) {
  return MppTest(MppTestBuilder(argc, argv)
                     .WithConfigEntry("PDESolverPlotting", 1)
                     .WithConfigEntry("PDESolverVerbose", 0)
                     .WithConfigEntry("GeneratorVerbose", 0)
                     .WithConfigEntry("NewtonVerbose", 0)
                     .WithConfigEntry("LinearVerbose", 0)
                     .WithConfigEntry("ConfigVerbose", 0)
                     .WithConfigEntry("MeshVerbose", 0)
                     .WithConfigEntry("MainVerbose", 0)
                     .WithConfigEntry("MLMCVerbose", 0)
                     .WithConfigEntry("SCVerbose", 2)
                     .WithScreenLogging()
                     .WithPPM())
      .RUN_ALL_MPP_TESTS();
}