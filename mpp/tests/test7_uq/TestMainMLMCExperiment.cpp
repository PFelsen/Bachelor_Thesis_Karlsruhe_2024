#include "Main.hpp"
#include "TestEnvironment.hpp"

TEST(TestMLMCExperiment, TestRun) { MLMCExperiment(); }

int main(int argc, char **argv) {
  return MppTest(MppTestBuilder(argc, argv)
                     .WithConfigEntry("Problem", "StochasticLaplace2D")
                     .WithConfigEntry("Model", "LagrangeElliptic")
                     .WithConfigEntry("Quantity", "L2")
                     .WithConfigEntry("degree", 1)
                     .WithConfigEntry("initSamples", "[32, 16, 8]")
                     .WithConfigEntry("initLevels", "[3, 4, 5]")
                     .WithConfigEntry("WTime", "00:02:00")
                     .WithConfigEntry("StochasticField", "LogNormal")
                     .WithConfigEntry("Mean", 0.0)
                     .WithConfigEntry("sigma", 1.0)
                     .WithConfigEntry("lambda", "[0.15, 0.15]")
                     .WithConfigEntry("smoothing", 1.8)
                     .WithConfigEntry("VtuPlot", 0)
                     .WithConfigEntry("ParallelPlotting", 0)
                     .WithConfigEntry("MainVerbose", 1)
                     .WithConfigEntry("MeshVerbose", 0)
                     .WithConfigEntry("ConfigVerbose", 1)
                     .WithConfigEntry("LinearVerbose", 0)
                     .WithConfigEntry("NewtonVerbose", 0)
                     .WithConfigEntry("MSFEMVerbose", 0)
                     .WithConfigEntry("PDESolverVerbose", 0)
                     .WithConfigEntry("SLEstimatorVerbose", 1)
                     .WithConfigEntry("MLEstimatorVerbose", 1)
                     .WithConfigEntry("CirculantEmbeddingVerbose", 0)
                     .WithScreenLogging()
                     .WithRandomInitialized()
                     .WithPPM())
      .RUN_ALL_MPP_TESTS();
}