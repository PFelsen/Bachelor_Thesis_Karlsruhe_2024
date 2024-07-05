#include "TestConformingInterpolation.hpp"

int main(int argc, char **argv) {
  return MppTest(MppTestBuilder(argc, argv)
                     .WithScreenLogging()
                     .WithConfigEntry("Overlap", "STCellsWithFaces")
                     .WithConfigEntry("OverlapVerbose", 0)
                     .WithConfigEntry("Distribution", "deformed_optimized")
                     .WithConfigEntry("ExcludedResults", "EE")
                     .WithConfigEntry("MeshVerbose", 100)
                     .WithConfigEntry("ConfigVerbose", 0)
                     .WithPPM())
      .RUN_ALL_MPP_TESTS();
}