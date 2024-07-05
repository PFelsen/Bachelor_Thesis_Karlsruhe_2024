#include "TestSpaceTimeMesh.hpp"

TEST_P(TestSpaceTimeMesh, TestCellPoints) { checkCellPoints(); }

TEST_P(TestSpaceTimeMesh, TestOverlapPoints) { checkOverlapPoints(); }

int main(int argc, char **argv) {
  return MppTest(MppTestBuilder(argc, argv)
                     .WithConfigEntry("DistributionVerbose", 1)
                     .WithConfigEntry("MeshesVerbose", 4)
                     .WithConfigEntry("MeshVerbose", 2)
                     .WithConfigEntry("Distribution", "deformed_optimized")
                     .WithConfigEntry("Overlap", "STCellsWithCorners")
                     .WithParallelListeners()
                     .WithScreenLogging()
                     .WithPPM())
      .RUN_ALL_MPP_TESTS();
}