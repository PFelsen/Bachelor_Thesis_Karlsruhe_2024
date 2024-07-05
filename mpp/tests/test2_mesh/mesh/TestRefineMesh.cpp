#include "TestRefineMesh.hpp"

#include "Mesh.hpp"

TEST_P(TestRefineMesh, TestRefineMesh) {
  const Mesh &m0 = meshes0->fine();
  const Mesh m1 = Mesh(m0, true);
  meshes0->PrintInfo();
  EXPECT_TRUE(areEqual(m1, meshes1->fine()));
}

int main(int argc, char **argv) {
  return MppTest(MppTestBuilder(argc, argv)
                     .WithConfigEntry("DistributionVerbose", 1)
                     .WithConfigEntry("MeshesVerbose", 4)
                     .WithConfigEntry("MeshVerbose", 2)
                     .WithParallelListeners()
                     .WithScreenLogging()
                     .WithPPM())
      .RUN_ALL_MPP_TESTS();
}