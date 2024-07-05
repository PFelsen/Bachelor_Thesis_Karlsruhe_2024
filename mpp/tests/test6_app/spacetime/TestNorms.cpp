#include "TestNorms.hpp"
#include "MeshesCreator.hpp"
#include "STAssemble.hpp"
#include "STDGViscoAcousticAssemble.hpp"

TEST_P(TestNorms, Test) {
  Config::PrintInfo();
  std::unique_ptr<Meshes> meshes =
      MeshesCreator("SpaceTimeSquare").WithLevel(4).WithPLevel(2).CreateUnique();
  auto zeroProblem = std::shared_ptr<AcousticProblem>(CreateAcousticProblemUnique("Zero"));
  auto assemble =
      std::make_unique<STDGViscoAcousticAssemble>(*meshes, DegreePair(1, 1), zeroProblem);
  Vector v(1.0, assemble->GetSharedDisc());

  EXPECT_DOUBLE_EQ(assemble->L2Error(v), assemble->L2Norm(v));

  EXPECT_NEAR(assemble->L2Norm(v), sqrt(3), 1e-10);

  fillSpaceTimeVector(v, assemble->GetProblem().Dim(),
                      [](const Point &x, const cell &, int) { return x[0] + x[1] + x.t(); });
  EXPECT_DOUBLE_EQ(assemble->L2Error(v), assemble->L2Norm(v));
  EXPECT_NEAR(assemble->L2Norm(v), sqrt(7.5), 1e-10);
}

int main(int argc, char **argv) {
  return MppTest(MppTestBuilder(argc, argv)
                     .WithScreenLogging()
                     .WithoutDefaultConfig()
                     .WithConfigEntry("Distribution", "deformed_optimized")
                     .WithConfigEntry("MatrixGraphVerbose", 0)
                     .WithConfigEntry("AssembleVerbose", 3)
                     .WithConfigEntry("ResultsVerbose", 1)
                     .WithConfigEntry("ConfigVerbose", 0)
                     .WithConfigEntry("LinearVerbose", 1)
                     .WithConfigEntry("MeshVerbose", 0)
                     .WithConfigEntry("MainVerbose", 0)
                     .WithConfigEntry("Overlap", "STCellsWithCorners")
                     .WithConfigEntry("vtkplot", 1)
                     .WithConfigEntry("Steps", 1000)
                     .WithConfigEntry("numL", 0)
                     .WithPPM())
      .RUN_ALL_MPP_TESTS();
}
