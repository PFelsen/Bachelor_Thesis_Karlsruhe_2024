#ifndef TESTMARMOUSI_HPP
#define TESTMARMOUSI_HPP

#include "AcousticProblems.hpp"
#include "TestEnvironment.hpp"

TEST(TestMarmousi, Test1) {
  static constexpr int distributeLevel = 2;
  static constexpr int fineLevel = 4;
  auto meshes = MeshesCreator()
                    .WithPLevel(distributeLevel)
                    .WithLevel(fineLevel)
                    .WithMeshName("ST_marmousi2_squares_abgeschnitten")
                    .CreateUnique();

  auto oldProblem = CreateAcousticProblem("Marmousi2_rhs");
  auto newProblem = CreateAcousticProblem("Marmousi2_rhsIMAGE");

  for (int i = distributeLevel; i <= fineLevel; i++) {
    const Mesh &mesh = meshes->operator[]({i, i});
    for (cell c = mesh.cells(); c != mesh.cells_end(); c++) {
      EXPECT_DOUBLE_EQ(oldProblem->Rho(*c, c()), newProblem->Rho(*c, c()))
          << "Error at " << c() << " on Lvl " << mesh.Level();
    }

    for (cell c = mesh.cells(); c != mesh.cells_end(); c++) {
      EXPECT_DOUBLE_EQ(oldProblem->Kappa(*c, c()), newProblem->Kappa(*c, c()))
          << "Error at " << c() << " on Lvl " << mesh.Level();
    }
  }
}

int main(int argc, char **argv) {
  return MppTest(MppTestBuilder(argc, argv)
                     .WithScreenLogging()
                     .WithConfigEntry("ConfigVerbose", 0)
                     .WithConfigEntry("ProblemVerbose", 10)
                     .WithConfigEntry("ProblemLevel", 2)
                     .WithConfigEntry("ModelImageKappa", "marmousi2-vp")
                     .WithConfigEntry("MaxKappa", 4.7)
                     .WithConfigEntry("MinKappa", 1.0279998779296875)
                     .WithConfigEntry("ModelImageRho", "marmousi2-density")
                     .WithConfigEntry("MaxRho", 2626.999855041504)
                     .WithConfigEntry("MinRho", 1009.9992752075195)
                     .WithPPM())
      .RUN_ALL_MPP_TESTS();
}


#endif // TESTMARMOUSI_HPP
