#include <map>
#include <STDGViscoAcousticAssemble.hpp>
#include "GMRES.hpp"
#include "SpaceTime.hpp"
#include "TestEnvironment.hpp"

class TestRiemannInterpolation : public Test {
protected:
  void SetUp() override {}

  void TearDown() override {}
};

TEST_F(TestRiemannInterpolation, TestInterpolation) {
  static constexpr int distributeLevel = 2;
  static constexpr int fineLevel = 5;
  Results results;
  auto meshes = MeshesCreator("QD").WithPLevel(distributeLevel).WithLevel(fineLevel).CreateUnique();
  for (short deg = 0; deg < 5; deg++) {
    DegreePair degree(deg, deg);
    auto assemble = std::make_unique<STDGViscoAcousticAssemble>(*meshes, degree, "DoubleRiemann");

    assemble->PrintInfo();
    for (int l = distributeLevel; l < fineLevel; l++) {
      mout.StartBlock("d" + std::to_string(deg) + "l" + std::to_string(l));
      assemble->AdaptQuadrature({l, l}, 6, true);

      Vector interpolated_solution(assemble->GetSharedDisc(), {l, l});

      assemble->get_exact_solution(interpolated_solution);
      results.append("deg", deg);
      results.append("l", l);
      results.append("L2", assemble->L2Error(interpolated_solution));
      results.append("L1", assemble->L1Error(interpolated_solution));

      interpolated_solution.Clear();

      assemble->get_projected_exact_solution(interpolated_solution);

      results.append("L2_proj", assemble->L2Error(interpolated_solution));
      results.append("L1_proj", assemble->L1Error(interpolated_solution));

      mout.EndBlock();
    }
  }
  results.PrintInfo();
}

int main(int argc, char **argv) {
  return MppTest(MppTestBuilder(argc, argv)
                     .WithScreenLogging()
                     .WithPPM()
                     .WithoutDefaultConfig()
                     .WithConfigEntry("Overlap", "STCellsWithCorners")
                     .WithConfigEntry("Distribution", "deformed_optimized")
                     .WithConfigEntry("Verbose", "1")
                     .WithConfigEntry("ConfigVerbose", "1")
                     .WithConfigEntry("LinearVerbose", "1")
                     .WithConfigEntry("Steps", "1000")
                     .WithConfigEntry("Epsilon", "1e-14")
                     .WithConfigEntry("Reduction", "1e-13")
                     .WithConfigEntry("adaptCellQuad", "1")
                     .WithConfigEntry("LinearVerbose", 0)
                     .WithConfigEntry("HighQuadDeg", "1")
                     .WithConfigEntry("adaptCellQuad", "1")
                     .WithConfigEntry("ResultsVerbose", 1)
                     .WithConfigEntry("useL2Projection", "1")
                     .WithConfigEntry("ResultsVerbose", "1"))
      .RUN_ALL_MPP_TESTS();
}
