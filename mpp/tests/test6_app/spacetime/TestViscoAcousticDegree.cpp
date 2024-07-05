#include "AdaptiveStudy.hpp"
#include "ConvergenceStudy.hpp"
#include "STAcousticMain.hpp"
#include "TestEnvironment.hpp"

#include <map>
#include <vector>

struct TestData {
  vector<double> expected_goal_functional;
};

class ViscoAcousticMarmousi : public TestWithParam<TestData> {
protected:
  vector<double> expected_goal_functional;

  AdaptiveStudySTAcousticPG adaptiveStudy;

  ViscoAcousticMarmousi() :
      expected_goal_functional(GetParam().expected_goal_functional),
      adaptiveStudy(PDESolverConfig()
                        .WithModel("STAcousticPG")
                        .WithPLevel(2)
                        .WithLevel(2)
                        .WithDegree({1, 1})
                        .WithProblem("RiemannST")
                        .WithMesh("ST_squ_simple"),
                    AdaptiveData{.refinement_steps = 2,
                                 .errorEstimator = "Dual",
                                 .theta = 0.1,
                                 .theta_min = 1e-3,
                                 .refine_by = "abs_value"}) {}

  void SetUp() override {}

  void TearDown() override {}
};

TEST_P(ViscoAcousticMarmousi, AdaptiveTest) {
  adaptiveStudy.Method();
  auto results = adaptiveStudy.GetResults();
  if (PPM->Master(0)) {
    std::vector<double> res = convertResults(results)["GoalFunctional"];
    EXPECT_EQ(res.size(), 3) << "Size of result has to match with refinements_steps + 1";
    for (int i = 0; i < 2; ++i) {
      double exp = expected_goal_functional[i];
      EXPECT_LT(std::abs((res[i] - exp) / exp), 1e-2)
          << "relative error has to be less than tolerance";
    }
  }
}

INSTANTIATE_TEST_SUITE_P(ViscoAcousticMarmousiTest, ViscoAcousticMarmousi,
                         Values(TestData{{0.502611, 0.504012, 0.503739}}));

int main(int argc, char **argv) {
  return MppTest(MppTestBuilder(argc, argv)
                     .WithScreenLogging()
                     .WithoutDefaultConfig()
                     .WithConfigEntry("AdaptQuadrature", "0")
                     .WithConfigEntry("vtkplot", "1")
                     .WithConfigEntry("GoalFunctional", "quadratic")
                     .WithConfigEntry("dual_functional", "quadratic")
                     .WithConfigEntry("roi_min", "-0.5, 0.0, 1.0")
                     .WithConfigEntry("roi_max", "0.5, 1.0, 1.0")
                     .WithConfigEntry("LinearSolver", "GMRES")
                     .WithConfigEntry("LinearReduction", "1e-8")
                     .WithConfigEntry("LinearEpsilon", "1e-10")
                     .WithConfigEntry("LinearSteps", "200")
                     .WithConfigEntry("LinearPrintSteps", "5")
                     .WithConfigEntry("MatrixGraphVerbose", 0)
                     .WithConfigEntry("AssembleVerbose", 3)
                     .WithConfigEntry("ResultsVerbose", 1)
                     .WithConfigEntry("ConfigVerbose", 5)
                     .WithConfigEntry("LinearVerbose", 1)
                     .WithConfigEntry("MeshVerbose", 50)
                     .WithConfigEntry("MeshesVerbose", 1)
                     .WithConfigEntry("MainVerbose", 0)
                     .WithConfigEntry("TimeLevel", 5)
                     .WithConfigEntry("Distribution", "deformed_optimized")
                     .WithConfigEntry("Overlap", "STCellsWithCorners")
                     .WithConfigEntry("Steps", 1000)
                     .WithConfigEntry("numL", 0)
                     .WithPPM())
      .RUN_ALL_MPP_TESTS();
}
