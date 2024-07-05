#include "ConvergenceStudy.hpp"
#include "STAcousticMain.hpp"
#include "TestEnvironment.hpp"

#include <map>
#include <vector>

class ViscoAcousticMarmousi : public Test {
protected:
  void SetUp() override {
    std::map<string, string> configMap({
        {"LinearPrintSteps", "1"},
        {"Overlap", "STCellsWithCorners"},
        {"LinearSolver", "GMRES"},
        {"ProblemLevel", "2"},
        {"Preconditioner", "PointBlockGaussSeidel"},
        {"ModelImageKappa", "marmousi2-vp"},
        {"pml", "1.0"},
        {"ModelImageRho", "marmousi2-density"},
        {"Distribution", "deformed_optimized"},
        {"LinearEpsilon", "1e-5"},
        {"LinearReduction", "1e-6"},
        {"LinearSteps", "1000"},
        {"truncateSTMesh", "1"},
        {"source_duration", "0.1"},
        {"ProblemMid", "3.0, 0.25, 0.15"},
        {"MinRho", "1009.9992752075195"},
        {"MaxRho", "2626.999855041504"},
        {"linear_weight", "1"},
        {"MinKappa", "1.0279998779296875"},
        {"MaxKappa", "4.7"},
        {"ModelImageMu", "marmousi2-vs"},
        {"mlLength", "1.0"},
        {"useL2Proj", "0"},
        {"roi_min", "4.75, 0.1, 4.0"},
        {"roi_max", "7.25, 0.4, 4.0"},
        {"GoalFunctional", "linear"},
        {"LinearVerbose", "1"},
        {"ConfigVerbose", "0"},
    });

    Config::Initialize(configMap);
  }

  void TearDown() override {}
};

TEST_F(ViscoAcousticMarmousi, BenchmarkTest) {
  PDESolverConfig mainConf;
  mainConf.pLevel = 2;
  mainConf.level = 3;
  mainConf.timeDegree = 1;
  mainConf.degree = 1;
  mainConf.mesh = "ST_marmousi2_squares_abgeschnitten";
  mainConf.problemName = "Marmousi2_rhs";
  mainConf.modelName = "STAcousticPG";
  auto convergenceMethod = ConvergenceStudySTAcousticPG(mainConf);
  convergenceMethod.Method();
  if (PPM->Master(0)) {
    std::vector<double> res = convertResults(convergenceMethod.GetResults())["GoalFunctional"];
    ASSERT_NEAR(res[1], 0.001928574092444563, 1e-4);
  }
}

int main(int argc, char **argv) {
  MppTest mppTest = MppTestBuilder(argc, argv).WithScreenLogging().WithPPM().WithoutDefaultConfig();
  return mppTest.RUN_ALL_MPP_TESTS();
}
