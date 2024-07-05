#include "ConvergenceStudy.hpp"
#include "STAcousticMain.hpp"
#include "TestEnvironment.hpp"

#include <map>
#include <vector>

class ViscoAcousticDGTMarmousi : public Test {
protected:
  void SetUp() override {
    std::map<std::string, string> configMap(
        {{"plevel", "2"},
         {"level", "3"},
         {"ConformingReconstruction", "0"},
         {"LinearPrintSteps", "1"},
         {"degree", "1"},
         {"time_degree", "1"},
         {"numL", "3"},
         {"AdaptQuadrature", "1"},
         {"load_balancing", "1"},
         {"ExcludedResults", "||u_proj-u_h||_DG, EE, Conf, DGNorm"},
         {"Overlap", "STCellsWithFaces"},
         {"PreconditionerVerbose", "1"},
         {"BaseSolverVerbose", "1"},
         {"Reduction", "1e-10"},
         {"LinearSolver", "GMRES"},
         {"ModelImageMu", "marmousi2-vs"},
         {"MinMu", "1e-3"},
         {"MaxMu", "2.802"},
         {"ProblemLevel", "2"},
         {"Preconditioner", "PointBlockGaussSeidel"},
         {"restart", "100"},
         {"ModelImageKappa", "marmousi2-vp"},
         {"MaxRho", "2626.999855041504"},
         {"pml", "1"},
         {"ModelImageRho", "marmousi2-density"},
         {"LinearVerbose", "1"},
         {"Distribution", "deformed_optimized"},
         {"Epsilon", "1e-10"},
         {"LinearEpsilon", "1e-5"},
         {"LinearReduction", "1e-6"},
         {"LinearSteps", "1000"},
         {"Verbose", "10"},
         {"roi_min", "4.75, 0.1, 4.0"},
         {"roi_max", "7.25, 0.4, 4.0"},
         {"truncateSTMesh", "1"},
         {"Mesh", "ST_marmousi2_squares_abgeschnitten"},
         {"Model", "STDGViscoAcousticAssemble"},
         {"Problem", "Marmousi2_rhs"},
         {"source_duration", "0.1"},
         {"ProblemMid", "3.0, 0.25, 0.15"},
         {"MinRho", "1009.9992752075195"},
         {"linear_weight", "1"},
         {"MinKappa", "1.0279998779296875"},
         {"ConfigVerbose", "1"},
         {"MaxKappa", "4.7"},
         {"DebugLevel", "0"},
         {"GoalFunctional", "linear"},
         {"ResultsVerbose", "100"}});

    Config::Initialize(configMap);
  }

  void TearDown() override {}
};

TEST_F(ViscoAcousticDGTMarmousi, BenchmarkTest) {
  PDESolverConfig mainConf;
  mainConf.modelName = "STAcoustic";
  ConvergenceStudySTAcoustic main(mainConf);
  main.Method();
  if (PPM->Master()) {
    vector<PrintInfoEntry<vector<double>>> values;
    for (auto &[key, value] : convertResults(main.GetResults())) {
      values.push_back({key, value, 0});
    }
    mout.PrintInfo("Values", 1, values);
    std::vector<double> res = convertResults(main.GetResults())["GoalFunctional"];
    EXPECT_VECTOR_NEAR(res, {0.010192, 0.010782}, 1e-4);
  }
}

int main(int argc, char **argv) {
  MppTest mppTest = MppTestBuilder(argc, argv).WithScreenLogging().WithPPM().WithoutDefaultConfig();
  return mppTest.RUN_ALL_MPP_TESTS();
}
