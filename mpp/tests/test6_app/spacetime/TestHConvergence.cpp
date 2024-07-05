#include "ConvergenceStudy.hpp"
#include "PrintUtil.hpp"
#include "STAcousticMain.hpp"
#include "TestEnvironment.hpp"

#include <map>
#include <memory>

#include <string>

struct ConvergenceParamType {
  std::string model;
  std::string problem;
  std::map<std::string, std::vector<double>> values;
};

std::ostream &operator<<(std::ostream &s, const ConvergenceParamType &param) {
  return s << "Model:" << param.model << " Problem:" << param.problem;
}

class TestHConvergence : public TestWithParam<ConvergenceParamType> {
protected:
  void SetUp() override {}

  void TearDown() override {}

  void run() {

    auto mainConfig = PDESolverConfig();

    mainConfig.level = 4;
    mainConfig.pLevel = 2;
    mainConfig.modelName = GetParam().model;
    mainConfig.degree = 1;
    mainConfig.timeDegree = 1;
    mainConfig.problemName = GetParam().problem;
    mainConfig.mesh = "ST_squ_simple";
    mainConfig.preconditioner = "PointBlockGaussSeidel";

    auto convergence = CreateConvergenceStudy(mainConfig);
    convergence->Method();
    auto results = convergence->GetResults();
    auto convertedResults = convertResults(results);
    for (auto &[name, values] : GetParam().values) {
      PrintValues(name, convertedResults[name], 100, 100);
      EXPECT_VECTOR_NEAR(convertedResults[name], GetParam().values.at(name), 1e-5);
    }
  }
};

INSTANTIATE_TEST_SUITE_P(
    TestHConvergence, TestHConvergence,
    Values(ConvergenceParamType{.model = "STAcoustic",
                                .problem = "SinCos",
                                .values = {{"L2_Error", {0.046296, 0.010290, 0.002409}}}},
           ConvergenceParamType{.model = "STAcousticGL",
                                .problem = "SinCos",
                                .values = {{"L2_Error", {0.192547, 0.060539, 0.016697}}}},
           ConvergenceParamType{.model = "STAcousticMF",
                                .problem = "SinCos",
                                .values = {{"L2_Error", {0.046296, 0.010290, 0.002409}}}},
           ConvergenceParamType{.model = "STAcoustic",
                                .problem = "QuadraticST",
                                .values = {{"L2_Error", {0.005981, 0.001497, 0.000375}}}},
           ConvergenceParamType{.model = "STAcousticGL",
                                .problem = "QuadraticST",
                                .values = {{"L2_Error", {0.048454, 0.012791, 0.003282}}}},
           ConvergenceParamType{.model = "STAcousticMF",
                                .problem = "QuadraticST",
                                .values = {{"L2_Error", {0.005981, 0.001497, 0.000375}}}},
           ConvergenceParamType{.model = "STAcoustic",
                                .problem = "RiemannST",
                                .values = {{"L2_Error", {0.033486, 0.026174, 0.020456}}}},
           ConvergenceParamType{.model = "STAcousticGL",
                                .problem = "RiemannST",
                                .values = {{"L2_Error", {0.043172, 0.035367, 0.028361}}}},
           ConvergenceParamType{.model = "STAcousticMF",
                                .problem = "RiemannST",
                                .values = {{"L2_Error", {0.033486, 0.026174, 0.020456}}}}));

TEST_P(TestHConvergence, runTest) { this->run(); }

int main(int argc, char **argv) {
  return MppTest(MppTestBuilder(argc, argv)
                     .WithConfigEntry("Distribution", "deformed_optimized")
                     .WithConfigEntry("AdaptQuadrature", 0)
                     .WithConfigEntry("DistributionVerbose", 0)
                     .WithConfigEntry("LinearReduction", 1e-6)
                     .WithConfigEntry("LinearEpsilon", 1e-8)
                     .WithConfigEntry("ConfigVerbose", 5)
                     .WithConfigEntry("LinearVerbose", 0)
                     .WithConfigEntry("LinearSteps", 200)
                     .WithConfigEntry("MeshesVerbose", 0)
                     .WithConfigEntry("ResultsVerbose", 0)
                     .WithConfigEntry("MeshVerbose", 0)
                     .WithConfigEntry("Overlap", "STCellsWithFaces")
                     .WithConfigEntry("ExcludedResults",
                                      "EE, DGNorm, DGError, DG_Int_Error, DGNorm_Conf, "
                                      "goal_functional, ||u_proj-u_h||_DG, Conf")
                     .WithScreenLogging()
                     .WithPPM())
      .RUN_ALL_MPP_TESTS();
}