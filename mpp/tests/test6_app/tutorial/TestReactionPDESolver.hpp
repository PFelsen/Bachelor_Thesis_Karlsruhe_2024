#ifndef TESTREACTIONMAIN_HPP
#define TESTREACTIONMAIN_HPP

#include "ReactionPDESolver.hpp"
#include "TestEnvironment.hpp"

class TestReactionPDESolver : public TestWithParam<std::pair<ConfigMap, ValueMap>> {
protected:
  std::unique_ptr<ReactionPDESolver> pdeSolver;

  double TEST_TOLERANCE = 1e-6;

  ConfigMap configMap{{"PDESolverPlotting", "1"}, {"PDESolverVerbose", "2"},
                      {"AssembleVerbose", "1"},   {"MeshesVerbose", "1"},
                      {"NewtonVerbose", "1"},     {"LinearVerbose", "0"},
                      {"ConfigVerbose", "1"},     {"MeshVerbose", "2"},
                      {"MainVerbose", "0"}};

  ValueMap refValues{};

  explicit TestReactionPDESolver(ConfigMap additionalMap1, ValueMap valueMap1 = {}) {
    additionalMap1.merge(configMap);
    ConfigMap additionalMap2 = GetParam().first;
    additionalMap2.merge(additionalMap1);
    Config::Initialize(additionalMap2);
    if (valueMap1.empty()) refValues = GetParam().second;
    else refValues = valueMap1;
    pdeSolver = std::make_unique<ReactionPDESolver>(PDESolverConfig());
  }

  void TestRun() {
    auto problem = CreateReactionProblemShared(PDESolverConfig().problemName);
    auto solution = pdeSolver->Run(problem);

    EXPECT_TRUE(solution.converged);

    pdeSolver->PrintValues(solution);

    for (auto &[name, refValue] : refValues) {
      EXPECT_NEAR(refValue, solution.values[name], TEST_TOLERANCE) << name;
    }
    EXPECT_NO_THROW(solution.vector.GetMesh().GetProcSets().CheckConsistency());
    EXPECT_NO_THROW(solution.vector.GetMatrixGraph().GetProcSets().CheckConsistency());
  }

  void TearDown() override {
    PPM->Barrier(0);
    Plotting::Instance().Clear();
    Config::Close();
  }
};

class TestProblemsPGReaction : public TestReactionPDESolver {
public:
  TestProblemsPGReaction() :
      TestReactionPDESolver({ConfigMap{{"Model", "PGReaction"}, {"rkorder", "-2"}}}) {}
};

class TestProblemsPGReactionFullCommSplit : public TestReactionPDESolver {};

#define TEST_PROBLEMS(TestClass)


INSTANTIATE_TEST_SUITE_P(
    TestReactionMainProgram, TestProblemsPGReaction,
    Values(std::pair{ConfigMap{{"Problem", "ExponentialReaction2D"},
                               {"level", "7"},
                               {"plevel", "2"},
                               {"degree", "1"},
                               {"Convection", "1.0"},
                               {"Diffusion", "0.01"},
                               {"Reaction", "1.0"},
                               {"T", "1.0"},
                               {"dt", "0.25"}},
                     ValueMap{
                        {"Energy", 0.},
                        {"InFlow", 0.},
                        {"Mass", 0.04055635},
                        {"Outflow", 0.140883},
                     }},

           std::pair{ConfigMap{{"Problem", "ExponentialReaction2D"},
                               {"level", "7"},
                               {"plevel", "2"},
                               {"degree", "1"},
                               {"Convection", "1.0"},
                               {"Diffusion", "0.01"},
                               {"Reaction", "0.0"},
                               {"T", "1.0"},
                               {"dt", "0.25"}},
                     ValueMap{
                        {"Energy", 0.},
                        {"InFlow", 0.},
                        {"Mass", 0.020919945},
                        {"Outflow", 0.060551511},
                     }},

           std::pair{ConfigMap{{"Problem", "ExponentialReaction2D"},
                               {"level", "7"},
                               {"plevel", "2"},
                               {"degree", "1"},
                               {"Convection", "0.5"},
                               {"Diffusion", "0.01"},
                               {"Reaction", "5.0"},
                               {"T", "1.0"},
                               {"dt", "0.25"}},
                     ValueMap{
                        {"Energy", 0.},
                        {"InFlow", 0.},
                        {"Mass", 0.79990508},
                        {"Outflow", 4.9827139},
                     }},

           std::pair{ConfigMap{{"Problem", "LogisticReaction2D"},
                               {"level", "7"},
                               {"plevel", "2"},
                               {"degree", "1"},
                               {"Convection", "1.0"},
                               {"Diffusion", "0.01"},
                               {"T", "1.0"},
                               {"dt", "0.25"},
                               {"Reaction0", "5"},
                               {"Reaction1", "5"}},
                     ValueMap{
                        {"Energy", 0.},
                        {"InFlow", 0.},
                        {"Mass", 0.20913255},
                        {"Outflow", 0.73577098},
                     }},

           std::pair{ConfigMap{{"Problem", "PollutionExponentialReactionSquare500"},
                               {"level", "0"},
                               {"plevel", "0"},
                               {"degree", "1"},
                               {"Convection", "1.0"},
                               {"Diffusion", "0.01"},
                               {"Reaction", "5.0"},
                               {"delta", "0.0"},
                               {"flux_alpha", "1"},
                               {"penalty", "25"},
                               {"sign", "1"},
                               {"T", "1.6"},
                               {"dt", "0.04"},
                               {"dt_min", "0.0001"}},
                     ValueMap{
                        {"Energy", 0.},
                        {"InFlow", 0.},
                        {"Mass", 0.016628757},
                        {"Outflow", 0.2066259},
                     }},

           std::pair{ConfigMap{{"Problem", "PollutionExponentialReactionSquare500"},
                               {"level", "0"},
                               {"plevel", "0"},
                               {"degree", "1"},
                               {"Convection", "1.0"},
                               {"Diffusion", "0.01"},
                               {"Reaction", "5.0"},
                               {"delta", "0.0"},
                               {"flux_alpha", "1"},
                               {"penalty", "25"},
                               {"sign", "1"},
                               {"T", "1.6"},
                               {"dt", "0.04"},
                               {"dt_min", "0.0001"},
                               {"Diffusion", "0.0001"}},
                     ValueMap{
                        {"Energy", 0.},
                        {"InFlow", 0.},
                        {"Mass", 0.016628757},
                        {"Outflow", 0.2066259},
                     }},

           std::pair{ConfigMap{{"Problem", "PollutionExponentialReactionSquare500"},
                               {"level", "0"},
                               {"plevel", "0"},
                               {"degree", "1"},
                               {"Convection", "1.0"},
                               {"Diffusion", "0.01"},
                               {"Reaction", "5.0"},
                               {"delta", "0.0"},
                               {"flux_alpha", "1"},
                               {"penalty", "25"},
                               {"sign", "1"},
                               {"T", "1.6"},
                               {"dt", "0.04"},
                               {"dt_min", "0.0001"},
                               {"Diffusion", "0.0001"},
                               {"delta", "1.0"}},
                     ValueMap{
                        {"Energy", 0.},
                        {"InFlow", 0.},
                        {"Mass", 0.016628757},
                        {"Outflow", 0.2066259},
                     }},

           std::pair{ConfigMap{{"Problem", "PollutionExponentialReactionSquare500"},
                               {"level", "0"},
                               {"plevel", "0"},
                               {"degree", "1"},
                               {"Convection", "1.0"},
                               {"Diffusion", "0.01"},
                               {"Reaction", "5.0"},
                               {"delta", "0.0"},
                               {"flux_alpha", "1"},
                               {"penalty", "25"},
                               {"sign", "1"},
                               {"T", "1.6"},
                               {"dt", "0.04"},
                               {"dt_min", "0.0001"},
                               {"Model", "DGReaction"}},
                     ValueMap{
// Todo this test has ambiguous results. Problem is probably caused in DGReactionAssemble
//                        {"Energy", 0.},
//                        {"InFlow", 0.},
//                        {"Mass", 0.3767317745189242},
//                        {"Outflow", 1.6762015201506151},
                     }},

           std::pair{ConfigMap{{"Problem", "PollutionLogisticReactionSquare500"},
                               {"level", "0"},
                               {"plevel", "0"},
                               {"degree", "1"},
                               {"Convection", "1.0"},
                               {"Diffusion", "0.01"},
                               {"delta", "0.0"},
                               {"flux_alpha", "1"},
                               {"penalty", "25"},
                               {"sign", "1"},
                               {"T", "1.6"},
                               {"dt", "0.04"},
                               {"dt_min", "0.0001"},
                               {"Reaction0", "5"},
                               {"Reaction1", "5"}},
                     ValueMap{
                        {"Energy", 0.},
                        {"InFlow", 0.},
                        {"Mass", 0.016628757},
                        {"Outflow", 0.2066259},
                     }}));

TEST_P(TestProblemsPGReaction, TestRun) { TestRun(); }

#endif // TESTREACTIONMAIN_HPP
