#ifndef TESTACOUSTICMAIN_HPP
#define TESTACOUSTICMAIN_HPP

#include "AcousticPDESolver.hpp"
#include "TestEnvironment.hpp"

class TestAcousticPDESolver : public TestWithParam<std::pair<ConfigMap, ValueMap>> {
protected:
  std::unique_ptr<AcousticPDESolver> pdeSolver;

  double TEST_TOLERANCE = 1e-6;

  ConfigMap configMap{{"TimeIntegratorVerbose", "1"}, {"PDESolverPlotting", "1"},
                      {"PDESolverVerbose", "2"},      {"AssembleVerbose", "1"},
                      {"MeshesVerbose", "1"},         {"NewtonVerbose", "1"},
                      {"LinearVerbose", "0"},         {"ConfigVerbose", "1"},
                      {"MeshVerbose", "2"},           {"MainVerbose", "0"}};

  ValueMap refValues{};

  int commSplit;

  TestAcousticPDESolver(ConfigMap additionalMap1, ValueMap valueMap1 = {}, int commSplit = 0) :
      commSplit(commSplit) {
    additionalMap1.merge(configMap);
    ConfigMap additionalMap2 = GetParam().first;
    additionalMap2.merge(additionalMap1);
    Config::Initialize(additionalMap2);
    if (valueMap1.empty()) refValues = GetParam().second;
    else refValues = valueMap1;
    pdeSolver = std::make_unique<AcousticPDESolver>(PDESolverConfig());
  }

  void TestRun() {
    auto problem = CreateAcousticProblemShared(PDESolverConfig().problemName);
    auto solution = pdeSolver->Run(problem, commSplit);
    EXPECT_TRUE(solution.converged);
    for (auto &refValue : refValues) {
      EXPECT_NEAR(refValue.second, solution.values[refValue.first], TEST_TOLERANCE);
    }
  }

  void TearDown() override {
    PPM->Barrier(0);
    Plotting::Instance().Clear();
    Config::Close();
  }
};

class TestProblems : public TestAcousticPDESolver {
public:
  TestProblems() : TestAcousticPDESolver({ConfigMap{{"Model", "DGAcoustic"}}}) {}
};

class TestProblemsFullCommSplit : public TestAcousticPDESolver {
public:
  TestProblemsFullCommSplit() :
      TestAcousticPDESolver({ConfigMap{{"Model", "DGAcoustic"}}}, {}, 1) {}
};

#define TEST_PROBLEMS(TestClass)                                                                   \
                                                                                                   \
  INSTANTIATE_TEST_SUITE_P(TestAcousticPDESolver, TestClass,                                       \
                           Values(std::pair{ConfigMap{{"Problem", "Linear"},                       \
                                                      {"degree", "1"},                             \
                                                      {"T", "0"},                                  \
                                                      {"rkorder", "2"}},                           \
                                            ValueMap{{"L2Error", 0.0}}},                           \
                                  std::pair{ConfigMap{{"Problem", "Quadratic"},                    \
                                                      {"degree", "2"},                             \
                                                      {"T", "0"},                                  \
                                                      {"rkorder", "2"}},                           \
                                            ValueMap{{"L2Error", 0.0}}},                           \
                                  std::pair{ConfigMap{{"Problem", "RiemannWave2D"},                \
                                                      {"rkorder", "-2"},                           \
                                                      {"level", "3"},                              \
                                                      {"plevel", "3"},                             \
                                                      {"degree", "1"},                             \
                                                      {"dt", "0.020"},                             \
                                                      {"T", "1.0"}},                               \
                                            ValueMap{{"L2Error", 1.1308499},                       \
                                                     {"Energy", 1.3502834}}},                      \
                                  std::pair{ConfigMap{{"Problem", "CRC"},                          \
                                                      {"rkorder", "-2"},                           \
                                                      {"level", "0"},                              \
                                                      {"degree", "1"},                             \
                                                      {"dt", "0.04"},                              \
                                                      {"T", "8.0"}},                               \
                                            ValueMap{{"NormCoefficients", 6.7790422021192223}}},   \
                                  std::pair{ConfigMap{{"Problem", "GaussHatAndRicker2D"},          \
                                                      {"rkorder", "-2"},                           \
                                                      {"level", "4"},                              \
                                                      {"plevel", "4"},                             \
                                                      {"dt", "0.03125"},                           \
                                                      {"degree", "2"}},                            \
                                            ValueMap{{"L2", 0.18414399}}}));                       \
                                                                                                   \
  TEST_P(TestClass, TestRun) { TestRun(); }

#endif // TESTACOUSTICMAIN_HPP
