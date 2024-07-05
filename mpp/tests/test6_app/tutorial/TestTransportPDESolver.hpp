#ifndef TESTTRANSPORTMAIN_HPP
#define TESTTRANSPORTMAIN_HPP

#include <utility>

#include <PDESolver.hpp>
#include <TransportPDESolver.hpp>

#include "TestEnvironment.hpp"

class TestTransportPDESolver : public TestWithParam<std::pair<ConfigMap, ValueMap>> {
protected:
  std::unique_ptr<TransportPDESolver> pdeSolver;

  double TEST_TOLERANCE = 1e-6;

  ConfigMap configMap{{"PDESolverPlotting", "1"}, {"PDESolverVerbose", "2"},
                      {"AssembleVerbose", "1"},   {"MeshesVerbose", "1"},
                      {"NewtonVerbose", "1"},     {"LinearVerbose", "0"},
                      {"ConfigVerbose", "1"},     {"MeshVerbose", "2"},
                      {"MainVerbose", "0"}};

  ValueMap refValues{};

  int commSplit;

  TestTransportPDESolver(ConfigMap additionalMap1, ValueMap valueMap1 = {}, int commSplit = 0) :
      commSplit(commSplit) {
    additionalMap1.merge(configMap);
    ConfigMap additionalMap2 = GetParam().first;
    additionalMap2.merge(additionalMap1);
    Config::Initialize(additionalMap2);
    if (valueMap1.empty()) refValues = GetParam().second;
    else refValues = valueMap1;
    pdeSolver = std::make_unique<TransportPDESolver>(PDESolverConfig());
  }

  void TestRun() {
    auto problem = CreateTransportProblemShared(PDESolverConfig().problemName);
    auto solution = pdeSolver->Run(problem, commSplit);
    EXPECT_TRUE(solution.converged);
    pdeSolver->PrintValues(solution);
    for (auto &[name, refValue] : refValues) {
      EXPECT_NEAR(refValue, solution.values[name], TEST_TOLERANCE) << name;
    }
  }

  void TearDown() override {
    PPM->Barrier(0);
    Plotting::Instance().Clear();
    Config::Close();
  }
};

class TestProblems : public TestTransportPDESolver {
public:
  TestProblems() :
      TestTransportPDESolver({ConfigMap{{"Model", "DGTransport"}, {"rkorder", "-2"}}}) {}
};

class TestProblemsCommSplit : public TestTransportPDESolver {
public:
  TestProblemsCommSplit() :
      TestTransportPDESolver({ConfigMap{{"Model", "DGTransport"}, {"rkorder", "-2"}}}, {}, 1) {}
};

#endif // TESTTRANSPORTMAIN_HPP
