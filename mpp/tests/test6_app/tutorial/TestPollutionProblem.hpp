#ifndef TESTREACTIONMAIN_HPP
#define TESTREACTIONMAIN_HPP

#include "TestEnvironment.hpp"

#include "STTransportMain.hpp"
#include "TransportPDESolver.hpp"

template<typename PDESolverType>
class TestPollutionProblem : public TestWithParam<std::pair<ConfigMap, ValueMap>> {
protected:
  std::unique_ptr<PDESolverType> pdeSolver;

  double TEST_TOLERANCE = 1e-5;

  ConfigMap configMap{{"PDESolverPlotting", "1"}, {"PDESolverVerbose", "2"},
                      {"AssembleVerbose", "1"},   {"MeshesVerbose", "1"},
                      {"NewtonVerbose", "1"},     {"LinearVerbose", "0"},
                      {"ConfigVerbose", "1"},     {"MeshVerbose", "2"},
                      {"MainVerbose", "0"}};

  ValueMap refValues{};

  // Todo constructor can be cleaned up. AddtionalMap not needed
  TestPollutionProblem(ConfigMap additionalMap1, ValueMap valueMap1 = {}) {
    additionalMap1.merge(configMap);
    ConfigMap additionalMap2 = GetParam().first;
    additionalMap2.merge(additionalMap1);
    Config::Initialize(additionalMap2);
    if (valueMap1.empty()) refValues = GetParam().second;
    else refValues = valueMap1;
    pdeSolver = std::make_unique<PDESolverType>(PDESolverConfig());
  }

  void TestRun() {
    auto problem = CreateTransportProblemShared(PDESolverConfig().problemName);
    auto solution = pdeSolver->Run(problem);
    for (auto &refValue : refValues) {
      EXPECT_NEAR(refValue.second, solution.values[refValue.first], TEST_TOLERANCE);
    }
    EXPECT_NO_THROW(solution.vector.GetMatrixGraph().GetProcSets().CheckConsistency());
    EXPECT_NO_THROW(solution.vector.GetMesh().GetProcSets().CheckConsistency());
  }

  void TearDown() override {
    PPM->Barrier(0);
    Plotting::Instance().Clear();
    Config::Close();
  }
};

class TestSpaceTimePollutionProblem : public TestPollutionProblem<STTransportMain> {
public:
  TestSpaceTimePollutionProblem() : TestPollutionProblem(ConfigMap{{"Model", "STTransport"}}) {}
};

class TestTimeSteppingPollutionProblem : public TestPollutionProblem<TransportPDESolver> {
public:
  TestTimeSteppingPollutionProblem() : TestPollutionProblem(ConfigMap{{"Model", "DGTransport"}}) {}
};

#endif // TESTREACTIONMAIN_HPP
