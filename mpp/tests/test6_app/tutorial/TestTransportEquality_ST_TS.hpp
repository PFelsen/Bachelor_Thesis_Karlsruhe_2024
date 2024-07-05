#ifndef MPP_TESTTRANSPORTEQUALITY_ST_TS_HPP
#define MPP_TESTTRANSPORTEQUALITY_ST_TS_HPP

#include "STTransportMain.hpp"
#include "TestEnvironment.hpp"
#include "TransportPDESolver.hpp"

class TestTransportEquality_ST_TS : public TestWithParam<std::pair<ConfigMap, ValueMap>> {
  std::unique_ptr<TransportPDESolver> mainTS;
  std::unique_ptr<STTransportMain> mainST;
  double TEST_TOLERANCE = 1e-10;
public:
  ValueMap refValues{};

  TestTransportEquality_ST_TS() {
    PDESolverConfig conf{};
    conf.degree = 2;
    conf.timeDegree = 0;
    conf.problemName = "PollutionSquare500";
    conf.modelName = "DGTransport";
    conf.rkorder = -1;
    mainTS = std::make_unique<TransportPDESolver>(conf);
    Plotting::Instance().Clear();
    conf.modelName = "STTransport";
    mainST = std::make_unique<STTransportMain>(conf);
  }

  void Run() {
    auto problem = CreateTransportProblemShared(mainTS->GetConfig().problemName);
    Solution solutionTS = mainTS->Run(problem);
    double last_mass_ts = solutionTS.values["Mass"];
    mout << "TS mass: " << last_mass_ts << endl;

    Plotting::Instance().Clear();

    problem = CreateTransportProblemShared(mainTS->GetConfig().problemName);
    Solution solutionST = mainST->Run(problem);
    mainST->PrintValues(solutionST);
    std::vector<double> massST = solutionST.multiValues["Mass"];
    double last_mass_st = massST.back();
    mout << "ST mass: " << last_mass_st << endl;

    EXPECT_NEAR(last_mass_ts, last_mass_st, TEST_TOLERANCE) << "Mass not equal.";
  }
};

TEST_P(TestTransportEquality_ST_TS, TestEquality) { Run(); }

INSTANTIATE_TEST_SUITE_P(TestTransportEquality, TestTransportEquality_ST_TS,
                         Values(std::pair{ConfigMap{}, ValueMap{{"Mass", 1.0},
                                                                {"Energy", 3.9999258206882993}}}));


#endif
