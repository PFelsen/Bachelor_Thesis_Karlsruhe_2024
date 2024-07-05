#ifndef TESTELLIPTICMAIN_HPP
#define TESTELLIPTICMAIN_HPP

#include "EllipticPDESolver.hpp"
#include "TestEnvironment.hpp"

class TestEllipticPDESolver : public TestWithParam<std::pair<ConfigMap, ValueMap>> {
protected:
  std::unique_ptr<EllipticPDESolver> pdeSolver;

  double TEST_PROBLEMS_TOLERANCE = 1e-5;

  bool checkExactSolution = false;

  double TEST_TOLERANCE = 1e-8;

  ConfigMap configMap{{"PDESolverPlotting", "1"}, {"PDESolverVerbose", "2"},
                      {"AssembleVerbose", "1"},   {"MeshesVerbose", "1"},
                      {"NewtonVerbose", "1"},     {"LinearVerbose", "1"},
                      {"ConfigVerbose", "1"},     {"MeshVerbose", "2"},
                      {"MainVerbose", "1"}};

  ValueMap refValues{};

  int commSplit;

  explicit TestEllipticPDESolver(ConfigMap additionalMap1, const ValueMap &valueMap1 = {},
                                 int commSplit = 0) : commSplit(commSplit) {
    additionalMap1.merge(configMap);
    ConfigMap additionalMap2 = GetParam().first;
    additionalMap2.merge(additionalMap1);
    Config::Initialize(additionalMap2);
    if (valueMap1.empty()) refValues = GetParam().second;
    else refValues = valueMap1;
    pdeSolver = std::make_unique<EllipticPDESolver>(PDESolverConfig());
  }

  void TestRun() {
    auto problem = CreateEllipticProblemShared(pdeSolver->GetConfig().problemName);
    auto solution = pdeSolver->Run(problem, commSplit);
    EXPECT_TRUE(solution.converged);
    pdeSolver->PrintValues(solution);
    for (auto &refValue : refValues) {
      EXPECT_NEAR(refValue.second, solution.values[refValue.first], TEST_TOLERANCE)
          << refValue.first << " failed!";
    }
    EXPECT_NO_THROW(solution.vector.GetMesh().GetProcSets().CheckConsistency());
    EXPECT_NO_THROW(solution.vector.GetMatrixGraph().GetProcSets().CheckConsistency());
  }

  void TestRunProblems() {
    auto problem = CreateEllipticProblemShared(pdeSolver->GetConfig().problemName);
    auto solution = pdeSolver->Run(problem, commSplit);
    EXPECT_TRUE(solution.converged);
    pdeSolver->PrintValues(solution);
    for (auto &refValue : refValues) {
      EXPECT_NEAR(refValue.second, solution.values[refValue.first], TEST_PROBLEMS_TOLERANCE)
          << refValue.first << " failed!";
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

class TestPolynomialsLagrange : public TestEllipticPDESolver {
public:
  TestPolynomialsLagrange() :
      TestEllipticPDESolver(ConfigMap{{"Model", "Lagrange"}, {"Overlap_Distribution", "0"}},
                            ValueMap{{"L2Error", 0.0},
                                     {"MaxError", 0.0},
                                     {"FaceError", 0.0},
                                     {"FluxError", 0.0},
                                     {"EnergyError", 0.0},
                                     {"L2CellAvgError", 0.0}}) {
    checkExactSolution = true;
  }
};

class TestProblemsLagrange : public TestEllipticPDESolver {
public:
  TestProblemsLagrange() : TestEllipticPDESolver({ConfigMap{{"Model", "Lagrange"}}}) {}
};

class TestProblemsLagrangeCommSplit : public TestEllipticPDESolver {
public:
  TestProblemsLagrangeCommSplit() :
      TestEllipticPDESolver({ConfigMap{{"Model", "Lagrange"}}}, {}, 1) {}
};

class TestProblemsDG : public TestEllipticPDESolver {
public:
  TestProblemsDG() : TestEllipticPDESolver({ConfigMap{{"Model", "DG"}}}) {}
};

class TestPolynomialsDG : public TestEllipticPDESolver {
public:
  TestPolynomialsDG() :
      TestEllipticPDESolver(ConfigMap{{"Model", "DG"},
                                      {"Overlap_Distribution", "1"},
                                      {"penalty", "9"}},
                            ValueMap{{"L2Error", 0.0},
                                     {"MaxError", 0.0},
                                     {"FaceError", 0.0},
                                     {"FluxError", 0.0},
                                     {"EnergyError", 0.0},
                                     {"L2CellAvgError", 0.0}}) {
    checkExactSolution = true;
  }
};

class TestProblemsDGCommSplit : public TestEllipticPDESolver {
public:
  TestProblemsDGCommSplit() : TestEllipticPDESolver({ConfigMap{{"Model", "DG"}}}, {}, 1) {}
};

class TestPolynomialsEG : public TestEllipticPDESolver {
public:
  TestPolynomialsEG() :
      TestEllipticPDESolver(ConfigMap{{"Model", "EG"},
                                      {"Overlap_Distribution", "1"},
                                      {"penalty", "3"}},
                            ValueMap{{"L2Error", 0.0},
                                     {"MaxError", 0.0},
                                     {"FaceError", 0.0},
                                     {"FluxError", 0.0},
                                     {"EnergyError", 0.0},
                                     {"L2CellAvgError", 0.0}}) {
    checkExactSolution = true;
  }
};

class TestProblemsEG : public TestEllipticPDESolver {
public:
  TestProblemsEG() : TestEllipticPDESolver({ConfigMap{{"Model", "EG"}}}) {}
};

class TestProblemsEGCommSplit : public TestEllipticPDESolver {
public:
  TestProblemsEGCommSplit() : TestEllipticPDESolver({ConfigMap{{"Model", "EG"}}}, {}, 1) {}
};

class TestProblemsMixed : public TestEllipticPDESolver {
public:
  TestProblemsMixed() : TestEllipticPDESolver({ConfigMap{{"Model", "Mixed"}}}) {}
};

class TestProblemsMixedCommSplit : public TestEllipticPDESolver {
public:
  TestProblemsMixedCommSplit() : TestEllipticPDESolver({ConfigMap{{"Model", "Mixed"}}}, {}, 1) {}
};

class TestProblemsHybrid : public TestEllipticPDESolver {
public:
  TestProblemsHybrid() : TestEllipticPDESolver({ConfigMap{{"Model", "Hybrid"}}}) {}
};

class TestProblemsHybridCommSplit : public TestEllipticPDESolver {
public:
  TestProblemsHybridCommSplit() : TestEllipticPDESolver({ConfigMap{{"Model", "Hybrid"}}}, {}, 1) {}
};

class TestLinearSolver : public TestEllipticPDESolver {
public:
  TestLinearSolver() :
      TestEllipticPDESolver(ConfigMap{{"Model", "Lagrange"},
                                      {"Problem", "Laplace2D"},
                                      {"level", "4"},
                                      {"plevel", "4"},
                                      {"degree", "1"},
                                      {"LinearSteps", "10000"},
                                      {"LinearSolver", "LS"}},
                            ValueMap{{"L2Error", 0.0},
                                     {"MaxError", 0.0},
                                     {"FaceError", 0.0},
                                     {"FluxError", 0.0},
                                     {"EnergyError", 0.0},
                                     {"L2CellAvgError", 0.0}}) {}
};

class TestLinearSolverCommSplit : public TestEllipticPDESolver {
public:
  TestLinearSolverCommSplit() :
      TestEllipticPDESolver(ConfigMap{{"Model", "Lagrange"},
                                      {"Problem", "Laplace2D"},
                                      {"level", "4"},
                                      {"plevel", "4"},
                                      {"degree", "1"},
                                      {"LinearSteps", "10000"},
                                      {"LinearSolver", "LS"}},
                            ValueMap{{"L2Error", 0.0},
                                     {"MaxError", 0.0},
                                     {"FaceError", 0.0},
                                     {"FluxError", 0.0},
                                     {"EnergyError", 0.0},
                                     {"L2CellAvgError", 0.0}},
                            1) {}
};

class TestCG : public TestEllipticPDESolver {
public:
  TestCG() :
      TestEllipticPDESolver(ConfigMap{{"Model", "Lagrange"},
                                      {"Problem", "Laplace2D"},
                                      {"level", "4"},
                                      {"plevel", "4"},
                                      {"degree", "1"},
                                      {"LinearSteps", "400"},
                                      {"LinearSolver", "CG"}},
                            ValueMap{{"L2Error", 0.0},
                                     {"MaxError", 0.0},
                                     {"FaceError", 0.0},
                                     {"FluxError", 0.0},
                                     {"EnergyError", 0.0},
                                     {"L2CellAvgError", 0.0}}) {}
};

class TestGMRES : public TestEllipticPDESolver {
public:
  TestGMRES() :
      TestEllipticPDESolver(ConfigMap{{"Model", "Lagrange"},
                                      {"Problem", "Laplace2D"},
                                      {"level", "4"},
                                      {"plevel", "4"},
                                      {"degree", "1"},
                                      {"LinearSteps", "200"},
                                      {"LinearSolver", "GMRES"}},
                            ValueMap{{"L2Error", 0.0},
                                     {"MaxError", 0.0},
                                     {"FaceError", 0.0},
                                     {"FluxError", 0.0},
                                     {"EnergyError", 0.0},
                                     {"L2CellAvgError", 0.0}}) {}
};

#define TEST_POLYNOMIALS_2D(TestClass)                                                             \
                                                                                                   \
  INSTANTIATE_TEST_SUITE_P(TestEllipticPDESolver, TestClass,                                       \
                           Values(std::pair{ConfigMap{{"degree", "1"},                             \
                                                      {"Problem", "P0Test1D"},                     \
                                                      {"level", "8"},                              \
                                                      {"plevel", "8"}},                            \
                                            ValueMap{}},                                           \
                                  std::pair{ConfigMap{{"degree", "1"},                             \
                                                      {"Problem", "P1Test1D"},                     \
                                                      {"level", "8"},                              \
                                                      {"plevel", "8"}},                            \
                                            ValueMap{}},                                           \
                                  std::pair{ConfigMap{{"degree", "2"},                             \
                                                      {"Problem", "P2Test1D"},                     \
                                                      {"level", "8"},                              \
                                                      {"plevel", "8"}},                            \
                                            ValueMap{}},                                           \
                                  std::pair{ConfigMap{{"degree", "3"},                             \
                                                      {"Problem", "P3Test1D"},                     \
                                                      {"level", "8"},                              \
                                                      {"plevel", "8"}},                            \
                                            ValueMap{}},                                           \
                                  std::pair{ConfigMap{{"degree", "4"},                             \
                                                      {"Problem", "P4Test1D"},                     \
                                                      {"level", "8"},                              \
                                                      {"plevel", "8"}},                            \
                                            ValueMap{}},                                           \
                                  std::pair{ConfigMap{{"degree", "1"},                             \
                                                      {"Problem", "P0Test2D"},                     \
                                                      {"level", "4"},                              \
                                                      {"plevel", "4"}},                            \
                                            ValueMap{}},                                           \
                                  std::pair{ConfigMap{{"degree", "1"},                             \
                                                      {"Problem", "P1Test2D"},                     \
                                                      {"level", "4"},                              \
                                                      {"plevel", "4"}},                            \
                                            ValueMap{}},                                           \
                                  std::pair{ConfigMap{{"degree", "2"},                             \
                                                      {"Problem", "P2Test2D"},                     \
                                                      {"level", "4"},                              \
                                                      {"plevel", "4"}},                            \
                                            ValueMap{}},                                           \
                                  std::pair{ConfigMap{{"degree", "3"},                             \
                                                      {"Problem", "P3Test2D"},                     \
                                                      {"level", "4"},                              \
                                                      {"plevel", "4"}},                            \
                                            ValueMap{}},                                           \
                                  std::pair{ConfigMap{{"degree", "4"},                             \
                                                      {"Problem", "P4Test2D"},                     \
                                                      {"level", "4"},                              \
                                                      {"plevel", "4"}},                            \
                                            ValueMap{}},                                           \
                                  std::pair{ConfigMap{{"degree", "1"},                             \
                                                      {"Problem", "P0Test2DTet"},                  \
                                                      {"level", "4"},                              \
                                                      {"plevel", "4"}},                            \
                                            ValueMap{}},                                           \
                                  std::pair{ConfigMap{{"degree", "1"},                             \
                                                      {"Problem", "P1Test2DTet"},                  \
                                                      {"level", "4"},                              \
                                                      {"plevel", "4"}},                            \
                                            ValueMap{}},                                           \
                                  std::pair{ConfigMap{{"degree", "2"},                             \
                                                      {"Problem", "P2Test2DTet"},                  \
                                                      {"level", "4"},                              \
                                                      {"plevel", "4"}},                            \
                                            ValueMap{}},                                           \
                                  std::pair{ConfigMap{{"degree", "3"},                             \
                                                      {"Problem", "P3Test2DTet"},                  \
                                                      {"level", "4"},                              \
                                                      {"plevel", "4"}},                            \
                                            ValueMap{}},                                           \
                                  std::pair{ConfigMap{{"degree", "4"},                             \
                                                      {"Problem", "P4Test2DTet"},                  \
                                                      {"level", "4"},                              \
                                                      {"plevel", "4"}},                            \
                                            ValueMap{}}));                                         \
                                                                                                   \
  TEST_P(TestClass, TestRun) { TestRun(); }

#define TEST_POLYNOMIALS_3D(TestClass)                                                             \
                                                                                                   \
  INSTANTIATE_TEST_SUITE_P(TestEllipticPDESolver, TestClass,                                       \
                           Values(std::pair{ConfigMap{{"degree", "1"},                             \
                                                      {"Problem", "P0Test1D"},                     \
                                                      {"level", "8"},                              \
                                                      {"plevel", "8"}},                            \
                                            ValueMap{}},                                           \
                                  std::pair{ConfigMap{{"degree", "1"},                             \
                                                      {"Problem", "P1Test1D"},                     \
                                                      {"level", "8"},                              \
                                                      {"plevel", "8"}},                            \
                                            ValueMap{}},                                           \
                                  std::pair{ConfigMap{{"degree", "2"},                             \
                                                      {"Problem", "P2Test1D"},                     \
                                                      {"level", "8"},                              \
                                                      {"plevel", "8"}},                            \
                                            ValueMap{}},                                           \
                                  std::pair{ConfigMap{{"degree", "3"},                             \
                                                      {"Problem", "P3Test1D"},                     \
                                                      {"level", "8"},                              \
                                                      {"plevel", "8"}},                            \
                                            ValueMap{}},                                           \
                                  std::pair{ConfigMap{{"degree", "4"},                             \
                                                      {"Problem", "P4Test1D"},                     \
                                                      {"level", "8"},                              \
                                                      {"plevel", "8"}},                            \
                                            ValueMap{}},                                           \
                                  std::pair{ConfigMap{{"degree", "1"},                             \
                                                      {"Problem", "P0Test2D"},                     \
                                                      {"level", "4"},                              \
                                                      {"plevel", "4"}},                            \
                                            ValueMap{}},                                           \
                                  std::pair{ConfigMap{{"degree", "1"},                             \
                                                      {"Problem", "P1Test2D"},                     \
                                                      {"level", "4"},                              \
                                                      {"plevel", "4"}},                            \
                                            ValueMap{}},                                           \
                                  std::pair{ConfigMap{{"degree", "2"},                             \
                                                      {"Problem", "P2Test2D"},                     \
                                                      {"level", "4"},                              \
                                                      {"plevel", "4"}},                            \
                                            ValueMap{}},                                           \
                                  std::pair{ConfigMap{{"degree", "3"},                             \
                                                      {"Problem", "P3Test2D"},                     \
                                                      {"level", "4"},                              \
                                                      {"plevel", "4"}},                            \
                                            ValueMap{}},                                           \
                                  std::pair{ConfigMap{{"degree", "4"},                             \
                                                      {"Problem", "P4Test2D"},                     \
                                                      {"level", "4"},                              \
                                                      {"plevel", "4"}},                            \
                                            ValueMap{}},                                           \
                                  std::pair{ConfigMap{{"degree", "1"},                             \
                                                      {"Problem", "P0Test3D"},                     \
                                                      {"level", "2"},                              \
                                                      {"plevel", "2"}},                            \
                                            ValueMap{}},                                           \
                                  std::pair{ConfigMap{{"degree", "1"},                             \
                                                      {"Problem", "P1Test3D"},                     \
                                                      {"level", "2"},                              \
                                                      {"plevel", "2"}},                            \
                                            ValueMap{}},                                           \
                                  std::pair{ConfigMap{{"degree", "2"},                             \
                                                      {"Problem", "P2Test3D"},                     \
                                                      {"level", "2"},                              \
                                                      {"plevel", "2"}},                            \
                                            ValueMap{}},                                           \
                                  std::pair{ConfigMap{{"degree", "3"},                             \
                                                      {"Problem", "P3Test3D"},                     \
                                                      {"level", "1"},                              \
                                                      {"plevel", "1"}},                            \
                                            ValueMap{}},                                           \
                                  std::pair{ConfigMap{{"degree", "4"},                             \
                                                      {"Problem", "P4Test3D"},                     \
                                                      {"level", "1"},                              \
                                                      {"plevel", "1"}},                            \
                                            ValueMap{}},                                           \
                                  std::pair{ConfigMap{{"degree", "1"},                             \
                                                      {"Problem", "P0Test2DTet"},                  \
                                                      {"level", "4"},                              \
                                                      {"plevel", "4"}},                            \
                                            ValueMap{}},                                           \
                                  std::pair{ConfigMap{{"degree", "1"},                             \
                                                      {"Problem", "P1Test2DTet"},                  \
                                                      {"level", "4"},                              \
                                                      {"plevel", "4"}},                            \
                                            ValueMap{}},                                           \
                                  std::pair{ConfigMap{{"degree", "2"},                             \
                                                      {"Problem", "P2Test2DTet"},                  \
                                                      {"level", "4"},                              \
                                                      {"plevel", "4"}},                            \
                                            ValueMap{}},                                           \
                                  std::pair{ConfigMap{{"degree", "3"},                             \
                                                      {"Problem", "P3Test2DTet"},                  \
                                                      {"level", "4"},                              \
                                                      {"plevel", "4"}},                            \
                                            ValueMap{}},                                           \
                                  std::pair{ConfigMap{{"degree", "4"},                             \
                                                      {"Problem", "P4Test2DTet"},                  \
                                                      {"level", "4"},                              \
                                                      {"plevel", "4"}},                            \
                                            ValueMap{}},                                           \
                                  std::pair{ConfigMap{{"degree", "1"},                             \
                                                      {"Problem", "P0Test3DTet"},                  \
                                                      {"level", "2"},                              \
                                                      {"plevel", "2"}},                            \
                                            ValueMap{}},                                           \
                                  std::pair{ConfigMap{{"degree", "1"},                             \
                                                      {"Problem", "P1Test3DTet"},                  \
                                                      {"level", "2"},                              \
                                                      {"plevel", "2"}},                            \
                                            ValueMap{}},                                           \
                                  std::pair{ConfigMap{{"degree", "2"},                             \
                                                      {"Problem", "P2Test3DTet"},                  \
                                                      {"level", "2"},                              \
                                                      {"plevel", "2"}},                            \
                                            ValueMap{}},                                           \
                                  std::pair{ConfigMap{{"degree", "3"},                             \
                                                      {"Problem", "P3Test3DTet"},                  \
                                                      {"level", "1"},                              \
                                                      {"plevel", "1"}},                            \
                                            ValueMap{}},                                           \
                                  std::pair{ConfigMap{{"degree", "4"},                             \
                                                      {"Problem", "P4Test3DTet"},                  \
                                                      {"level", "1"},                              \
                                                      {"plevel", "1"}},                            \
                                            ValueMap{}}));                                         \
                                                                                                   \
  TEST_P(TestClass, TestRun) { TestRun(); }

#define TEST_PROBLEMS_LAGRANGE(TestClass)                                                          \
                                                                                                   \
  INSTANTIATE_TEST_SUITE_P(TestEllipticPDESolver, TestClass,                                       \
                           Values(std::pair{ConfigMap{{"level", "8"},                              \
                                                      {"plevel", "8"},                             \
                                                      {"degree", "1"},                             \
                                                      {"Problem", "Laplace1D"}},                   \
                                            ValueMap{{"L2Error", 0.0},                             \
                                                     {"MaxError", 0.0},                            \
                                                     {"FaceError", 0.0},                           \
                                                     {"FluxError", 0.0},                           \
                                                     {"EnergyError", 0.0},                         \
                                                     {"L2CellAvgError", 0.0},                      \
                                                     {"Inflow", -1.0},                             \
                                                     {"Outflow", 1.0}}},                           \
                                                                                                   \
                                  std::pair{ConfigMap{{"level", "4"},                              \
                                                      {"plevel", "4"},                             \
                                                      {"degree", "1"},                             \
                                                      {"Problem", "Laplace2D"}},                   \
                                            ValueMap{{"L2Error", 0.0},                             \
                                                     {"MaxError", 0.0},                            \
                                                     {"FaceError", 0.0},                           \
                                                     {"FluxError", 0.0},                           \
                                                     {"EnergyError", 0.0},                         \
                                                     {"L2CellAvgError", 0.0}}},                    \
                                                                                                   \
                                  std::pair{ConfigMap{{"level", "8"},                              \
                                                      {"plevel", "8"},                             \
                                                      {"degree", "1"},                             \
                                                      {"Problem", "Discontinuous1D"}},             \
                                            ValueMap{{"FluxError", 0.0}, {"FluxLoss", 0.0}}},      \
                                                                                                   \
                                  std::pair{ConfigMap{{"level", "4"},                              \
                                                      {"plevel", "4"},                             \
                                                      {"degree", "1"},                             \
                                                      {"Problem", "Discontinuous2D"}},             \
                                            ValueMap{{"H1", 0.98441779}, {"L2", 0.45283722}}},     \
                                                                                                   \
                                  std::pair{ConfigMap{{"level", "0"},                              \
                                                      {"plevel", "0"},                             \
                                                      {"degree", "1"},                             \
                                                      {"Problem", "LaplaceSquare500"}},            \
                                            ValueMap{{"H1", 1.2170148}, {"L2", 0.6145404}}},       \
                                                                                                   \
                                  std::pair{ConfigMap{{"level", "4"},                              \
                                                      {"plevel", "4"},                             \
                                                      {"degree", "1"},                             \
                                                      {"Problem", "Divergent"}},                   \
                                            ValueMap{{"H1", 8.1237594e-05},                        \
                                                     {"L2", 3.184605e-05}}},                       \
                                                                                                   \
                                  std::pair{ConfigMap{{"level", "4"},                              \
                                                      {"plevel", "4"},                             \
                                                      {"degree", "1"},                             \
                                                      {"Problem", "Kellogg"}},                     \
                                            ValueMap{{"H1", 0.26443412}, {"L2", 0.11422282}}},     \
                                                                                                   \
                                  std::pair{ConfigMap{{"level", "4"},                              \
                                                      {"plevel", "4"},                             \
                                                      {"degree", "1"},                             \
                                                      {"Problem", "Rock"}},                        \
                                            ValueMap{{"OutflowLeft", 1.5},                         \
                                                     {"OutflowRight", 1.5}}}));                    \
                                                                                                   \
  TEST_P(TestClass, TestRunProblems) { TestRunProblems(); }

#define TEST_PROBLEMS_DG(TestClass)                                                                \
                                                                                                   \
  INSTANTIATE_TEST_SUITE_P(TestEllipticPDESolver, TestClass,                                       \
                           Values(std::pair{ConfigMap{{"level", "8"},                              \
                                                      {"plevel", "8"},                             \
                                                      {"degree", "1"},                             \
                                                      {"Problem", "Laplace1D"}},                   \
                                            ValueMap{{"L2Error", 0.0},                             \
                                                     {"MaxError", 0.0},                            \
                                                     {"FaceError", 0.0},                           \
                                                     {"FluxError", 0.0},                           \
                                                     {"EnergyError", 0.0},                         \
                                                     {"L2CellAvgError", 0.0}}},                    \
                                                                                                   \
                                  std::pair{ConfigMap{{"level", "4"},                              \
                                                      {"plevel", "4"},                             \
                                                      {"degree", "1"},                             \
                                                      {"Problem", "Laplace2D"}},                   \
                                            ValueMap{{"L2Error", 0.0},                             \
                                                     {"MaxError", 0.0},                            \
                                                     {"FaceError", 0.0},                           \
                                                     {"FluxError", 0.0},                           \
                                                     {"EnergyError", 0.0},                         \
                                                     {"L2CellAvgError", 0.0}}},                    \
                                                                                                   \
                                  std::pair{ConfigMap{{"level", "8"},                              \
                                                      {"plevel", "8"},                             \
                                                      {"degree", "1"},                             \
                                                      {"Problem", "Discontinuous1D"}},             \
                                            ValueMap{{"FluxError", 0.0}, {"FluxLoss", 0.0}}},      \
                                                                                                   \
                                  std::pair{ConfigMap{{"level", "4"},                              \
                                                      {"plevel", "4"},                             \
                                                      {"degree", "1"},                             \
                                                      {"Problem", "Discontinuous2D"}},             \
                                            ValueMap{{"H1", 0.9859994}, {"L2", 0.45385792}}},      \
                                                                                                   \
                                  std::pair{ConfigMap{{"level", "0"},                              \
                                                      {"plevel", "0"},                             \
                                                      {"degree", "1"},                             \
                                                      {"Problem", "LaplaceSquare500"}},            \
                                            ValueMap{{"H1", 1.2206765}, {"L2", 0.62824253}}},      \
                                                                                                   \
                                  std::pair{ConfigMap{{"level", "4"},                              \
                                                      {"plevel", "4"},                             \
                                                      {"degree", "1"},                             \
                                                      {"Problem", "Divergent"}},                   \
                                            ValueMap{{"H1", 0.019451009}, {"L2", 0.0051816658}}},  \
                                                                                                   \
                                  std::pair{ConfigMap{{"level", "4"},                              \
                                                      {"plevel", "4"},                             \
                                                      {"degree", "1"},                             \
                                                      {"Problem", "Kellogg"}},                     \
                                            ValueMap{{"H1", 0.78891974}, {"L2", 0.12141242}}},     \
                                                                                                   \
                                  std::pair{ConfigMap{{"level", "4"},                              \
                                                      {"plevel", "4"},                             \
                                                      {"degree", "1"},                             \
                                                      {"Problem", "Rock"}},                        \
                                            ValueMap{{"OutflowLeft", 1.5},                         \
                                                     {"OutflowRight", 1.5}}}));                    \
                                                                                                   \
  TEST_P(TestClass, TestRunProblems) { TestRunProblems(); }

#define TEST_PROBLEMS_EG(TestClass)                                                                \
                                                                                                   \
  INSTANTIATE_TEST_SUITE_P(TestEllipticPDESolver, TestClass,                                       \
                           Values(std::pair{ConfigMap{{"level", "8"},                              \
                                                      {"plevel", "8"},                             \
                                                      {"degree", "1"},                             \
                                                      {"Problem", "Laplace1D"}},                   \
                                            ValueMap{{"L2Error", 0.0},                             \
                                                     {"MaxError", 0.0},                            \
                                                     {"FaceError", 0.0},                           \
                                                     {"FluxError", 0.0},                           \
                                                     {"EnergyError", 0.0},                         \
                                                     {"L2CellAvgError", 0.0},                      \
                                                     {"Inflow", -1.0},                             \
                                                     {"Outflow", 1.0}}},                           \
                                                                                                   \
                                  std::pair{ConfigMap{{"level", "4"},                              \
                                                      {"plevel", "4"},                             \
                                                      {"degree", "1"},                             \
                                                      {"Problem", "Laplace2D"}},                   \
                                            ValueMap{{"L2Error", 0.0},                             \
                                                     {"MaxError", 0.0},                            \
                                                     {"FaceError", 0.0},                           \
                                                     {"FluxError", 0.0},                           \
                                                     {"EnergyError", 0.0},                         \
                                                     {"L2CellAvgError", 0.0}}},                    \
                                                                                                   \
                                  std::pair{ConfigMap{{"level", "8"},                              \
                                                      {"plevel", "8"},                             \
                                                      {"degree", "1"},                             \
                                                      {"Problem", "Discontinuous1D"}},             \
                                            ValueMap{{"FluxError", 0.0}, {"FluxLoss", 0.0}}},      \
                                                                                                   \
                                  std::pair{ConfigMap{{"level", "4"},                              \
                                                      {"plevel", "4"},                             \
                                                      {"degree", "1"},                             \
                                                      {"Problem", "Discontinuous2D"},              \
                                                      {"Overlap_Distribution", "1"},               \
                                                      {"sign", "-1"},                              \
                                                      {"penalty", "40"}},                          \
                                            ValueMap{{"H1", 0.98418176}, {"L2", 0.45288666}}},     \
                                                                                                   \
                                  std::pair{ConfigMap{{"level", "0"},                              \
                                                      {"plevel", "0"},                             \
                                                      {"degree", "1"},                             \
                                                      {"Problem", "LaplaceSquare500"}},            \
                                            ValueMap{{"H1", 1.2133397}, {"L2", 0.61840207}}},      \
                                                                                                   \
                                  std::pair{ConfigMap{{"level", "4"},                              \
                                                      {"plevel", "4"},                             \
                                                      {"degree", "1"},                             \
                                                      {"Problem", "Divergent"},                    \
                                                      {"Overlap_Distribution", "1"},               \
                                                      {"sign", "-1"},                              \
                                                      {"penalty", "40"}},                          \
                                            ValueMap{{"H1", 0.00010488862},                        \
                                                     {"L2", 7.0864056e-05}}},                      \
                                                                                                   \
                                  std::pair{ConfigMap{{"level", "4"},                              \
                                                      {"plevel", "4"},                             \
                                                      {"degree", "1"},                             \
                                                      {"Problem", "Kellogg"}},                     \
                                            ValueMap{{"H1", 0.26332864}, {"L2", 0.11452042}}},     \
                                                                                                   \
                                  std::pair{ConfigMap{{"level", "4"},                              \
                                                      {"plevel", "4"},                             \
                                                      {"degree", "1"},                             \
                                                      {"Problem", "Rock"},                         \
                                                      {"Overlap_Distribution", "1"},               \
                                                      {"sign", "-1"},                              \
                                                      {"penalty", "40"}},                          \
                                            ValueMap{{"OutflowLeft", 1.5},                         \
                                                     {"OutflowRight", 1.5}}}));                    \
                                                                                                   \
  TEST_P(TestClass, TestRunProblems) { TestRunProblems(); }

#define TEST_PROBLEMS_MIXED(TestClass)                                                             \
                                                                                                   \
  INSTANTIATE_TEST_SUITE_P(TestEllipticPDESolver, TestClass,                                       \
                           Values(std::pair{ConfigMap{{"level", "8"},                              \
                                                      {"plevel", "8"},                             \
                                                      {"DualPrimal", "1"},                         \
                                                      {"degree", "1"},                             \
                                                      {                                            \
                                                          "Problem",                               \
                                                          "Laplace1D",                             \
                                                      }},                                          \
                                            ValueMap{{"FluxError", 0.0},                           \
                                                     {"EnergyError", 0.0},                         \
                                                     {"L2CellAvgError", 0.0},                      \
                                                     {"Inflow", -1.0},                             \
                                                     {"Outflow", 1.0},                             \
                                                     {"DualPrimal", 0.0}}},                        \
                                                                                                   \
                                  std::pair{ConfigMap{{"level", "8"},                              \
                                                      {"plevel", "8"},                             \
                                                      {"DualPrimal", "1"},                         \
                                                      {"degree", "1"},                             \
                                                      {"Problem", "P0Test1D"}},                    \
                                            ValueMap{{"L2Error", 0.0},                             \
                                                     {"MaxError", 0.0},                            \
                                                     {"FaceError", 0.0},                           \
                                                     {"FluxError", 0.0},                           \
                                                     {"EnergyError", 0.0},                         \
                                                     {"L2CellAvgError", 0.0},                      \
                                                     {"DualPrimal", 0.0}}},                        \
                                                                                                   \
                                  std::pair{ConfigMap{{"level", "4"},                              \
                                                      {"plevel", "4"},                             \
                                                      {"DualPrimal", "1"},                         \
                                                      {"degree", "1"},                             \
                                                      {"Problem", "Laplace2D"}},                   \
                                            ValueMap{{"FluxError", 0.0},                           \
                                                     {"EnergyError", 0.0},                         \
                                                     {"L2CellAvgError", 0.0},                      \
                                                     {"DualPrimal", 0.0}}},                        \
                                                                                                   \
                                  std::pair{ConfigMap{{"level", "4"},                              \
                                                      {"plevel", "4"},                             \
                                                      {"DualPrimal", "1"},                         \
                                                      {"degree", "2"},                             \
                                                      {"Problem", "Laplace2D"}},                   \
                                            ValueMap{{"FluxError", 0.0},                           \
                                                     {"EnergyError", 0.0},                         \
                                                     {"L2CellAvgError", 0.0},                      \
                                                     {"DualPrimal", 0.0}}},                        \
                                                                                                   \
                                  std::pair{ConfigMap{{"level", "4"},                              \
                                                      {"plevel", "4"},                             \
                                                      {"DualPrimal", "1"},                         \
                                                      {"degree", "3"},                             \
                                                      {"Problem", "Laplace2D"}},                   \
                                            ValueMap{{"FluxError", 0.0},                           \
                                                     {"EnergyError", 0.0},                         \
                                                     {"L2CellAvgError", 0.0},                      \
                                                     {"DualPrimal", 0.0}}},                        \
                                                                                                   \
                                  std::pair{ConfigMap{{"level", "4"},                              \
                                                      {"plevel", "4"},                             \
                                                      {"DualPrimal", "1"},                         \
                                                      {"degree", "4"},                             \
                                                      {"Problem", "Laplace2D"}},                   \
                                            ValueMap{{"FluxError", 0.0},                           \
                                                     {"EnergyError", 0.0},                         \
                                                     {"L2CellAvgError", 0.0},                      \
                                                     {"DualPrimal", 0.0}}},                        \
                                                                                                   \
                                  std::pair{ConfigMap{{"level", "4"},                              \
                                                      {"plevel", "4"},                             \
                                                      {"DualPrimal", "1"},                         \
                                                      {"degree", "1"},                             \
                                                      {"Problem", "P0Test2D"}},                    \
                                            ValueMap{{"L2Error", 0.0},                             \
                                                     {"MaxError", 0.0},                            \
                                                     {"FaceError", 0.0},                           \
                                                     {"FluxError", 0.0},                           \
                                                     {"EnergyError", 0.0},                         \
                                                     {"L2CellAvgError", 0.0},                      \
                                                     {"DualPrimal", 0.0}}},                        \
                                                                                                   \
                                  std::pair{ConfigMap{{"level", "8"},                              \
                                                      {"plevel", "8"},                             \
                                                      {"DualPrimal", "1"},                         \
                                                      {"degree", "1"},                             \
                                                      {"Problem", "Discontinuous1D"}},             \
                                            ValueMap{{"FluxError", 0.0},                           \
                                                     {"FluxLoss", 0.0},                            \
                                                     {"DualPrimal", 0.0}}},                        \
                                                                                                   \
                                  std::pair{ConfigMap{{"level", "4"},                              \
                                                      {"plevel", "4"},                             \
                                                      {"DualPrimal", "1"},                         \
                                                      {"degree", "1"},                             \
                                                      {"Problem", "Discontinuous2D"}},             \
                                            ValueMap{{"H1", 1.0056944},                            \
                                                     {"L2", 0.46567983},                           \
                                                     {"DualPrimal", 0.14534581}}},                 \
                                                                                                   \
                                  std::pair{ConfigMap{{"level", "0"},                              \
                                                      {"plevel", "0"},                             \
                                                      {"DualPrimal", "1"},                         \
                                                      {"degree", "1"},                             \
                                                      {"Problem", "LaplaceSquare500"}},            \
                                            ValueMap{{"H1", 1.3096751},                            \
                                                     {"L2", 0.6884561},                            \
                                                     {"DualPrimal", 0.37123006}}},                 \
                                                                                                   \
                                  std::pair{ConfigMap{{"level", "4"},                              \
                                                      {"plevel", "4"},                             \
                                                      {"DualPrimal", "1"},                         \
                                                      {"degree", "1"},                             \
                                                      {"Problem", "Divergent"}},                   \
                                            ValueMap{{"H1", 0.3086837},                            \
                                                     {"L2", 0.069566362},                          \
                                                     {"DualPrimal", 0.35108471}}},                 \
                                                                                                   \
                                  std::pair{ConfigMap{{"level", "4"},                              \
                                                      {"plevel", "4"},                             \
                                                      {"DualPrimal", "1"},                         \
                                                      {"degree", "1"},                             \
                                                      {"Problem", "Kellogg"}},                     \
                                            ValueMap{{"H1", 0.38900275},                           \
                                                     {"L2", 0.12353347},                           \
                                                     {"DualPrimal", 0.73041508}}},                 \
                                                                                                   \
                                  std::pair{ConfigMap{{"level", "4"},                              \
                                                      {"plevel", "4"},                             \
                                                      {"DualPrimal", "1"},                         \
                                                      {"degree", "1"},                             \
                                                      {"Problem", "Rock"}},                        \
                                            ValueMap{{"OutflowLeft", 1.5},                         \
                                                     {"OutflowRight", 1.5},                        \
                                                     {"DualPrimal", 0.12389593}}}));               \
                                                                                                   \
  TEST_P(TestClass, TestRunProblems) { TestRunProblems(); }

#define TEST_SOLVERS(TestClass)                                                                    \
                                                                                                   \
  INSTANTIATE_TEST_SUITE_P(TestEllipticPDESolver, TestClass,                                       \
                           Values(std::pair{ConfigMap{{"Preconditioner", "Jacobi"}}, ValueMap{}},  \
                                  std::pair{ConfigMap{{"Preconditioner", "SSOR"}}, ValueMap{}},    \
                                  std::pair{ConfigMap{{"Preconditioner", "SuperLU"}}, ValueMap{}}, \
                                  std::pair{ConfigMap{{"Preconditioner", "PS"}}, ValueMap{}},      \
                                  std::pair{ConfigMap{{"Preconditioner", "Multigrid"}},            \
                                            ValueMap{}}));                                         \
                                                                                                   \
  TEST_P(TestClass, TestRunProblems) { TestRunProblems(); }

#endif // TESTELLIPTICMAIN_HPP
