//
// Created by lstengel on 25.07.22.
//

#ifndef TESTVECTORVALUEDMAIN_HPP
#define TESTVECTORVALUEDMAIN_HPP

#include "SingleExecution.hpp"
#include "TestEnvironment.hpp"
#include "VectorValuedMain.hpp"

typedef std::map<std::string, std::string> ConfigMap;

struct TestParameter {
  std::map<std::string, std::string> ConfigMap;
  int maxLevel;
  int degree;
};

class TestVectorValuedMainRates : public TestWithParam<TestParameter> {
protected:
  std::unique_ptr<VectorValuedMain> main;

  ConfigMap configMap{{"PDESolverVerbose", "2"}, {"AssembleVerbose", "1"}, {"MeshesVerbose", "1"},
                      {"NewtonVerbose", "1"},    {"LinearVerbose", "1"},   {"ConfigVerbose", "5"},
                      {"MeshVerbose", "2"},      {"MainVerbose", "1"}};

  ConfigMap additionalMap1{};

  TestVectorValuedMainRates() = default;
};

class TestVectorValuedMain : public TestWithParam<std::pair<ConfigMap, ValueMap>> {
protected:
  std::unique_ptr<VectorValuedMain> main;

  double TEST_TOLERANCE = 1e-9;

  ConfigMap configMap{{"PDESolverVerbose", "2"},    {"AssembleVerbose", "1"},
                      {"MeshesVerbose", "1"},       {"NewtonVerbose", "1"},
                      {"LinearVerbose", "1"},       {"ConfigVerbose", "5"},
                      {"MeshVerbose", "2"},         {"MainVerbose", "1"},
                      {"LinearReduction", "1e-15"}, {"LinearEpsilon", "1e-15"}};

  ValueMap refValues{};

  explicit TestVectorValuedMain(ConfigMap additionalMap1, ValueMap valueMap1 = {}) {
    additionalMap1.merge(configMap);
    ConfigMap additionalMap2 = GetParam().first;
    additionalMap2.merge(additionalMap1);
    Config::Initialize(additionalMap2);
    if (valueMap1.empty()) refValues = GetParam().second;
    else refValues = valueMap1;
    main = std::make_unique<VectorValuedMain>(PDESolverConfig());
  }

  void TestRun() {
    auto problem = CreateVectorValuedProblemShared(main->GetConfig().problemName);
    Solution solution = main->Run(problem);
    ValueMap values = solution.values;
    for (auto &refValue : refValues) {
      EXPECT_NEAR(refValue.second, values[refValue.first], TEST_TOLERANCE);
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

TEST_P(TestVectorValuedMainRates, convergenceRate) {
  std::vector<double> errL2(GetParam().maxLevel);
  std::vector<double> errEnergy(GetParam().maxLevel);

  for (int l = 0; l < GetParam().maxLevel; ++l) {
    additionalMap1.merge(configMap);
    ConfigMap params = {};
    params = {{"level", std::to_string(l)}, {"degree", std::to_string(GetParam().degree)}};
    additionalMap1.merge(params);
    ConfigMap additionalMap2 = GetParam().ConfigMap;
    additionalMap2.merge(additionalMap1);
    Config::Initialize(additionalMap2);
    main = std::make_unique<VectorValuedMain>(PDESolverConfig());
    auto problem = CreateVectorValuedProblemShared(main->GetConfig().problemName);
    Solution s = main->Run(problem);
    errL2[l] = main->EvaluateQuantity(s.vector, "L2");
    errEnergy[l] = main->EvaluateQuantity(s.vector, "Energy");

    Config::Close();
  }
  for (int i = 0; i < errL2.size() - 1; ++i) {
    EXPECT_LE(errL2[i + 1], errL2[i]);
    EXPECT_LE(errEnergy[i + 1], errEnergy[i]);
  }
  for (int i = 1; i < errL2.size(); ++i) {
    mout << "errL2[" << i - 1 << "] " << errL2[i - 1] << endl << endl;
    mout << "errL2[" << i << "] " << errL2[i] << endl << endl;
    mout << "errL2[" << i - 1 << "]/errL2[" << i << "] " << errL2[i - 1] / errL2[i] << endl << endl;
    mout << "EOC L2 " << log2(errL2[i - 1] / errL2[i]) << endl << endl;
    mout << "EOC Energy " << log2(errEnergy[i - 1] / errEnergy[i]) << endl << endl;
    if (GetParam().degree == 1) {
      // EOC
      EXPECT_NEAR(log2(errL2[i - 1] / errL2[i]), 2.0, 0.32);
      EXPECT_NEAR(log2(errEnergy[i - 1] / errEnergy[i]), 1.0, 0.2);
      // convergence rate
      /*EXPECT_NEAR(errL2[i - 1] / errL2[i], 4.0, 0.3);
      EXPECT_NEAR(errEnergy[i - 1] / errEnergy[i], 2.0, 0.3);*/
    } else if (GetParam().degree == 2) {
      // EOC
      EXPECT_NEAR(log2(errL2[i - 1] / errL2[i]), 3.0, 0.2);
      EXPECT_NEAR(log2(errEnergy[i - 1] / errEnergy[i]), 2.0, 0.2);
      // convergence rate
      /*EXPECT_NEAR(errL2[i - 1] / errL2[i], 8.0, 0.3);
      EXPECT_NEAR(errEnergy[i - 1] / errEnergy[i], 4.0, 0.3);*/
    } else if (GetParam().degree == 3) {
      // EOC
      EXPECT_NEAR(log2(errL2[i - 1] / errL2[i]), 4.0, 0.2);
      EXPECT_NEAR(log2(errEnergy[i - 1] / errEnergy[i]), 3.0, 0.2);
      // convergence rate
      /*EXPECT_NEAR(errL2[i - 1] / errL2[i], 16.0, 0.3);
      EXPECT_NEAR(errEnergy[i - 1] / errEnergy[i], 8.0, 0.3);*/
    }
  }
}

class TestVectorPolynomialsLagrange : public TestVectorValuedMain {
public:
  TestVectorPolynomialsLagrange() :
      TestVectorValuedMain(ConfigMap{{"Model", "VectorValuedLagrange"}},
                           ValueMap{{"L2Error", 0.0},
                                    {"MaxError", 0.0},
                                    {"FaceError", 0.0},
                                    {"L2CellAvgError", 0.0}}) {}
};

class TestVectorPolynomialsDG : public TestVectorValuedMain {
public:
  TestVectorPolynomialsDG() :
      TestVectorValuedMain(ConfigMap{{"Model", "DGVectorValuedAssemble"},
                                     {"sign", "-1"},
                                     {"Overlap_Distribution", "1"},
                                     {"penalty", "2"}},
                           ValueMap{{"L2Error", 0.0},
                                    {"MaxError", 0.0},
                                    {"FaceError", 0.0},
                                    {"L2CellAvgError", 0.0}}) {}
};

class TestVectorPolynomialsEG : public TestVectorValuedMain {
public:
  TestVectorPolynomialsEG() :
      TestVectorValuedMain(ConfigMap{{"Model", "EGVectorValuedAssemble"},
                                     {"sign", "-1"},
                                     {"Overlap_Distribution", "1"},
                                     {"penalty", "10"}},
                           ValueMap{{"L2Error", 0.0},
                                    {"MaxError", 0.0},
                                    {"FaceError", 0.0},
                                    {"L2CellAvgError", 0.0}}) {}
};

class TestProblemRatesLagrange : public TestVectorValuedMainRates {
public:
  TestProblemRatesLagrange() : TestVectorValuedMainRates() {}
};

std::vector<TestParameter> Parameters() {
  return std::vector<TestParameter>{
#if SpaceDimension >= 2
      TestParameter(
          {ConfigMap{{"Model", "VectorValuedLagrange"}, {"plevel", "0"}, {"Problem", "P2Test2D"}},
           4, 1}),
      TestParameter(
          {ConfigMap{{"Model", "VectorValuedLagrange"}, {"plevel", "0"}, {"Problem", "P3Test2D"}},
           4, 2}),
      TestParameter(
          {ConfigMap{{"Model", "VectorValuedLagrange"}, {"plevel", "0"}, {"Problem", "P4Test2D"}},
           4, 3}),
      TestParameter({ConfigMap{{"Model", "DGVectorValuedAssemble"},
                               {"plevel", "0"},
                               {"Problem", "P2Test2D"},
                               {"sign", "-1"},
                               {"penalty", "8.0"}},
                     4, 1}),
      TestParameter({ConfigMap{{"Model", "DGVectorValuedAssemble"},
                               {"plevel", "0"},
                               {"Problem", "P3Test2D"},
                               {"sign", "-1"},
                               {"penalty", "14.0"}},
                     4, 2}),
      TestParameter({ConfigMap{{"Model", "DGVectorValuedAssemble"},
                               {"plevel", "0"},
                               {"Problem", "P4Test2D"},
                               {"sign", "-1"},
                               {"penalty", "3.0"}},
                     4, 3}),
      TestParameter({ConfigMap{{"Model", "EGVectorValuedAssemble"},
                               {"plevel", "0"},
                               {"Problem", "P2Test2D"},
                               {"Overlap_Distribution", "1"},
                               {"sign", "-1"},
                               {"penalty", "4.0"}},
                     4, 1}),
      TestParameter({ConfigMap{{"Model", "EGVectorValuedAssemble"},
                               {"plevel", "0"},
                               {"Problem", "P3Test2D"},
                               {"Overlap_Distribution", "1"},
                               {"sign", "-1"},
                               {"penalty", "4.0"}},
                     4, 2}),
      TestParameter({ConfigMap{{"Model", "EGVectorValuedAssemble"},
                               {"plevel", "0"},
                               {"Problem", "P4Test2D"},
                               {"Overlap_Distribution", "1"},
                               {"sign", "-1"},
                               {"penalty", "4.0"}},
                     4, 3}),
#endif
#if SpaceDimension >= 3
      TestParameter(
          {ConfigMap{{"Model", "VectorValuedLagrange"}, {"plevel", "0"}, {"Problem", "P2Test3D"}},
           4, 1}),
      TestParameter(
          {ConfigMap{{"Model", "VectorValuedLagrange"}, {"plevel", "0"}, {"Problem", "P3Test3D"}},
           3, 2}),
      TestParameter(
          {ConfigMap{{"Model", "VectorValuedLagrange"}, {"plevel", "0"}, {"Problem", "P4Test3D"}},
           3, 2}),
      TestParameter({ConfigMap{{"Model", "DGVectorValuedAssemble"},
                               {"plevel", "0"},
                               {"Problem", "P2Test3D"},
                               {"Overlap_Distribution", "1"},
                               {"sign", "-1"},
                               {"penalty", "8.0"}},
                     4, 1}),
      TestParameter({ConfigMap{{"Model", "DGVectorValuedAssemble"},
                               {"plevel", "0"},
                               {"Problem", "P3Test3D"},
                               {"Overlap_Distribution", "1"},
                               {"sign", "-1"},
                               {"penalty", "30.0"}},
                     3, 2}),
      TestParameter({ConfigMap{{"Model", "DGVectorValuedAssemble"},
                               {"plevel", "0"},
                               {"Problem", "P4Test3D"},
                               {"Overlap_Distribution", "1"},
                               {"sign", "-1"},
                               {"penalty", "4.0"}},
                     3, 3}),
      TestParameter({ConfigMap{{"Model", "EGVectorValuedAssemble"},
                               {"plevel", "0"},
                               {"Problem", "P2Test3D"},
                               {"Overlap_Distribution", "1"},
                               {"sign", "-1"},
                               {"penalty", "8.0"}},
                     3, 1}),
      TestParameter({ConfigMap{{"Model", "EGVectorValuedAssemble"},
                               {"plevel", "0"},
                               {"Problem", "P3Test3D"},
                               {"Overlap_Distribution", "1"},
                               {"sign", "-1"},
                               {"penalty", "3.0"}},
                     3, 2}),
      TestParameter({ConfigMap{{"Model", "EGVectorValuedAssemble"},
                               {"plevel", "0"},
                               {"Problem", "P4Test3D"},
                               {"Overlap_Distribution", "1"},
                               {"sign", "-1"},
                               {"penalty", "5.0"}},
                     3, 3}),
#endif
  };
}

#define TEST_VECTOR_POLYNOMIALS1D(TestClass)                                                       \
  INSTANTIATE_TEST_SUITE_P(TestVectorValuedMain##1D, TestClass,                                    \
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
                                                      {"level", "7"},                              \
                                                      {"plevel", "7"}},                            \
                                            ValueMap{}},                                           \
                                  std::pair{ConfigMap{{"degree", "4"},                             \
                                                      {"Problem", "P4Test1D"},                     \
                                                      {"level", "7"},                              \
                                                      {"plevel", "7"}},                            \
                                            ValueMap{}}));


#define TEST_VECTOR_POLYNOMIALS2D(TestClass)                                                       \
  INSTANTIATE_TEST_SUITE_P(TestVectorValuedMain##2D, TestClass,                                    \
                           Values(std::pair{ConfigMap{{"degree", "1"},                             \
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
                                            ValueMap{}}));

#define TEST_VECTOR_POLYNOMIALS3D(TestClass)                                                       \
  INSTANTIATE_TEST_SUITE_P(TestVectorValuedMain##3D, TestClass,                                    \
                           Values(std::pair{ConfigMap{{"degree", "1"},                             \
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
                                            ValueMap{}}));

#define INITIALIZE_TEST_P(TestClass)                                                               \
                                                                                                   \
  TEST_P(TestClass, TestRun) { TestRun(); }


#endif // TESTVECTORVALUEDMAIN_HPP
