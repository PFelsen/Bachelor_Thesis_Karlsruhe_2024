#include "GMRES.hpp"
#include "STAcousticMain.hpp"
#include "SingleExecution.hpp"
#include "TestEnvironment.hpp"


#include <map>
#include <memory>
/**
 * TestParam:
 * std::string : Model
 * int : degree
 * int : level
 * std::string : Dirichlet, "Dirichlet" or ""
 * std::string : mesh
 * int : numL
 */

using TestParam = std::tuple<std::string, int, int, std::string, std::string, int>;

std::ostream &operator<<(std::ostream &os, const TestParam &param) {
  return os << "Model:  " << std::get<0>(param) << "\n"
            << "Degree:    " << std::get<1>(param) << "\n"
            << "Meshlevel: " << std::get<2>(param) << "\n"
            << "Boundary:  " << std::get<3>(param) << "\n"
            << "Mesh:      " << std::get<4>(param) << "\n"
            << "numL:      " << std::get<5>(param) << "\n";
}

std::string CreateProblemName(const std::string &basename, int degree, int numL,
                              const std::string &dirichlet) {
  std::stringstream ss;
  ss << basename;
  ss << SpaceDimension << "D";
  ss << "Degree" << degree;
  if (numL > 0) { ss << "numL" << numL; }
  ss << dirichlet;
  return ss.str();
}

class TestPolynomialSTProblem : public TestWithParam<TestParam> {
protected:
  const double TOL = 1e-12;
  std::string modelName;
  int degree;
  int level;
  std::string dirichlet;
  std::string mesh;
  int numL;
  std::string problemName;

  TestPolynomialSTProblem() :
      modelName(std::get<0>(GetParam())), degree(std::get<1>(GetParam())),
      level(std::get<2>(GetParam())), dirichlet(std::get<3>(GetParam())),
      mesh(std::get<4>(GetParam())), numL(std::get<5>(GetParam())),
      problemName(CreateProblemName("Polynomial", degree, numL, dirichlet)) {
    mout << problemName << endl;
  }

  void SetUp() override {}

  void TearDown() override { Plotting::Instance().Clear(); }

  void testSolvingAndExactSolution() {
    if (degree == 0 && modelName == "STAcousticPG") {
      mout << "degree == 0 and STPGViscoAcousticAssemble: Ok!";
      return;
    }
    PDESolverConfig mainConf;
    mainConf.mesh = mesh;
    mainConf.level = level;
    mainConf.pLevel = 0;
    mainConf.modelName = modelName;
    mainConf.problemName = problemName;
    mainConf.degree = degree;
    mainConf.timeDegree = degree;
    SingleExecution main(mainConf);

    /*STMainBuilder()
                               .WithMesh(mesh)
                               .WithPLevel(0)
                               .WithLevel(level)
                               .WithModel(assembleName)
                               .WithProblem(problemName)
                               .WithSpaceDegree(degree)
                               .WithTimeDegree(degree));*/
    main.Method();
    auto results = main.GetResults();
    double l1Error = results["L1_Error"];
    double l1intError = results["L1_int_Error"];
    double l2Error = results["L2_Error"];
    double l2intError = results["L2_int_Error"];
    double lInfError = results["LInf_Error"];
    double lInfintError = results["LInf_int_Error"];
    mout.PrintInfo(modelName, 1,
                   std::vector<PrintInfoEntry<double>>{
                       {"L1Error(solution)", l1Error},
                       {"L1Error(exact_solution)", l1intError},
                       {"L2Error(solution)", l2Error},
                       {"L2Error(exact_solution)", l2intError},
                       {"LInfError(solution)", lInfError},
                       {"LInfError(exact_solution)", lInfintError},
                   });

    EXPECT_LT(l1Error, TOL);
    EXPECT_LT(l1intError, TOL);

    EXPECT_LT(l2Error, TOL);
    EXPECT_LT(l2intError, TOL);

    EXPECT_LT(lInfError, TOL);
    EXPECT_LT(lInfintError, TOL);
  }
};

#if SpaceDimension == 1
INSTANTIATE_TEST_SUITE_P(TestPolynomialSTProblem, TestPolynomialSTProblem,
                         testing::Combine(Values("STAcousticMF", "STAcoustic"),
                                          testing::Range(1, 6),            // Degree
                                          Values(3),                       // Mesh Level
                                          Values("", "Dirichlet"),         // Boundary Type
                                          Values("SpaceTimeUnitInterval"), // Mesh
                                          Values(0)                        // NumL
                                          ));
#elif SpaceDimension == 2
INSTANTIATE_TEST_SUITE_P(TestPolynomialSTProblemWithDampingOnSquares, TestPolynomialSTProblem,
                         testing::Combine(Values("STAcousticMF", "STAcoustic", "STAcousticPG"),
                                          testing::Range(0, 3),      // Degree
                                          testing::Range(1, 3),      // Mesh Level
                                          Values("", "Dirichlet"),   // Boundary Type
                                          Values("SpaceTimeSquare"), // Mesh
                                          Values(0, 1, 2, 3)         // NumL

                                          ));

INSTANTIATE_TEST_SUITE_P(TestPolynomialSTProblemWithDampingOnTriangles, TestPolynomialSTProblem,
                         testing::Combine(Values("STAcousticMF", "STAcoustic"),
                                          testing::Range(0, 3),               // Degree
                                          testing::Range(1, 3),               // Mesh Level
                                          Values("", "Dirichlet"),            // Boundary Type
                                          Values("SpaceTimeSquareTriangles"), // Mesh
                                          Values(0, 1, 2, 3)                  // NumL

                                          ));
#elif SpaceDimension == 3
INSTANTIATE_TEST_SUITE_P(TestPolynomialSTProblem, TestPolynomialSTProblem,
                         testing::Combine(Values("STAcousticMF", "STAcoustic"),
                                          testing::Range(0, 3),        // Degree
                                          Values(0, 1),                // Mesh Level
                                          Values("", "Dirichlet"),     // Boundary Type
                                          Values("SpaceTimeUnitCube"), // Mesh
                                          Values(0)                    // NumL
                                          ));
#endif

TEST_P(TestPolynomialSTProblem, TestExactSolution) { this->testSolvingAndExactSolution(); }

int main(int argc, char **argv) {
  return MppTest(MppTestBuilder(argc, argv)
                     .WithScreenLogging()
                     .WithoutDefaultConfig()
                     .WithConfigEntry("Overlap", "STCellsWithCorners")
                     .WithConfigEntry("Distribution", "deformed_optimized")
                     .WithConfigEntry("Verbose", -1)
                     .WithConfigEntry("MainVerbose", 0)
                     .WithConfigEntry("ConfigVerbose", -1)
                     .WithConfigEntry("AssembleVerbose", 1)
                     .WithConfigEntry("MeshesVerbose", -1)
                     .WithConfigEntry("LinearVerbose", -1)
                     .WithConfigEntry("LinearPrintSteps", 1000)
                     .WithConfigEntry("MeshVerbose", -1)
                     .WithConfigEntry("ExcludedResults",
                                      "Conf, projError, exact, EE, DiscNorm, DGError, DGNorm, "
                                      "L2SpaceAtT, DGNorm_Conf, goal_functional, ||u_proj-u_h||_DG")
                     .WithConfigEntry("ResultsVerbose", 1)
                     .WithConfigEntry("LinearVerbose", 1)
                     .WithConfigEntry("LinearEpsilon", "1e-20")
                     .WithConfigEntry("LinearReduction", "1e-16")
                     .WithConfigEntry("ResultsVerbose", -1)
                     .WithPPM())
      .RUN_ALL_MPP_TESTS();
}
