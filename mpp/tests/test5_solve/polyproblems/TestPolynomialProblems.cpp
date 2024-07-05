#include "TestPolynomialProblems.hpp"
#include "Newton.hpp"

#include "ArgyrisDiscretization.hpp"
#include "LagrangeDiscretization.hpp"
#include "MeshesCreator.hpp"
#include "TestEnvironment.hpp"

class TestMonomialProblems : public TestWithParam<int> {
protected:
  int degree;
  CELLTYPE cellType;

  std::unique_ptr<Meshes> meshes;
  std::shared_ptr<IDiscretization> disc;
  std::unique_ptr<TestProblem> problem;
  std::unique_ptr<TestAssemble> assemble;

  TestMonomialProblems(CELLTYPE cellType) : cellType(cellType) { degree = GetParam(); }

  void RunTestLaplace_X() { runTestLaplaceT<0>(); }

  void RunTestLaplace_Y() { runTestLaplaceT<1>(); }

  void RunTestLaplace_Z() { runTestLaplaceT<2>(); }

  void RunTestBiLaplace() {
    if (degree > 0) {
      mout << "Test for Argyris element is independent of degree and therefore only performed for "
              "degree 0"
           << endl;
      return;
    }
    runTestBiLaplaceT<false>();
    mout << endl;
    runTestBiLaplaceT<true>();
  }
private:
  void initMesh(std::string name, int level) {
    meshes = MeshesCreator(name).WithPLevel(0).WithLevel(level).CreateUnique();
  }

  template<int i>
  void runTestLaplaceT() {
    if (degree < 2) {
      mout << "No test implemented for degree " + std::to_string(degree) << endl;
      return;
    }
    switch (cellType) {
    case INTERVAL:
      initMesh("Interval", 0);
      break;
    case TRIANGLE:
      initMesh("Triangle", 0);
      break;
    case QUADRILATERAL:
      initMesh("Square", 0);
      break;
    case TETRAHEDRON:
      initMesh("Tetrahedron", 0);
      break;
    case HEXAHEDRON:
      initMesh("Hexahedron", 0);
      break;
    }
    disc = std::make_unique<LagrangeDiscretization>(*meshes, degree);
    problem = std::make_unique<TestProblemLaplaceMonomialT<i>>(degree);
    assemble = std::make_unique<LaplaceAssemble>(*problem);

    Vector u(0.0, disc);
    NewtonMethod(*assemble, u);
    ASSERT_NEAR(assemble->L2Error(u), 0.0, 1e-9);
  }

  template<bool turned>
  void runTestBiLaplaceT() {
    if (cellType != TRIANGLE) THROW("No test implemented for celltype " + std::to_string(cellType));

    if constexpr (!turned) {
      initMesh("Square4Triangles", 4);
      mout << "Running BiLaplace on mesh 'Square4Triangles'" << endl;
    } else {
      initMesh("Square4TrianglesTurned", 4);
      mout << "Running BiLaplace on mesh 'Square4TrianglesTurned'" << endl;
    }
    disc = std::make_shared<ArgyrisDiscretization>(*meshes);
    std::list<Point> corners;
    if constexpr (!turned) {
      problem = std::make_unique<TestProblemBiLaplace>();
      corners = {Point(0.0, 0.0), Point(1.0, 0.0), Point(1.0, 1.0), Point(0.0, 1.0)};
    } else {
      problem = std::make_unique<TestProblemBiLaplaceTurned>();
      corners = {Point(0.70710678118, 0.0), Point(0.0, 0.70710678118), Point(-0.70710678118, 0.0),
                 Point(0.0, -0.70710678118)};
    }
    assemble = std::make_unique<BiLaplaceAssemble>(*problem, corners);

    Vector u(0.0, disc);
    NewtonMethod(*assemble, u);
    ASSERT_NEAR(assemble->L2Error(u), 0.0, 2e-3); // TODO: Check why the error is still that
                                                  // large!!!
  }
};

class TestInterval : public TestMonomialProblems {
protected:
  TestInterval() : TestMonomialProblems(INTERVAL) {}
};

TEST_P(TestInterval, MonomialLaplaceTest) { RunTestLaplace_X(); }

INSTANTIATE_TEST_SUITE_P(TestPolynomialProblems, TestInterval, Values(0, 1, 2, 3, 4, 5, 6, 7, 8));


#if SpaceDimension >= 2

class TestTriangle : public TestMonomialProblems {
protected:
  TestTriangle() : TestMonomialProblems(TRIANGLE) {}
};

TEST_P(TestTriangle, MonomialLaplaceTest_X) { RunTestLaplace_X(); }

TEST_P(TestTriangle, MonomialLaplaceTest_Y) { RunTestLaplace_Y(); }

TEST_P(TestTriangle, BiLaplaceTest) { RunTestBiLaplace(); }

INSTANTIATE_TEST_SUITE_P(TestPolynomialProblems, TestTriangle, Values(0, 1, 2, 3, 4, 5, 6));

class TestQuadrilateral : public TestMonomialProblems {
protected:
  TestQuadrilateral() : TestMonomialProblems(QUADRILATERAL) {}
};

TEST_P(TestQuadrilateral, MonomialLaplaceTest_X) { RunTestLaplace_X(); }

TEST_P(TestQuadrilateral, MonomialLaplaceTest_Y) { RunTestLaplace_Y(); }

INSTANTIATE_TEST_SUITE_P(TestPolynomialProblems, TestQuadrilateral,
                         Values(0, 1, 2, 3, 4, 5, 6, 7, 8));

#endif


#if SpaceDimension >= 3

class TestTetrahedron : public TestMonomialProblems {
protected:
  TestTetrahedron() : TestMonomialProblems(TETRAHEDRON) {}
};

TEST_P(TestTetrahedron, MonomialLaplaceTest_X) { RunTestLaplace_X(); }

TEST_P(TestTetrahedron, MonomialLaplaceTest_Y) { RunTestLaplace_Y(); }

TEST_P(TestTetrahedron, MonomialLaplaceTest_Z) { RunTestLaplace_Z(); }

INSTANTIATE_TEST_SUITE_P(TestPolynomialProblems, TestTetrahedron, Values(0, 1, 2, 3, 4));

class TestHexahedron : public TestMonomialProblems {
protected:
  TestHexahedron() : TestMonomialProblems(HEXAHEDRON) {}
};

TEST_P(TestHexahedron, MonomialLaplaceTest_X) { RunTestLaplace_X(); }

TEST_P(TestHexahedron, MonomialLaplaceTest_Y) { RunTestLaplace_Y(); }

TEST_P(TestHexahedron, MonomialLaplaceTest_Z) { RunTestLaplace_Z(); }

INSTANTIATE_TEST_SUITE_P(TestPolynomialProblems, TestHexahedron, Values(0, 1, 2, 3, 4, 5, 6, 7, 8));

#endif

int main(int argc, char **argv) {
  return MppTest(MppTestBuilder(argc, argv)
                     .WithConfigEntry("NewtonVerbose", 1)
                     .WithConfigEntry("LinearVerbose", 1)
                     .WithScreenLogging()
                     .WithPPM())
      .RUN_ALL_MPP_TESTS();
}