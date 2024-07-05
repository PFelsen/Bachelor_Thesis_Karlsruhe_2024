#include "TestEigenSolver.hpp"
#include "ArgyrisDiscretization.hpp"
#include "ArgyrisElement.hpp"
#include "LagrangeDiscretization.hpp"
#include "ScalarElement.hpp"
#include "SerendipityDiscretization.hpp"

using namespace std;

class EigenSolverTest : public TestWithParam<EigenSolverTestParameter> {
protected:
  std::shared_ptr<const Meshes> meshes;
  std::shared_ptr<const IDiscretization> disc;
  IEigenSolver *esolver = nullptr;

  vector<double> expectedLambda;

  EigenSolverTest(string eigensolverName) {
    meshes =
        MeshesCreator(GetParam().meshName).WithLevel(GetParam().level).WithPLevel(1).CreateShared();

    if (GetParam().discName == "Lagrange") {
      disc = std::make_shared<const LagrangeDiscretization>(*meshes, GetParam().degree);
    } else if (GetParam().discName == "Serendipity") {
      disc = std::make_shared<const SerendipityDiscretization>(*meshes);
    } else if (GetParam().discName == "Argyris") {
      disc = std::make_shared<const ArgyrisDiscretization>(*meshes);
    } else Exit("Not implemented") esolver = EigenSolverCreator(eigensolverName).Create();

    expectedLambda.resize(6);
    if (GetParam().discName == "Lagrange" || GetParam().discName == "Serendipity") {
      expectedLambda[0] = 2 * Pi * Pi;
      expectedLambda[1] = expectedLambda[2] = 5 * Pi * Pi;
      expectedLambda[3] = 8 * Pi * Pi;
      expectedLambda[4] = expectedLambda[5] = 10 * Pi * Pi;
    } else if (GetParam().discName == "Argyris") {
      expectedLambda[0] = 1294.93394;
      expectedLambda[1] = expectedLambda[2] = 5386.6564;
      expectedLambda[3] = 11710.8103;
      expectedLambda[4] = 17313.49969;
      expectedLambda[5] = 17478.104;
    }
  }

  ~EigenSolverTest() override {
    if (esolver) delete esolver;
  }

  void Dirichlet(Vectors &U) {
    for (int n = 0; n < U.size(); ++n) {
      if (GetParam().discName == "Lagrange" || GetParam().discName == "Serendipity") {
        ScalarElement::H10BC(U[n]);
      } else if (GetParam().discName == "Argyris") {
        ArgyrisElement::H20BC(GetParam().corners, U[n]);
      }
    }
  }

  void MatrixentriesScalar(Matrix &A, Matrix &B) {
    A = 0;
    B = 0;
    for (cell c = A.cells(); c != A.cells_end(); ++c) {
      ScalarElement E(A, *c);
      RowEntries A_c(A, *c);
      RowEntries B_c(B, *c);

      for (int q = 0; q < E.nQ(); ++q)
        for (int i = 0; i < E.size(); ++i)
          for (int j = 0; j < E.size(); ++j) {
            A_c(i, j) += E.QWeight(q) * E.Derivative(q, i) * E.Derivative(q, j);
            B_c(i, j) += E.QWeight(q) * E.Value(q, i) * E.Value(q, j);
          }
    }
    A.ClearDirichletValues();
    B.ClearDirichletValues();
  }

  void MatrixentriesArgyris(Matrix &A, Matrix &B) {
    A = 0;
    B = 0;
    for (cell c = A.cells(); c != A.cells_end(); ++c) {
      ArgyrisElement E(A, *c);
      RowEntries A_c(A, *c);
      RowEntries B_c(B, *c);

      for (int q = 0; q < E.nQ(); ++q)
        for (int i = 0; i < E.size(); ++i)
          for (int j = 0; j < E.size(); ++j)
            for (int k = 0; k < E.get_maxk(i); ++k)
              for (int l = 0; l < E.get_maxk(j); ++l) {
                A_c(i, j, k, l) += E.QWeight(q) * E.Laplace(q, i, k) * E.Laplace(q, j, l);
                B_c(i, j, k, l) += E.QWeight(q) * E.Value(q, i, k) * E.Value(q, j, l);
              }
    }
    A.ClearDirichletValues();
    B.ClearDirichletValues();
  }

  void test() {
    Eigenfcts U(6, 0.0, disc);
    Dirichlet(U);

    Matrix A(U[0]);
    Matrix B(U[0]);
    if (GetParam().discName == "Lagrange" || GetParam().discName == "Serendipity") {
      MatrixentriesScalar(A, B);
    } else if (GetParam().discName == "Argyris") {
      MatrixentriesArgyris(A, B);
    }

    Eigenvalues lambda;
    (*esolver)(U, lambda, A, B);

    for (int i = 0; i < U.size(); ++i) {
      mout << lambda[i] << endl;
      EXPECT_NEAR(lambda[i], expectedLambda[i], 3e-3);
    }
  }
};

class RitzGalerkinTest : public EigenSolverTest {
protected:
  RitzGalerkinTest() : EigenSolverTest("RitzGalerkin") {}
};

TEST_P(RitzGalerkinTest, Problem) { test(); }

INSTANTIATE_TEST_SUITE_P(EigenSolverTest, RitzGalerkinTest, ValuesIn(TestCases));

class LOBPCGTest : public EigenSolverTest {
protected:
  LOBPCGTest() : EigenSolverTest("LOBPCG") {}
};

TEST_P(LOBPCGTest, Problem) { test(); }

INSTANTIATE_TEST_SUITE_P(EigenSolverTest, LOBPCGTest, ValuesIn(TestCases));

class LOBPCGExtendedTest : public EigenSolverTest {
protected:
  LOBPCGExtendedTest() : EigenSolverTest("LOBPCGExtended") {}
};

TEST_P(LOBPCGExtendedTest, Problem) { test(); }

INSTANTIATE_TEST_SUITE_P(EigenSolverTest, LOBPCGExtendedTest, ValuesIn(TestCases));

class LOBPCGSelectiveTest : public EigenSolverTest {
protected:
  LOBPCGSelectiveTest() : EigenSolverTest("LOBPCGSelective") {}
};

TEST_P(LOBPCGSelectiveTest, Problem) { test(); }

INSTANTIATE_TEST_SUITE_P(EigenSolverTest, LOBPCGSelectiveTest, ValuesIn(TestCases));

class LOBPCGSelectiveFastTest : public EigenSolverTest {
protected:
  LOBPCGSelectiveFastTest() : EigenSolverTest("LOBPCGSelectiveFast") {}
};

TEST_P(LOBPCGSelectiveFastTest, Problem) { test(); }

INSTANTIATE_TEST_SUITE_P(EigenSolverTest, LOBPCGSelectiveFastTest, ValuesIn(TestCases));

class LOBPCGSelectiveVeryFastTest : public EigenSolverTest {
protected:
  LOBPCGSelectiveVeryFastTest() : EigenSolverTest("LOBPCGSelectiveVeryFast") {}
};

TEST_P(LOBPCGSelectiveVeryFastTest, Problem) { test(); }

INSTANTIATE_TEST_SUITE_P(EigenSolverTest, LOBPCGSelectiveVeryFastTest, ValuesIn(TestCases));

int main(int argc, char **argv) {
  MppTest mppTest = MppTestBuilder(argc, argv).WithPPM().WithScreenLogging();
  return mppTest.RUN_ALL_MPP_TESTS();
}
