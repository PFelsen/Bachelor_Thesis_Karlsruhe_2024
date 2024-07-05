#include "BasicSolver.hpp"
#include "TestEnvironment.hpp"

RMatrix rMatrix() {
  return RMatrix{{1, 2, 3, 4, 5, 7},  {-1, 2, 0, 5, 7, -1}, {9, -3, 6, 1, 0, -4},
                 {3, -4, 2, 0, 5, 8}, {-1, 1, -1, 2, 9, 2}, {3, 0, -4, 2, 1, -1}};
}

SymRMatrix sMatrix() {
  SymRMatrix A(6);
  A(1, 0, 0);
  A(-3, 1, 0);
  A(5, 1, 1);
  A(4, 2, 0);
  A(-9, 2, 1);
  A(-2, 2, 2);
  A(-4, 3, 0);
  A(2, 3, 1);
  A(3, 3, 2);
  A(1, 3, 3);
  A(2, 4, 0);
  A(-6, 4, 1);
  A(-1, 4, 2);
  A(0, 4, 3);
  A(-3, 4, 4);
  A(-1, 5, 0);
  A(0, 5, 1);
  A(3, 5, 2);
  A(7, 5, 3);
  A(-9, 5, 4);
  A(2, 5, 5);
  return A;
}

AntisymRMatrix aMatrix() {
  AntisymRMatrix A(6);
  A(-3, 1, 0);
  A(4, 2, 0);
  A(-9, 2, 1);
  A(-4, 3, 0);
  A(2, 3, 1);
  A(3, 3, 2);
  A(2, 4, 0);
  A(-6, 4, 1);
  A(-1, 4, 2);
  A(0, 4, 3);
  A(-1, 5, 0);
  A(0, 5, 1);
  A(3, 5, 2);
  A(7, 5, 3);
  A(-9, 5, 4);
  return A;
}

CMatrix cMatrix() { return CMatrix(rMatrix()); };

HermCMatrix hMatrix() { return HermCMatrix(sMatrix()); };

template<template<template<typename> class MATRIX, typename T> class SOLVER>
class BasicSolverTest : public ::testing::Test {
protected:
  RMatrix A_r;
  SymRMatrix A_s;
  AntisymRMatrix A_a;
  CMatrix A_c;
  HermCMatrix A_h;
  std::vector<RVector> sol_r;
  RMatrix SOL_r;
  std::vector<CVector> sol_c;
  CMatrix SOL_c;

  BasicSolverTest() :
      A_r(rMatrix()), A_s(sMatrix()), A_a(aMatrix()), A_c(cMatrix()), A_h(hMatrix()),
      sol_r(10, RVector(6)), SOL_r(6, 10), sol_c(10, CVector(6)), SOL_c(6, 10) {
    for (int k = 0; k < sol_r.size(); ++k)
      for (int i = 0; i < sol_r[k].size(); ++i) {
        sol_r[k][i] = RandomDouble();
        sol_c[k][i] = RandomComplex();
      }
    for (int k = 0; k < SOL_r.rows(); ++k)
      for (int i = 0; i < SOL_r.cols(); ++i) {
        SOL_r(k, i) = RandomDouble();
        SOL_c(k, i) = RandomComplex();
      }
  }

  void checkVectorSolveR() {
    for (int k = 0; k < sol_r.size(); ++k) {
      RVector rhs = A_r * sol_r[k];
      SOLVER solver(A_r);
      RVector x = solver.Solve(rhs);
      EXPECT_EQ(x, sol_r[k]);
    }
  }

  void checkVectorSolveSymR() {
    for (int k = 0; k < sol_r.size(); ++k) {
      RVector rhs = A_s * sol_r[k];
      SOLVER solver(A_s);
      RVector x = solver.Solve(rhs);
      EXPECT_EQ(x, sol_r[k]);
    }
  }

  void checkVectorSolveAntisymR() {
    for (int k = 0; k < sol_r.size(); ++k) {
      RVector rhs = A_a * sol_r[k];
      SOLVER solver(A_a);
      RVector x = solver.Solve(rhs);
      EXPECT_EQ(x, sol_r[k]);
    }
  }

  void checkVectorSolveC() {
    for (int k = 0; k < sol_c.size(); ++k) {
      CVector rhs = A_c * sol_c[k];
      SOLVER solver(A_c);
      CVector x = solver.Solve(rhs);
      EXPECT_EQ(x, sol_c[k]);
    }
  }

  void checkVectorSolveHermC() {
    for (int k = 0; k < sol_c.size(); ++k) {
      CVector rhs = A_h * sol_c[k];
      SOLVER solver(A_h);
      CVector x = solver.Solve(rhs);
      EXPECT_EQ(x, sol_c[k]);
    }
  }

  void checkMatrixSolveR() {
    RMatrix RHS = A_r * SOL_r;
    SOLVER solver(A_r);
    RMatrix X = solver.Solve(RHS);
    EXPECT_EQ(X, SOL_r);
  }

  void checkMatrixSolveSymR() {
    RMatrix RHS = A_s * SOL_r;
    SOLVER solver(A_s);
    RMatrix X = solver.Solve(RHS);
    EXPECT_EQ(X, SOL_r);
  }

  void checkMatrixSolveAntisymR() {
    RMatrix RHS = A_a * SOL_r;
    SOLVER solver(A_a);
    RMatrix X = solver.Solve(RHS);
    EXPECT_EQ(X, SOL_r);
  }

  void checkMatrixSolveC() {
    CMatrix RHS = A_c * SOL_c;
    SOLVER solver(A_c);
    CMatrix X = solver.Solve(RHS);
    EXPECT_EQ(X, SOL_c);
  }

  void checkMatrixSolveHermC() {
    CMatrix RHS = A_h * SOL_c;
    SOLVER solver(A_h);
    CMatrix X = solver.Solve(RHS);
    EXPECT_EQ(X, SOL_c);
  }
};

/// To avoid the same code in each test
#define BASICSOLVER_TESTS(testClass)                                                               \
                                                                                                   \
  TEST_F(testClass, VectorSolveRTest) { checkVectorSolveR(); }                                     \
                                                                                                   \
  TEST_F(testClass, VectorSolveSymRTest) { checkVectorSolveSymR(); }                               \
                                                                                                   \
  TEST_F(testClass, VectorSolveAntisymRTest) { checkVectorSolveAntisymR(); }                       \
                                                                                                   \
  TEST_F(testClass, VectorSolveCTest) { checkVectorSolveC(); }                                     \
                                                                                                   \
  TEST_F(testClass, VectorSolveHermCTest) { checkVectorSolveHermC(); }                             \
                                                                                                   \
  TEST_F(testClass, MatrixSolveRTest) { checkMatrixSolveR(); }                                     \
                                                                                                   \
  TEST_F(testClass, MatrixSolveSymRTest) { checkMatrixSolveSymR(); }                               \
                                                                                                   \
  TEST_F(testClass, MatrixSolveAntisymRTest) { checkMatrixSolveAntisymR(); }                       \
                                                                                                   \
  TEST_F(testClass, MatrixSolveCTest) { checkMatrixSolveC(); }                                     \
                                                                                                   \
  TEST_F(testClass, MatrixSolveHermCTest) { checkMatrixSolveHermC(); }

using DirectBasicSolverTest = BasicSolverTest<DirectLinearSolverT>;

BASICSOLVER_TESTS(DirectBasicSolverTest)

using LAPACKBasicSolverTest = BasicSolverTest<LAPACKLinearSolverT>;

BASICSOLVER_TESTS(LAPACKBasicSolverTest)

int main(int argc, char **argv) {
  MppTest mppTest = MppTestBuilder(argc, argv);
  return mppTest.RUN_ALL_MPP_TESTS();
}
