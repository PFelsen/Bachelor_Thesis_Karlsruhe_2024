#include "gtest/gtest.h"

#include "Cholesky.hpp"
#include "TestEnvironment.hpp"
#include "basicalgebra/RMatrix.hpp"

#include <vector>
using namespace std;

const double tol = 1e-14;

struct BasicAlgebraCholeskyTestParameter {
  int dim;
  vector<vector<vector<double>>> A;
  vector<vector<vector<double>>> expectedL;
  vector<vector<double>> b;
  vector<vector<double>> expected_x;
};

static vector<BasicAlgebraCholeskyTestParameter> TestCases_real =
    {BasicAlgebraCholeskyTestParameter{2,
                                       {{{1}, {2}}, {{2}, {13}}},
                                       {{{1}}, {{2}, {3}}},
                                       {{-3}, {-24}},
                                       {{1}, {-2}}},
     BasicAlgebraCholeskyTestParameter{2,
                                       {{{9}, {-3}}, {{-3}, {5}}},
                                       {{{3}}, {{-1}, {2}}},
                                       {{6}, {2}},
                                       {{1}, {1}}},
     BasicAlgebraCholeskyTestParameter{3,
                                       {{{4}, {2}, {0}}, {{2}, {5}, {2}}, {{0}, {2}, {5}}},
                                       {{{2}}, {{1}, {2}}, {{0}, {1}, {2}}},
                                       {{0}, {-2}, {11}},
                                       {{1}, {-2}, {3}}},
     BasicAlgebraCholeskyTestParameter{3,
                                       {{{4}, {6}, {8}}, {{6}, {10}, {17}}, {{8}, {17}, {45}}},
                                       {{{2}}, {{3}, {1}}, {{4}, {5}, {2}}},
                                       {{38}, {73}, {168}},
                                       {{2}, {1}, {3}}},
     BasicAlgebraCholeskyTestParameter{4,
                                       {{{9}, {6}, {-15}, {3}},
                                        {{6}, {5}, {-12}, {2}},
                                        {{-15}, {-12}, {45}, {-13}},
                                        {{3}, {2}, {-13}, {6}}},
                                       {{{3}},
                                        {{2}, {1}},
                                        {{-5}, {-2}, {4}},
                                        {{1}, {0}, {-2}, {1}}},
                                       {{0}, {-1}, {10}, {-5}},
                                       {{1}, {-1}, {0}, {-1}}},
     BasicAlgebraCholeskyTestParameter{4,
                                       {{{9}, {3}, {-6}, {12}},
                                        {{3}, {26}, {-7}, {-11}},
                                        {{-6}, {-7}, {9}, {7}},
                                        {{12}, {-11}, {7}, {65}}},
                                       {{{3}},
                                        {{1}, {5}},
                                        {{-2}, {-1}, {2}},
                                        {{4}, {-3}, {6}, {2}}},
                                       {{72}, {34}, {22}, {326}},
                                       {{2}, {4}, {3}, {5}}},
     BasicAlgebraCholeskyTestParameter{4,
                                       {{{4}, {2}, {4}, {4}},
                                        {{2}, {10}, {5}, {2}},
                                        {{4}, {5}, {9}, {6}},
                                        {{4}, {2}, {6}, {9}}},
                                       {{{2}}, {{1}, {3}}, {{2}, {1}, {2}}, {{2}, {0}, {1}, {2}}},
                                       {{2}, {-5}, {2}, {-1}},
                                       {{1}, {-1}, {1}, {-1}}}};

static vector<BasicAlgebraCholeskyTestParameter> TestCases_complex = {
    BasicAlgebraCholeskyTestParameter{2,
                                      {{{1, 0}, {0, -1}}, {{0, 1}, {5, 0}}},
                                      {{{1, 0}}, {{0, 1}, {2, 0}}},
                                      {{1, -1}, {5, 1}},
                                      {{1, 0}, {1, 0}}}};

static vector<BasicAlgebraCholeskyTestParameter> TestCases_all;

void combineTestCases() {
  TestCases_all = TestCases_real;
  TestCases_all.insert(TestCases_all.end(), TestCases_complex.begin(), TestCases_complex.end());
}

class BasicAlgebraCholeskyTest : public TestWithParam<BasicAlgebraCholeskyTestParameter> {
protected:
  int dim;

  BasicAlgebraCholeskyTest() { dim = GetParam().dim; }
};

class BasicAlgebraRCholeskyTest : public BasicAlgebraCholeskyTest {
protected:
  RMatrix expectedL;
  RVector b, expected_x;

  BasicAlgebraRCholeskyTest() {
    expectedL.resize(dim);
    b.resize(dim);
    expected_x.resize(dim);
    for (int i = 0; i < dim; ++i) {
      for (int j = 0; j <= i; ++j)
        expectedL[i][j] = (GetParam().expectedL)[i][j][0];
      b[i] = (GetParam().b)[i][0];
      expected_x[i] = (GetParam().expected_x)[i][0];
    }
  }
};

class BasicAlgebraSymRMatrixCholeskyTest : public BasicAlgebraRCholeskyTest {
protected:
  SymRMatrix A;
  RMatrix L;

  BasicAlgebraSymRMatrixCholeskyTest() : BasicAlgebraRCholeskyTest() {
    A.resize(dim);
    for (int i = 0; i < dim; ++i)
      for (int j = 0; j <= i; ++j)
        A(i, j) = (GetParam().A)[i][j][0];
    L = Cholesky(A);
  }
};

class BasicAlgebraRMatrixCholeskyTest : public BasicAlgebraRCholeskyTest {
protected:
  RMatrix A;
  RMatrix L;

  BasicAlgebraRMatrixCholeskyTest() : BasicAlgebraRCholeskyTest() {
    A.resize(dim);
    for (int i = 0; i < dim; ++i)
      for (int j = 0; j < dim; ++j)
        A(i, j) = (GetParam().A)[i][j][0];
    L = Cholesky(A);
  }
};

TEST_P(BasicAlgebraSymRMatrixCholeskyTest, CholeskyLTest) {
  for (int i = 0; i < dim; ++i)
    for (int j = 0; j < dim; ++j)
      EXPECT_NEAR(L[i][j], expectedL[i][j], tol);
}

TEST_P(BasicAlgebraSymRMatrixCholeskyTest, CholeskyMultiplicationTest) {
  RMatrix B = L * transpose(L);
  for (int i = 0; i < dim; ++i)
    for (int j = 0; j < dim; ++j)
      EXPECT_NEAR(A(i, j), B[i][j], tol);
}

TEST_P(BasicAlgebraSymRMatrixCholeskyTest, LinearSystemTest) {
  RVector x = linearSystem(A, b);
  for (int i = 0; i < dim; ++i)
    EXPECT_NEAR(x[i], expected_x[i], tol);
}

INSTANTIATE_TEST_CASE_P(BasicAlgebraCholeskyTest, BasicAlgebraSymRMatrixCholeskyTest,
                        ValuesIn(TestCases_real));

TEST_P(BasicAlgebraRMatrixCholeskyTest, CholeskyLTest) {
  for (int i = 0; i < dim; ++i)
    for (int j = 0; j < dim; ++j)
      EXPECT_NEAR(L[i][j], expectedL[i][j], tol);
}

TEST_P(BasicAlgebraRMatrixCholeskyTest, CholeskyMultiplicationTest) {
  RMatrix B = L * transpose(L);
  for (int i = 0; i < dim; ++i)
    for (int j = 0; j < dim; ++j)
      EXPECT_NEAR(A(i, j), B[i][j], tol);
}

TEST_P(BasicAlgebraRMatrixCholeskyTest, LinearSystemTest) {
  RVector x = linearSystem(A, b);
  for (int i = 0; i < dim; ++i)
    EXPECT_NEAR(x[i], expected_x[i], tol);
}

INSTANTIATE_TEST_CASE_P(BasicAlgebraCholeskyTest, BasicAlgebraRMatrixCholeskyTest,
                        ValuesIn(TestCases_real));

class BasicAlgebraCCholeskyTest : public BasicAlgebraCholeskyTest {
protected:
  CMatrix expectedL;
  CVector b, expected_x;

  BasicAlgebraCCholeskyTest() {
    expectedL.resize(dim);
    b.resize(dim);
    expected_x.resize(dim);
    for (int i = 0; i < dim; ++i) {
      for (int j = 0; j <= i; ++j) {
        vector<double> entry = (GetParam().expectedL)[i][j];
        if (entry.size() == 1) expectedL[i][j] = entry[0];
        else expectedL[i][j] = std::complex<double>(entry[0], entry[1]);
      }
      vector<double> entry_b = (GetParam().b)[i];
      vector<double> entry_x = (GetParam().expected_x)[i];
      if (entry_b.size() == 1) b[i] = entry_b[0];
      else b[i] = std::complex<double>(entry_b[0], entry_b[1]);
      if (entry_x.size() == 1) expected_x[i] = entry_x[0];
      else expected_x[i] = std::complex<double>(entry_x[0], entry_x[1]);
    }
  }
};

class BasicAlgebraHermCMatrixCholeskyTest : public BasicAlgebraCCholeskyTest {
protected:
  HermCMatrix A;
  CMatrix L;

  BasicAlgebraHermCMatrixCholeskyTest() : BasicAlgebraCCholeskyTest() {
    A.resize(dim);
    for (int i = 0; i < dim; ++i)
      for (int j = 0; j <= i; ++j) {
        vector<double> entry = (GetParam().A)[i][j];
        if (entry.size() == 1) A(entry[0], i, j);
        else A(std::complex<double>(entry[0], entry[1]), i, j);
      }
    L = Cholesky(A);
  }
};

class BasicAlgebraCMatrixCholeskyTest : public BasicAlgebraCCholeskyTest {
protected:
  CMatrix A;
  CMatrix L;

  BasicAlgebraCMatrixCholeskyTest() : BasicAlgebraCCholeskyTest() {
    A.resize(dim);
    for (int i = 0; i < dim; ++i)
      for (int j = 0; j < dim; ++j) {
        vector<double> entry = (GetParam().A)[i][j];
        if (entry.size() == 1) A[i][j] = entry[0];
        else A[i][j] = std::complex<double>(entry[0], entry[1]);
      }
    L = Cholesky(A);
  }
};

TEST_P(BasicAlgebraHermCMatrixCholeskyTest, CholeskyLTest) {
  for (int i = 0; i < dim; ++i)
    for (int j = 0; j < dim; ++j) {
      EXPECT_NEAR(std::real(L[i][j]), std::real(expectedL[i][j]), tol);
      EXPECT_NEAR(std::imag(L[i][j]), std::imag(expectedL[i][j]), tol);
    }
}

TEST_P(BasicAlgebraHermCMatrixCholeskyTest, CholeskyMultiplicationTest) {
  CMatrix B = L * adjoint(L);
  for (int i = 0; i < dim; ++i)
    for (int j = 0; j < dim; ++j) {
      EXPECT_NEAR(std::real(A(i, j)), std::real(B[i][j]), tol);
      EXPECT_NEAR(std::imag(A(i, j)), std::imag(B[i][j]), tol);
    }
}

TEST_P(BasicAlgebraHermCMatrixCholeskyTest, LinearSystemTest) {
  CVector x = linearSystem(A, b);
  for (int i = 0; i < dim; ++i) {
    EXPECT_NEAR(std::real(x[i]), std::real(expected_x[i]), tol);
    EXPECT_NEAR(std::imag(x[i]), std::imag(expected_x[i]), tol);
  }
}

INSTANTIATE_TEST_CASE_P(BasicAlgebraCholeskyTest, BasicAlgebraHermCMatrixCholeskyTest,
                        ValuesIn(TestCases_all));

TEST_P(BasicAlgebraCMatrixCholeskyTest, CholeskyLTest) {
  for (int i = 0; i < dim; ++i)
    for (int j = 0; j < dim; ++j) {
      EXPECT_NEAR(std::real(L[i][j]), std::real(expectedL[i][j]), tol);
      EXPECT_NEAR(std::imag(L[i][j]), std::imag(expectedL[i][j]), tol);
    }
}

TEST_P(BasicAlgebraCMatrixCholeskyTest, CholeskyMultiplicationTest) {
  CMatrix B = L * adjoint(L);
  for (int i = 0; i < dim; ++i)
    for (int j = 0; j < dim; ++j) {
      EXPECT_NEAR(std::real(A(i, j)), std::real(B[i][j]), tol);
      EXPECT_NEAR(std::imag(A(i, j)), std::imag(B[i][j]), tol);
    }
}

TEST_P(BasicAlgebraCMatrixCholeskyTest, LinearSystemTest) {
  CVector x = linearSystem(A, b);
  for (int i = 0; i < dim; ++i) {
    EXPECT_NEAR(std::real(x[i]), std::real(expected_x[i]), tol);
    EXPECT_NEAR(std::imag(x[i]), std::imag(expected_x[i]), tol);
  }
}

INSTANTIATE_TEST_CASE_P(BasicAlgebraCholeskyTest, BasicAlgebraCMatrixCholeskyTest,
                        ValuesIn(TestCases_all));

int main(int argc, char **argv) {
  combineTestCases();
  MppTest mppTest = MppTestBuilder(argc, argv).WithPPM();
  return mppTest.RUN_ALL_MPP_TESTS();
}
