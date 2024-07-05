#include "CMatrix.hpp"
#include "TestEnvironment.hpp"

#include <vector>

template<typename T>
struct TestCases {
  int dim;
  std::vector<std::vector<T>> A;
  std::vector<std::vector<T>> Ainv;
};

class RMatrixTest : public TestWithParam<TestCases<double>> {
protected:
  int dim;
  RMatrix A;
  RMatrix Ainv;

  RMatrixTest() : dim(GetParam().dim), A(GetParam().dim), Ainv(GetParam().dim) {
    for (int i = 0; i < dim; ++i)
      for (int j = 0; j < dim; ++j) {
        A(i, j) = GetParam().A[i][j];
        Ainv(i, j) = GetParam().Ainv[i][j];
      }
  }
};

TEST_P(RMatrixTest, InvertTest) {
  mpp_ba::SetTolerance(1e-14);
  RMatrix B(A);
  EXPECT_EQ(B.Invert(), Ainv);
  RMatrix C(Ainv);
  EXPECT_EQ(C.Invert(), A);
}

INSTANTIATE_TEST_CASE_P(
    BasicAlgebraTest, RMatrixTest,
    Values(TestCases<double>{3,
                             {{-1, 0, 2}, {-2, -1, 5}, {3, 1, -8}},
                             {{-3, -2, -2}, {1, -2, -1}, {-1, -1, -1}}},
           TestCases<double>{3,
                             {{1, 0, -2}, {0, -1, 3}, {0, 0, -3}},
                             {{1, 0, -2.0 / 3.0}, {0, -1, -1}, {0, 0, -1.0 / 3.0}}},
           TestCases<double>{3,
                             {{1, -1, 2}, {-1, 0, -2}, {2, -2, 2}},
                             {{-2, -1, 1}, {-1, -1, 0}, {1, 0, -0.5}}},
           TestCases<double>{3,
                             {{1, 2, -1}, {2, 3, -2}, {-1, -2, 2}},
                             {{-2, 2, 1}, {2, -1, 0}, {1, 0, 1}}}));

class SymRMatrixTest : public TestWithParam<TestCases<double>> {
protected:
  int dim;
  SymRMatrix A;
  SymRMatrix Ainv;

  SymRMatrixTest() : dim(GetParam().dim), A(GetParam().dim), Ainv(GetParam().dim) {
    for (int i = 0; i < dim; ++i)
      for (int j = 0; j < dim; ++j) {
        A(i, j) = GetParam().A[i][j];
        Ainv(i, j) = GetParam().Ainv[i][j];
      }
  }
};

TEST_P(SymRMatrixTest, InvertTest) {
  mpp_ba::SetTolerance(1e-14);
  SymRMatrix B(A);
  EXPECT_EQ(B.Invert(), Ainv);
  SymRMatrix C(Ainv);
  EXPECT_EQ(C.Invert(), A);
}

INSTANTIATE_TEST_CASE_P(BasicAlgebraTest, SymRMatrixTest,
                        Values(TestCases<double>{3,
                                                 {{1, -1, 2}, {-1, 0, -2}, {2, -2, 2}},
                                                 {{-2, -1, 2}, {-1, -1, 0}, {1, 0, -0.5}}},
                               TestCases<double>{3,
                                                 {{1, 2, -1}, {2, 3, -2}, {-1, -2, 2}},
                                                 {{-2, 2, 1}, {2, -1, 0}, {1, 0, 1}}}));

class CMatrixTest : public TestWithParam<TestCases<std::complex<double>>> {
protected:
  int dim;
  CMatrix A;
  CMatrix Ainv;

  CMatrixTest() : dim(GetParam().dim), A(GetParam().dim), Ainv(GetParam().dim) {
    for (int i = 0; i < dim; ++i)
      for (int j = 0; j < dim; ++j) {
        A(i, j) = GetParam().A[i][j];
        Ainv(i, j) = GetParam().Ainv[i][j];
      }
  }
};

TEST_P(CMatrixTest, InvertTest) {
  mpp_ba::SetTolerance(1e-14);
  CMatrix B(A);
  EXPECT_EQ(B.Invert(), Ainv);
  CMatrix C(Ainv);
  EXPECT_EQ(C.Invert(), A);
}

INSTANTIATE_TEST_CASE_P(
    BasicAlgebraTest, CMatrixTest,
    Values(TestCases<std::complex<double>>{3,
                                           {{-1, 0, 2}, {-2, -1, 5}, {3, 1, -8}},
                                           {{-3, -2, -2}, {1, -2, -1}, {-1, -1, -1}}},
           TestCases<std::complex<double>>{3,
                                           {{1, 0, -2}, {0, -1, 3}, {0, 0, -3}},
                                           {{1, 0, -2.0 / 3.0}, {0, -1, -1}, {0, 0, -1.0 / 3.0}}},
           TestCases<std::complex<double>>{3,
                                           {{1, -1, 2}, {-1, 0, -2}, {2, -2, 2}},
                                           {{-2, -1, 1}, {-1, -1, 0}, {1, 0, -0.5}}},
           TestCases<std::complex<double>>{3,
                                           {{1, 2, -1}, {2, 3, -2}, {-1, -2, 2}},
                                           {{-2, 2, 1}, {2, -1, 0}, {1, 0, 1}}}));

class HermCMatrixTest : public TestWithParam<TestCases<std::complex<double>>> {
protected:
  int dim;
  HermCMatrix A;
  HermCMatrix Ainv;

  HermCMatrixTest() : dim(GetParam().dim), A(GetParam().dim), Ainv(GetParam().dim) {
    for (int i = 0; i < dim; ++i)
      for (int j = 0; j < dim; ++j) {
        A(GetParam().A[i][j], i, j);
        Ainv(GetParam().Ainv[i][j], i, j);
      }
  }
};

TEST_P(HermCMatrixTest, InvertTest) {
  mpp_ba::SetTolerance(1e-14);
  HermCMatrix B(A);
  EXPECT_EQ(B.Invert(), Ainv);
  HermCMatrix C(Ainv);
  EXPECT_EQ(C.Invert(), A);
}

INSTANTIATE_TEST_CASE_P(
    BasicAlgebraTest, HermCMatrixTest,
    Values(TestCases<std::complex<double>>{3,
                                           {{1, -1, 2}, {-1, 0, -2}, {2, -2, 2}},
                                           {{-2, -1, 2}, {-1, -1, 0}, {1, 0, -0.5}}},
           TestCases<std::complex<double>>{3,
                                           {{1, 2, -1}, {2, 3, -2}, {-1, -2, 2}},
                                           {{-2, 2, 1}, {2, -1, 0}, {1, 0, 1}}}));

int main(int argc, char **argv) {
  MppTest mppTest = MppTestBuilder(argc, argv).WithPPM();
  return mppTest.RUN_ALL_MPP_TESTS();
}