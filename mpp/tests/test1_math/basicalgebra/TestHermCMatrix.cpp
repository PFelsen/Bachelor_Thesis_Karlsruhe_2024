#include "gtest/gtest.h"

#include "TestBasicAlgebra.hpp"
#include "basicalgebra/HermCMatrix.hpp"

class HermCMatrixTest : public BasicAlgebraMatrixTest {
protected:
  HermCMatrix A, B;
  SymRMatrix Sym;
  CVector U;
  RVector U_d;
  double value_d;
  int value_i;

  int dim;
public:
  HermCMatrixTest() : BasicAlgebraMatrixTest() {
    dim = cols = rows;
    A.resize(dim);
    B.resize(dim);
    Sym.resize(dim);
    U.resize(dim);
    U_d.resize(dim);

    for (int n = 0; n < dim; ++n) {
      U[n] = RandomComplex();
      U_d[n] = RandomDouble();
      A(RandomDouble(), n, n);
      B(RandomDouble(), n, n);
      Sym(RandomDouble(), n, n);
      for (int m = n + 1; m < dim; ++m) {
        A(RandomComplex(), n, m);
        B(RandomComplex(), n, m);
        Sym(RandomDouble(), n, m);
      }
    }
    value_d = RandomNonZeroDouble();
    value_i = RandomNonZeroInt();
  }
};

TEST_P(HermCMatrixTest, ConstructorTest) {
  HermCMatrix W0(dim);
  HermCMatrix W1(A);
  HermCMatrix W2(Sym);
  for (int n = 0; n < dim; ++n)
    for (int m = 0; m < dim; ++m) {
      EXPECT_COMPLEX_EQ(W0(n, m), std::complex<double>(0.0));
      EXPECT_COMPLEX_EQ(W1(n, m), A(n, m));
      EXPECT_COMPLEX_EQ(W2(n, m), std::complex<double>(Sym(n, m)));
    }
}

TEST_P(HermCMatrixTest, AssignmentTest) {
  HermCMatrix Z0 = A;
  EXPECT_EQ(Z0, A);
  HermCMatrix Z1;
  Z1 = Sym;
  EXPECT_EQ(Z1, HermCMatrix(Sym));
}

TEST_P(HermCMatrixTest, SummationTest) {
  HermCMatrix AB = A + B;
  HermCMatrix W1 = A + Sym;
  HermCMatrix W2 = Sym + B;
  for (int n = 0; n < dim; ++n)
    for (int m = 0; m < dim; ++m) {
      EXPECT_COMPLEX_EQ(AB(n, m), A(n, m) + B(n, m));
      EXPECT_COMPLEX_EQ(W1(n, m), A(n, m) + Sym(n, m));
      EXPECT_COMPLEX_EQ(W2(n, m), Sym(n, m) + B(n, m));
    }
}

TEST_P(HermCMatrixTest, DifferenceTest) {
  HermCMatrix AB = A - B;
  HermCMatrix W1 = A - Sym;
  HermCMatrix W2 = Sym - B;
  for (int n = 0; n < dim; ++n)
    for (int m = 0; m < dim; ++m) {
      EXPECT_COMPLEX_EQ(AB(n, m), A(n, m) - B(n, m));
      EXPECT_COMPLEX_EQ(W1(n, m), A(n, m) - Sym(n, m));
      EXPECT_COMPLEX_EQ(W2(n, m), Sym(n, m) - B(n, m));
    }
}

TEST_P(HermCMatrixTest, AdditiveInversTest) {
  HermCMatrix W = -A;
  for (int n = 0; n < dim; ++n)
    for (int m = 0; m < dim; ++m)
      EXPECT_COMPLEX_EQ(W(n, m), -A(n, m));
}

TEST_P(HermCMatrixTest, ScalarMultiplicationTest) {
  HermCMatrix W1 = value_d * A;
  HermCMatrix W2 = A * value_d;
  HermCMatrix W3 = value_i * A;
  HermCMatrix W4 = A * value_i;
  for (int n = 0; n < dim; ++n) {
    EXPECT_COMPLEX_EQ(W1[n], value_d * A[n]);
    EXPECT_COMPLEX_EQ(W2[n], value_d * A[n]);
    EXPECT_COMPLEX_EQ(W3[n], value_i * A[n]);
    EXPECT_COMPLEX_EQ(W4[n], value_i * A[n]);
  }
}

TEST_P(HermCMatrixTest, ScalarDivisionTest) {
  HermCMatrix W1 = A / value_d;
  HermCMatrix W2 = A / value_i;
  for (int n = 0; n < dim; ++n) {
    EXPECT_COMPLEX_EQ(W1[n], A[n] / value_d);
    EXPECT_COMPLEX_EQ(W2[n], A[n] / value_i);
  }
}

TEST_P(HermCMatrixTest, MatrixVectorProductTest) {
  CVector W0(dim);
  CVector W1(dim);
  for (int n = 0; n < dim; ++n)
    for (int k = 0; k < dim; ++k) {
      W0[n] += A(n, k) * U[k];
      W1[n] += A(n, k) * U_d[k];
    }
  EXPECT_EQ(W0, A * U);
  EXPECT_EQ(W1, A * U_d);
}

TEST_P(HermCMatrixTest, ConjTest) {
  HermCMatrix B = conj(A);
  for (int n = 0; n < rows; ++n)
    for (int m = 0; m < cols; ++m)
      EXPECT_COMPLEX_EQ(B(n, m), conj(A(n, m)));
}

TEST_P(HermCMatrixTest, RealTest) {
  SymRMatrix B = real(A);
  for (int n = 0; n < rows; ++n)
    for (int m = 0; m < cols; ++m)
      EXPECT_DOUBLE_EQ(B(n, m), real(A(n, m)));
}

TEST_P(HermCMatrixTest, ImagTest) {
  AntisymRMatrix B = imag(A);
  for (int n = 0; n < rows; ++n)
    for (int m = 0; m < cols; ++m)
      EXPECT_DOUBLE_EQ(B(n, m), imag(A(n, m)));
}

TEST_P(HermCMatrixTest, TransposeTest) {
  HermCMatrix B = transpose(A);
  for (int n = 0; n < cols; ++n)
    for (int m = 0; m < rows; ++m)
      EXPECT_COMPLEX_EQ(B(n, m), A(m, n));
}

TEST_P(HermCMatrixTest, AdjointTest) {
  HermCMatrix B = adjoint(A);
  for (int n = 0; n < cols; ++n)
    for (int m = 0; m < rows; ++m)
      EXPECT_COMPLEX_EQ(B(n, m), conj(A(m, n)));
}

TEST_P(HermCMatrixTest, EqualTest) {
  EXPECT_EQ(A, A);
  EXPECT_EQ(B, B);
}

TEST_P(HermCMatrixTest, NonEqualTest) {
  EXPECT_NE(A, B);
  EXPECT_NE(B, A);
}

TEST_P(HermCMatrixTest, SaveLoadTest) {
  Saver saver("SaveLoadTest");
  saver << A;
  saver.close();
  Loader loader("SaveLoadTest");
  HermCMatrix Z;
  loader >> Z;
  loader.close();
  EXPECT_EQ(A, Z);
}

INSTANTIATE_TEST_CASE_P(BasicAlgebraMatrixTest, HermCMatrixTest,
                        Values(BasicAlgebraMatrixTestParameter{1, 1},
                               BasicAlgebraMatrixTestParameter{2, 2},
                               BasicAlgebraMatrixTestParameter{3, 3},
                               BasicAlgebraMatrixTestParameter{4, 4},
                               BasicAlgebraMatrixTestParameter{5, 5}));

int main(int argc, char **argv) {
  MppTest mppTest = MppTestBuilder(argc, argv).WithPPM();
  return mppTest.RUN_ALL_MPP_TESTS();
}
