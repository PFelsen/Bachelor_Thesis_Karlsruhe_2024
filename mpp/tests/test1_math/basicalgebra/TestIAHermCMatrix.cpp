#include "gtest/gtest.h"

#include "HermCMatrix.hpp"
#include "TestIABasicAlgebra.hpp"

class IAHermCMatrixTest : public IABasicAlgebraMatrixTest {
protected:
  IAHermCMatrix A, B;
  IASymRMatrix Sym;
  IACVector U;
  IARVector U_ia;
  IAInterval value_ia;
  double value_d;
  int value_i;

  int dim;
public:
  IAHermCMatrixTest() : IABasicAlgebraMatrixTest() {
    dim = cols = rows;
    A.resize(dim);
    B.resize(dim);
    Sym.resize(dim);
    U.resize(dim);
    U_ia.resize(dim);

    for (int n = 0; n < dim; ++n) {
      U[n] = RandomIACInterval();
      U_ia[n] = RandomIAInterval();
      A(RandomIAInterval(), n, n);
      B(RandomIAInterval(), n, n);
      Sym(RandomIAInterval(), n, n);
      for (int m = n + 1; m < dim; ++m) {
        A(RandomIACInterval(), n, m);
        B(RandomIACInterval(), n, m);
        Sym(RandomIAInterval(), n, m);
      }
    }
    value_ia = RandomNonZeroIAInterval();
    value_d = RandomNonZeroDouble();
    value_i = RandomNonZeroInt();
  }
};

TEST_P(IAHermCMatrixTest, ConstructorTest) {
  IAHermCMatrix W0(dim);
  IAHermCMatrix W1(A);
  IAHermCMatrix W2(mid(A));
  IAHermCMatrix W3(Sym);
  IAHermCMatrix W4(mid(Sym));

  for (int n = 0; n < dim; ++n)
    for (int m = 0; m < dim; ++m) {
      EXPECT_EQ(W0(n, m), IACInterval());
      EXPECT_EQ(W1(n, m), A(n, m));
      EXPECT_EQ(W2(n, m), mid(A(n, m)));
      EXPECT_EQ(W3(n, m), IACInterval(Sym(n, m)));
      EXPECT_EQ(W4(n, m), IACInterval(mid(Sym(n, m))));
    }
}

TEST_P(IAHermCMatrixTest, AssignmentTest) {
  IAHermCMatrix Z0 = A;
  EXPECT_EQ(Z0, A);
  IAHermCMatrix Z1;
  Z1 = mid(A);
  EXPECT_EQ(Z1, IAHermCMatrix(mid(A)));
  IAHermCMatrix Z2;
  Z2 = Sym;
  EXPECT_EQ(Z2, IAHermCMatrix(Sym));
  IAHermCMatrix Z3;
  Z3 = mid(Sym);
  EXPECT_EQ(Z3, IAHermCMatrix(mid(Sym)));
}

TEST_P(IAHermCMatrixTest, SummationTest) {
  IAHermCMatrix AB = A + B;
  IAHermCMatrix W1 = A + mid(B);
  IAHermCMatrix W2 = mid(A) + B;
  IAHermCMatrix W3 = A + Sym;
  IAHermCMatrix W4 = Sym + B;
  IAHermCMatrix W5 = A + mid(Sym);
  IAHermCMatrix W6 = mid(Sym) + B;
  for (int n = 0; n < dim; ++n)
    for (int m = 0; m < dim; ++m) {
      EXPECT_EQ(AB(n, m), A(n, m) + B(n, m));
      EXPECT_EQ(W1(n, m), A(n, m) + mid(B(n, m)));
      EXPECT_EQ(W2(n, m), mid(A(n, m)) + B(n, m));
      EXPECT_EQ(W3(n, m), A(n, m) + Sym(n, m));
      EXPECT_EQ(W4(n, m), Sym(n, m) + B(n, m));
      EXPECT_EQ(W5(n, m), A(n, m) + mid(Sym(n, m)));
      EXPECT_EQ(W6(n, m), mid(Sym(n, m)) + B(n, m));
    }
}

TEST_P(IAHermCMatrixTest, DifferenceTest) {
  IAHermCMatrix AB = A - B;
  IAHermCMatrix W1 = A - mid(B);
  IAHermCMatrix W2 = mid(A) - B;
  IAHermCMatrix W3 = A - Sym;
  IAHermCMatrix W4 = Sym - B;
  IAHermCMatrix W5 = A - mid(Sym);
  IAHermCMatrix W6 = mid(Sym) - B;
  for (int n = 0; n < dim; ++n)
    for (int m = 0; m < dim; ++m) {
      EXPECT_EQ(AB(n, m), A(n, m) - B(n, m));
      EXPECT_EQ(W1(n, m), A(n, m) - mid(B(n, m)));
      EXPECT_EQ(W2(n, m), mid(A(n, m)) - B(n, m));
      EXPECT_EQ(W3(n, m), A(n, m) - Sym(n, m));
      EXPECT_EQ(W4(n, m), Sym(n, m) - B(n, m));
      EXPECT_EQ(W5(n, m), A(n, m) - mid(Sym(n, m)));
      EXPECT_EQ(W6(n, m), mid(Sym(n, m)) - B(n, m));
    }
}

TEST_P(IAHermCMatrixTest, AdditiveInversTest) {
  IAHermCMatrix W = -A;
  for (int n = 0; n < dim; ++n)
    for (int m = 0; m < dim; ++m)
      EXPECT_EQ(W(n, m), -A(n, m));
}

TEST_P(IAHermCMatrixTest, ScalarMultiplicationTest) {
  IAHermCMatrix W1 = value_ia * A;
  IAHermCMatrix W2 = A * value_ia;
  IAHermCMatrix W3 = value_d * A;
  IAHermCMatrix W4 = A * value_d;
  IAHermCMatrix W5 = value_i * A;
  IAHermCMatrix W6 = A * value_i;
  for (int n = 0; n < dim; ++n) {
    EXPECT_EQ(W1[n], value_ia * A[n]);
    EXPECT_EQ(W2[n], value_ia * A[n]);
    EXPECT_EQ(W3[n], value_d * A[n]);
    EXPECT_EQ(W4[n], value_d * A[n]);
    EXPECT_EQ(W5[n], value_i * A[n]);
    EXPECT_EQ(W6[n], value_i * A[n]);
  }
}

TEST_P(IAHermCMatrixTest, ScalarDivisionTest) {
  IAHermCMatrix W1 = A / value_ia;
  IAHermCMatrix W2 = A / value_d;
  IAHermCMatrix W3 = A / value_i;
  for (int n = 0; n < dim; ++n) {
    EXPECT_EQ(W1[n], A[n] / value_ia);
    EXPECT_EQ(W2[n], A[n] / value_d);
    EXPECT_EQ(W3[n], A[n] / value_i);
  }
}

TEST_P(IAHermCMatrixTest, MatrixVectorProductTest) {
  IACVector W0(dim);
  IACVector W1(dim);
  IACVector W2(dim);
  IACVector W3(dim);
  IACVector W4(dim);
  IACVector W5(dim);
  for (int n = 0; n < dim; ++n)
    for (int k = 0; k < dim; ++k) {
      W0[n] += A(n, k) * U[k];
      W1[n] += A(n, k) * mid(U[k]);
      W2[n] += mid(A(n, k)) * U[k];
      W3[n] += A(n, k) * U_ia[k];
      W4[n] += A(n, k) * mid(U_ia[k]);
      W5[n] += IACInterval(mid(A(n, k))) * U_ia[k];
    }
  EXPECT_EQ(W0, A * U);
  EXPECT_EQ(W1, A * mid(U));
  EXPECT_EQ(W2, mid(A) * U);
  EXPECT_EQ(W3, A * U_ia);
  EXPECT_EQ(W4, A * mid(U_ia));
  EXPECT_EQ(W5, mid(A) * U_ia);
}

TEST_P(IAHermCMatrixTest, ConjTest) {
  IAHermCMatrix B = conj(A);
  for (int n = 0; n < rows; ++n)
    for (int m = 0; m < cols; ++m)
      EXPECT_EQ(B(n, m), conj(A(n, m)));
}

TEST_P(IAHermCMatrixTest, RealTest) {
  IASymRMatrix B = real(A);
  for (int n = 0; n < rows; ++n)
    for (int m = 0; m < cols; ++m)
      EXPECT_EQ(B(n, m), real(A(n, m)));
}

TEST_P(IAHermCMatrixTest, ImagTest) {
  IAAntisymRMatrix B = imag(A);
  for (int n = 0; n < rows; ++n)
    for (int m = 0; m < cols; ++m)
      EXPECT_EQ(B(n, m), imag(A(n, m)));
}

TEST_P(IAHermCMatrixTest, TransposeTest) {
  IAHermCMatrix B = transpose(A);
  for (int n = 0; n < cols; ++n)
    for (int m = 0; m < rows; ++m)
      EXPECT_EQ(B(n, m), A(m, n));
}

TEST_P(IAHermCMatrixTest, AdjointTest) {
  IAHermCMatrix B = adjoint(A);
  for (int n = 0; n < cols; ++n)
    for (int m = 0; m < rows; ++m)
      EXPECT_EQ(B(n, m), conj(A(m, n)));
}

TEST_P(IAHermCMatrixTest, EqualTest) {
  EXPECT_EQ(A, A);
  EXPECT_EQ(B, B);
}

TEST_P(IAHermCMatrixTest, NonEqualTest) {
  EXPECT_NE(A, B);
  EXPECT_NE(B, A);
}

TEST_P(IAHermCMatrixTest, SaveLoadTest) {
  Saver saver("SaveLoadTest");
  saver << A;
  saver.close();
  Loader loader("SaveLoadTest");
  IAHermCMatrix Z;
  loader >> Z;
  loader.close();
  EXPECT_EQ(A, Z);
}

TEST_P(IAHermCMatrixTest, MidTest) {
  HermCMatrix W = mid(A);
  for (int n = 0; n < dim; ++n)
    for (int m = 0; m < dim; ++m)
      EXPECT_EQ(W(n, m), mid(A(n, m)));
}

INSTANTIATE_TEST_CASE_P(IABasicAlgebraMatrixTest, IAHermCMatrixTest,
                        Values(IABasicAlgebraMatrixTestParameter{1, 1},
                               IABasicAlgebraMatrixTestParameter{2, 2},
                               IABasicAlgebraMatrixTestParameter{3, 3},
                               IABasicAlgebraMatrixTestParameter{4, 4},
                               IABasicAlgebraMatrixTestParameter{5, 5}));

int main(int argc, char **argv) {
  MppTest mppTest = MppTestBuilder(argc, argv).WithPPM();
  return mppTest.RUN_ALL_MPP_TESTS();
}
