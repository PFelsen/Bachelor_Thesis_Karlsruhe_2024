#include "gtest/gtest.h"

#include "RMatrix.hpp"
#include "TestIABasicAlgebra.hpp"

class IARMatrixTest : public IABasicAlgebraMatrixTest {
protected:
  IARMatrix A, B;
  IASymRMatrix Sym;
  IAAntisymRMatrix Antisym;
  IARVector U;
  IAInterval value;
  double value_d;
  int value_i;
public:
  IARMatrixTest() : IABasicAlgebraMatrixTest() {
    A.resize(rows, cols);
    B.resize(rows, cols);
    U.resize(cols);

    for (int n = 0; n < cols; ++n) {
      U[n] = RandomIAInterval();
    }
    for (int n = 0; n < rows; ++n) {
      for (int m = 0; m < cols; ++m) {
        A(RandomIAInterval(), n, m);
        B[n][m] = RandomIAInterval();
      }
    }

    if (rows == cols) {
      Sym.resize(rows);
      Antisym.resize(rows);
      for (int n = 0; n < rows; ++n) {
        for (int m = n; m < cols; ++m) {
          Sym(RandomIAInterval(), n, m);
          Antisym(RandomIAInterval(), n, m);
        }
      }
    }

    value = RandomNonZeroIAInterval();
    value_d = RandomNonZeroDouble();
    value_i = RandomNonZeroInt();
  }
};

TEST_P(IARMatrixTest, ConstructorTest) {
  IARMatrix W0(rows, cols);
  IARMatrix W1(value, rows, cols);
  IARMatrix W2(value_d, rows, cols);
  IARMatrix W3(value_i, rows, cols);
  IARMatrix W4(A);
  IARMatrix W5(mid(A));
  IARMatrix W6(B);
  IARMatrix W7(Sym);
  IARMatrix W8(mid(Sym));
  IARMatrix W9(Antisym);
  IARMatrix W10(mid(Antisym));
  for (int n = 0; n < rows; ++n)
    for (int m = 0; m < cols; ++m) {
      EXPECT_EQ(W0(n, m), IAInterval());
      EXPECT_EQ(W1(n, m), value);
      EXPECT_EQ(W2(n, m), IAInterval(value_d));
      EXPECT_EQ(W3(n, m), IAInterval(value_i));
      EXPECT_EQ(W4(n, m), A(n, m));
      EXPECT_EQ(W5(n, m), mid(A(n, m)));
      EXPECT_EQ(W6[n][m], B[n][m]);
      if (rows == cols) {
        EXPECT_EQ(W7(n, m), Sym(n, m));
        EXPECT_EQ(W8(n, m), IAInterval(mid(Sym(n, m))));
        EXPECT_EQ(W9(n, m), Antisym(n, m));
        EXPECT_EQ(W10(n, m), IAInterval(mid(Antisym(n, m))));
      }
    }

  if (rows == cols) {
    IARMatrix W(rows);
    for (int n = 0; n < rows; ++n)
      for (int m = 0; m < cols; ++m)
        EXPECT_EQ(W[n][m], IAInterval());
  }
}

TEST_P(IARMatrixTest, AssignmentTest) {
  IARMatrix Z0 = A;
  EXPECT_EQ(Z0, A);
  IARMatrix Z1;
  Z1 = mid(A);
  EXPECT_EQ(Z1, IARMatrix(mid(A)));
  if (rows == cols) {
    Z0 = Sym;
    EXPECT_EQ(Z0, Sym);
    Z0 = Antisym;
    EXPECT_EQ(Z0, Antisym);
    Z1 = mid(Sym);
    EXPECT_EQ(Z1, IARMatrix(mid(Sym)));
    Z1 = mid(Antisym);
    EXPECT_EQ(Z1, IARMatrix(mid(Antisym)));
  }
  IARMatrix W(rows, cols);
  W = value;
  EXPECT_EQ(W, IARMatrix(value, rows, cols));
  W = value_d;
  EXPECT_EQ(W, IARMatrix(value_d, rows, cols));
  W = value_i;
  EXPECT_EQ(W, IARMatrix(value_i, rows, cols));
}

TEST_P(IARMatrixTest, SummationTest) {
  IARMatrix AB = A + B;
  IARMatrix W1 = A + mid(B);
  IARMatrix W2 = mid(A) + B;
  for (int n = 0; n < rows; ++n)
    for (int m = 0; m < cols; ++m) {
      EXPECT_EQ(AB(n, m), A(n, m) + B(n, m));
      EXPECT_EQ(W1(n, m), A(n, m) + mid(B(n, m)));
      EXPECT_EQ(W2(n, m), mid(A(n, m)) + B(n, m));
    }

  if (rows == cols) {
    IARMatrix W1 = A + Sym;
    IARMatrix W2 = A + mid(Sym);
    IARMatrix W3 = Sym + B;
    IARMatrix W4 = mid(Sym) + B;
    IARMatrix W5 = A + Antisym;
    IARMatrix W6 = A + mid(Antisym);
    IARMatrix W7 = Antisym + B;
    IARMatrix W8 = mid(Antisym) + B;
    for (int n = 0; n < rows; ++n)
      for (int m = 0; m < cols; ++m) {
        EXPECT_EQ(W1(n, m), A(n, m) + Sym(n, m));
        EXPECT_EQ(W2(n, m), A(n, m) + mid(Sym(n, m)));
        EXPECT_EQ(W3(n, m), Sym(n, m) + B(n, m));
        EXPECT_EQ(W4(n, m), mid(Sym(n, m)) + B(n, m));
        EXPECT_EQ(W5(n, m), A(n, m) + Antisym(n, m));
        EXPECT_EQ(W6(n, m), A(n, m) + mid(Antisym(n, m)));
        EXPECT_EQ(W7(n, m), Antisym(n, m) + B(n, m));
        EXPECT_EQ(W8(n, m), mid(Antisym(n, m)) + B(n, m));
      }
  }
}

TEST_P(IARMatrixTest, DifferenceTest) {
  IARMatrix AB = A - B;
  IARMatrix W1 = A - mid(B);
  IARMatrix W2 = mid(A) - B;
  for (int n = 0; n < rows; ++n)
    for (int m = 0; m < cols; ++m) {
      EXPECT_EQ(AB(n, m), A(n, m) - B(n, m));
      EXPECT_EQ(W1(n, m), A(n, m) - mid(B(n, m)));
      EXPECT_EQ(W2(n, m), mid(A(n, m)) - B(n, m));
    }

  if (rows == cols) {
    IARMatrix W1 = A - Sym;
    IARMatrix W2 = A - mid(Sym);
    IARMatrix W3 = Sym - B;
    IARMatrix W4 = mid(Sym) - B;
    IARMatrix W5 = A - Antisym;
    IARMatrix W6 = A - mid(Antisym);
    IARMatrix W7 = Antisym - B;
    IARMatrix W8 = mid(Antisym) - B;
    for (int n = 0; n < rows; ++n)
      for (int m = 0; m < cols; ++m) {
        EXPECT_EQ(W1(n, m), A(n, m) - Sym(n, m));
        EXPECT_EQ(W2(n, m), A(n, m) - mid(Sym(n, m)));
        EXPECT_EQ(W3(n, m), Sym(n, m) - B(n, m));
        EXPECT_EQ(W4(n, m), mid(Sym(n, m)) - B(n, m));
        EXPECT_EQ(W5(n, m), A(n, m) - Antisym(n, m));
        EXPECT_EQ(W6(n, m), A(n, m) - mid(Antisym(n, m)));
        EXPECT_EQ(W7(n, m), Antisym(n, m) - B(n, m));
        EXPECT_EQ(W8(n, m), mid(Antisym(n, m)) - B(n, m));
      }
  }
}

TEST_P(IARMatrixTest, AdditiveInversTest) {
  IARMatrix W = -A;
  for (int n = 0; n < rows; ++n)
    for (int m = 0; m < cols; ++m)
      EXPECT_EQ(W(n, m), -A(n, m));
}

TEST_P(IARMatrixTest, ScalarMultiplicationTest) {
  IARMatrix W1 = value * A;
  IARMatrix W2 = A * value;
  IARMatrix W3 = value_d * A;
  IARMatrix W4 = A * value_d;
  IARMatrix W5 = value_i * A;
  IARMatrix W6 = A * value_i;
  for (int n = 0; n < rows; ++n)
    for (int m = 0; m < cols; ++m) {
      EXPECT_EQ(W1(n, m), value * A(n, m));
      EXPECT_EQ(W2(n, m), value * A(n, m));
      EXPECT_EQ(W3(n, m), value_d * A(n, m));
      EXPECT_EQ(W4(n, m), value_d * A(n, m));
      EXPECT_EQ(W5(n, m), value_i * A(n, m));
      EXPECT_EQ(W6(n, m), value_i * A(n, m));
    }
}

TEST_P(IARMatrixTest, ScalarDivisionTest) {
  IARMatrix W1 = A / value;
  IARMatrix W2 = A / value_d;
  IARMatrix W3 = A / value_i;
  for (int n = 0; n < rows; ++n)
    for (int m = 0; m < cols; ++m) {
      EXPECT_EQ(W1(n, m), A(n, m) / value);
      EXPECT_EQ(W2(n, m), A(n, m) / value_d);
      EXPECT_EQ(W3(n, m), A(n, m) / value_i);
    }
}

TEST_P(IARMatrixTest, MatrixVectorProductTest) {
  IARVector W0(rows);
  IARVector W1(rows);
  IARVector W2(rows);
  for (int n = 0; n < rows; ++n)
    for (int k = 0; k < cols; ++k) {
      W0[n] += A(n, k) * U[k];
      W1[n] += A(n, k) * mid(U[k]);
      W2[n] += mid(A(n, k)) * U[k];
    }
  EXPECT_EQ(W0, A * U);
  EXPECT_EQ(W1, A * mid(U));
  EXPECT_EQ(W2, mid(A) * U);
}

TEST_P(IARMatrixTest, MatrixMatrixProductTest) {
  IARMatrix transB = transpose(B);
  IARMatrix W0(rows, rows);
  IARMatrix W1(rows, rows);
  IARMatrix W2(rows, rows);
  for (int n = 0; n < rows; ++n)
    for (int m = 0; m < rows; ++m)
      for (int k = 0; k < cols; ++k) {
        W0[n][m] += A[n][k] * transB[k][m];
        W1[n][m] += A[n][k] * mid(transB[k][m]);
        W2[n][m] += mid(A[n][k]) * transB[k][m];
      }
  EXPECT_EQ(W0, A * transpose(B));
  EXPECT_EQ(W1, A * mid(transpose(B)));
  EXPECT_EQ(W2, mid(A) * transpose(B));

  if (rows == cols) {
    IARMatrix W0(rows);
    IARMatrix W1(rows);
    IARMatrix W2(rows);
    IARMatrix W3(rows);
    IARMatrix W4(rows);
    IARMatrix W5(rows);
    IARMatrix V0(rows);
    IARMatrix V1(rows);
    IARMatrix V2(rows);
    IARMatrix V3(rows);
    IARMatrix V4(rows);
    IARMatrix V5(rows);
    for (int n = 0; n < rows; ++n)
      for (int m = 0; m < rows; ++m)
        for (int k = 0; k < cols; ++k) {
          W0[n][m] += A[n][k] * Sym(k, m);
          W1[n][m] += Sym(n, k) * B[k][m];
          W2[n][m] += A[n][k] * mid(Sym(k, m));
          W3[n][m] += mid(Sym(n, k)) * B[k][m];
          W4[n][m] += mid(A[n][k]) * Sym(k, m);
          W5[n][m] += Sym(n, k) * mid(B[k][m]);
          V0[n][m] += A[n][k] * Antisym(k, m);
          V1[n][m] += Antisym(n, k) * B[k][m];
          V2[n][m] += A[n][k] * mid(Antisym(k, m));
          V3[n][m] += mid(Antisym(n, k)) * B[k][m];
          V4[n][m] += mid(A[n][k]) * Antisym(k, m);
          V5[n][m] += Antisym(n, k) * mid(B[k][m]);
        }
    EXPECT_EQ(W0, A * Sym);
    EXPECT_EQ(W1, Sym * B);
    EXPECT_EQ(W2, A * mid(Sym));
    EXPECT_EQ(W3, mid(Sym) * B);
    EXPECT_EQ(W4, mid(A) * Sym);
    EXPECT_EQ(W5, Sym * mid(B));
    EXPECT_EQ(V0, A * Antisym);
    EXPECT_EQ(V1, Antisym * B);
    EXPECT_EQ(V2, A * mid(Antisym));
    EXPECT_EQ(V3, mid(Antisym) * B);
    EXPECT_EQ(V4, mid(A) * Antisym);
    EXPECT_EQ(V5, Antisym * mid(B));
  }
}

TEST_P(IARMatrixTest, TransposeTest) {
  IARMatrix B = transpose(A);
  for (int n = 0; n < cols; ++n)
    for (int m = 0; m < rows; ++m)
      EXPECT_EQ(B[n][m], A[m][n]);
}

TEST_P(IARMatrixTest, DiagonalTest) {
  IARVector W = A.diag();
  for (int n = 0; n < std::min(rows, cols); ++n)
    EXPECT_EQ(W[n], A[n][n]);

  if (rows = cols) {
    IARMatrix W0(A);
    W0.diag(U);
    IARMatrix W1(A);
    W1.diag(mid(U));
    for (int n = 0; n < U.size(); ++n)
      for (int m = 0; m < U.size(); ++m) {
        if (n == m) {
          EXPECT_EQ(W0(n, m), U[n]);
          EXPECT_EQ(W1(n, m), mid(U[n]));
        } else {
          EXPECT_EQ(W0(n, m), IAInterval());
          EXPECT_EQ(W1(n, m), IAInterval());
        }
      }
  }
}

TEST_P(IARMatrixTest, EqualTest) {
  EXPECT_EQ(A, A);
  EXPECT_EQ(B, B);
  IARMatrix W1(Sym);
  EXPECT_EQ(W1, Sym);
  IARMatrix W2(Antisym);
  EXPECT_EQ(W2, Antisym);
}

TEST_P(IARMatrixTest, NonEqualTest) {
  EXPECT_NE(A, B);
  EXPECT_NE(B, A);
}

TEST_P(IARMatrixTest, SaveLoadTest) {
  Saver saver("SaveLoadTest");
  saver << A;
  saver.close();
  Loader loader("SaveLoadTest");
  IARMatrix Z;
  loader >> Z;
  loader.close();
  EXPECT_EQ(A, Z);
}

TEST_P(IARMatrixTest, InfTest) {
  RMatrix W = inf(A);
  for (int n = 0; n < rows; ++n)
    for (int m = 0; m < cols; ++m)
      EXPECT_EQ(W(n, m), inf(A(n, m)));
}

TEST_P(IARMatrixTest, SupTest) {
  RMatrix W = sup(A);
  for (int n = 0; n < rows; ++n)
    for (int m = 0; m < cols; ++m)
      EXPECT_EQ(W(n, m), sup(A(n, m)));
}

TEST_P(IARMatrixTest, MidTest) {
  RMatrix W = mid(A);
  for (int n = 0; n < rows; ++n)
    for (int m = 0; m < cols; ++m)
      EXPECT_EQ(W(n, m), mid(A(n, m)));
}

INSTANTIATE_TEST_CASE_P(
    IABasicAlgebraMatrixTest, IARMatrixTest,
    Values(IABasicAlgebraMatrixTestParameter{1, 1}, IABasicAlgebraMatrixTestParameter{1, 2},
           IABasicAlgebraMatrixTestParameter{2, 1}, IABasicAlgebraMatrixTestParameter{1, 3},
           IABasicAlgebraMatrixTestParameter{2, 2}, IABasicAlgebraMatrixTestParameter{3, 1},
           IABasicAlgebraMatrixTestParameter{1, 4}, IABasicAlgebraMatrixTestParameter{2, 3},
           IABasicAlgebraMatrixTestParameter{3, 2}, IABasicAlgebraMatrixTestParameter{4, 1},
           IABasicAlgebraMatrixTestParameter{1, 5}, IABasicAlgebraMatrixTestParameter{2, 4},
           IABasicAlgebraMatrixTestParameter{3, 3}, IABasicAlgebraMatrixTestParameter{4, 2},
           IABasicAlgebraMatrixTestParameter{5, 1}));

int main(int argc, char **argv) {
  MppTest mppTest = MppTestBuilder(argc, argv).WithPPM();
  return mppTest.RUN_ALL_MPP_TESTS();
}
