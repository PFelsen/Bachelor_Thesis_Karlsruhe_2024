#include "gtest/gtest.h"

#include "CMatrix.hpp"
#include "TestIABasicAlgebra.hpp"

class IACMatrixTest : public IABasicAlgebraMatrixTest {
protected:
  IACMatrix A, B;
  IARMatrix A_ia;
  IAHermCMatrix Herm;
  IASymRMatrix Sym;
  IAAntisymRMatrix Antisym;
  IACVector U;
  IARVector U_ia;
  IACInterval value;
  IAInterval value_ia;
  std::complex<double> value_c;
  double value_d;
  int value_i;
public:
  IACMatrixTest() : IABasicAlgebraMatrixTest() {
    A.resize(rows, cols);
    B.resize(rows, cols);
    A_ia.resize(rows, cols);
    U.resize(cols);
    U_ia.resize(cols);

    for (int n = 0; n < cols; ++n) {
      U[n] = RandomIACInterval();
      U_ia[n] = RandomIAInterval();
    }
    for (int n = 0; n < rows; ++n) {
      for (int m = 0; m < cols; ++m) {
        A(RandomIACInterval(), n, m);
        B[n][m] = RandomIACInterval();
        A_ia[n][m] = RandomIAInterval();
      }
    }

    if (rows == cols) {
      Herm.resize(rows);
      Sym.resize(rows);
      Antisym.resize(rows);
      for (int n = 0; n < rows; ++n) {
        Herm(RandomIAInterval(), n, n);
        Sym(RandomIAInterval(), n, n);
        for (int m = n + 1; m < cols; ++m) {
          Herm(RandomIACInterval(), n, m);
          Sym(RandomIAInterval(), n, m);
          Antisym(RandomIAInterval(), n, m);
        }
      }
    }

    value = RandomNonZeroIACInterval();
    value_ia = RandomNonZeroIAInterval();
    value_c = RandomNonZeroComplex();
    value_d = RandomNonZeroDouble();
    value_i = RandomNonZeroInt();
  }
};

TEST_P(IACMatrixTest, ConstructorTest) {
  IACMatrix W0(rows, cols);
  IACMatrix W1(value, rows, cols);
  IACMatrix W2(value_ia, rows, cols);
  IACMatrix W3(value_c, rows, cols);
  IACMatrix W4(value_d, rows, cols);
  IACMatrix W5(value_i, rows, cols);
  IACMatrix W6(A);
  IACMatrix W7(mid(A));
  IACMatrix W8(A_ia);
  IACMatrix W9(mid(A_ia));
  IACMatrix W10(B);
  IACMatrix W11(Herm);
  IACMatrix W12(mid(Herm));
  IACMatrix W13(Sym);
  IACMatrix W14(mid(Sym));
  IACMatrix W15(Antisym);
  IACMatrix W16(mid(Antisym));
  for (int n = 0; n < rows; ++n)
    for (int m = 0; m < cols; ++m) {
      EXPECT_EQ(W0(n, m), IACInterval());
      EXPECT_EQ(W1(n, m), value);
      EXPECT_EQ(W2(n, m), IACInterval(value_ia));
      EXPECT_EQ(W3(n, m), IACInterval(value_c));
      EXPECT_EQ(W4(n, m), IACInterval(value_d));
      EXPECT_EQ(W5(n, m), IAInterval(value_i));
      EXPECT_EQ(W6(n, m), A(n, m));
      EXPECT_EQ(W7(n, m), IACInterval(mid(A(n, m))));
      EXPECT_EQ(W8(n, m), IACInterval(A_ia(n, m)));
      EXPECT_EQ(W9(n, m), IACInterval(mid(A_ia(n, m))));
      EXPECT_EQ(W10[n][m], B[n][m]);
      if (rows == cols) {
        EXPECT_EQ(W11(n, m), Herm(n, m));
        EXPECT_EQ(W12(n, m), IACInterval(mid(Herm(n, m))));
        EXPECT_EQ(W13(n, m), IACInterval(Sym(n, m)));
        EXPECT_EQ(W14(n, m), IACInterval(mid(Sym(n, m))));
        EXPECT_EQ(W15(n, m), IACInterval(Antisym(n, m)));
        EXPECT_EQ(W16(n, m), IACInterval(mid(Antisym(n, m))));
      }
    }

  if (rows == cols) {
    IACMatrix W(rows);
    for (int n = 0; n < rows; ++n)
      for (int m = 0; m < cols; ++m)
        EXPECT_EQ(W[n][m], IACInterval());
  }
}

TEST_P(IACMatrixTest, AssignmentTest) {
  IACMatrix Z0 = A;
  EXPECT_EQ(Z0, A);
  IACMatrix Z1;
  Z1 = mid(A);
  EXPECT_EQ(Z1, IACMatrix(mid(A)));
  Z0 = A_ia;
  EXPECT_EQ(Z0, A_ia);
  Z1 = mid(A_ia);
  EXPECT_EQ(Z1, IACMatrix(mid(A_ia)));
  if (rows == cols) {
    Z0 = Herm;
    EXPECT_EQ(Z0, Herm);
    Z1 = mid(Herm);
    EXPECT_EQ(Z1, IACMatrix(mid(Herm)));
    Z0 = Sym;
    EXPECT_EQ(Z0, Sym);
    Z1 = mid(Sym);
    EXPECT_EQ(Z1, IACMatrix(mid(Sym)));
    Z0 = Antisym;
    EXPECT_EQ(Z0, Antisym);
    Z1 = mid(Antisym);
    EXPECT_EQ(Z1, IACMatrix(mid(Antisym)));
  }
  IACMatrix W(rows, cols);
  W = value;
  EXPECT_EQ(W, IACMatrix(value, rows, cols));
  W = value_ia;
  EXPECT_EQ(W, IACMatrix(value_ia, rows, cols));
  W = value_c;
  EXPECT_EQ(W, IACMatrix(value_c, rows, cols));
  W = value_d;
  EXPECT_EQ(W, IACMatrix(value_d, rows, cols));
  W = value_i;
  EXPECT_EQ(W, IACMatrix(value_i, rows, cols));
}

TEST_P(IACMatrixTest, SummationTest) {
  IACMatrix AB = A + B;
  IACMatrix W1 = A + mid(B);
  IACMatrix W2 = mid(A) + B;
  IACMatrix W3 = A + A_ia;
  IACMatrix W4 = A_ia + B;
  IACMatrix W5 = A + mid(A_ia);
  IACMatrix W6 = mid(A_ia) + B;
  for (int n = 0; n < rows; ++n)
    for (int m = 0; m < cols; ++m) {
      EXPECT_EQ(AB(n, m), A(n, m) + B(n, m));
      EXPECT_EQ(W1(n, m), A(n, m) + mid(B(n, m)));
      EXPECT_EQ(W2(n, m), mid(A(n, m)) + B(n, m));
      EXPECT_EQ(W3(n, m), A(n, m) + A_ia(n, m));
      EXPECT_EQ(W4(n, m), A_ia(n, m) + B(n, m));
      EXPECT_EQ(W5(n, m), A(n, m) + mid(A_ia(n, m)));
      EXPECT_EQ(W6(n, m), mid(A_ia(n, m)) + B(n, m));
    }

  if (rows == cols) {
    IACMatrix W7 = A + Herm;
    IACMatrix W8 = Herm + B;
    IACMatrix W9 = A + mid(Herm);
    IACMatrix W10 = mid(Herm) + B;
    IACMatrix W11 = A + Sym;
    IACMatrix W12 = Sym + B;
    IACMatrix W13 = A + mid(Sym);
    IACMatrix W14 = mid(Sym) + B;
    IACMatrix W15 = A + Antisym;
    IACMatrix W16 = Antisym + B;
    IACMatrix W17 = A + mid(Antisym);
    IACMatrix W18 = mid(Antisym) + B;
    for (int n = 0; n < rows; ++n)
      for (int m = 0; m < cols; ++m) {
        EXPECT_EQ(W7(n, m), A(n, m) + Herm(n, m));
        EXPECT_EQ(W8(n, m), Herm(n, m) + B(n, m));
        EXPECT_EQ(W9(n, m), A(n, m) + mid(Herm(n, m)));
        EXPECT_EQ(W10(n, m), mid(Herm(n, m)) + B(n, m));
        EXPECT_EQ(W11(n, m), A(n, m) + Sym(n, m));
        EXPECT_EQ(W12(n, m), Sym(n, m) + B(n, m));
        EXPECT_EQ(W13(n, m), A(n, m) + mid(Sym(n, m)));
        EXPECT_EQ(W14(n, m), mid(Sym(n, m)) + B(n, m));
        EXPECT_EQ(W15(n, m), A(n, m) + Antisym(n, m));
        EXPECT_EQ(W16(n, m), Antisym(n, m) + B(n, m));
        EXPECT_EQ(W17(n, m), A(n, m) + mid(Antisym(n, m)));
        EXPECT_EQ(W18(n, m), mid(Antisym(n, m)) + B(n, m));
      }
  }
}

TEST_P(IACMatrixTest, DifferenceTest) {
  IACMatrix AB = A - B;
  IACMatrix W1 = A - mid(B);
  IACMatrix W2 = mid(A) - B;
  IACMatrix W3 = A - A_ia;
  IACMatrix W4 = A_ia - B;
  IACMatrix W5 = A - mid(A_ia);
  IACMatrix W6 = mid(A_ia) - B;
  for (int n = 0; n < rows; ++n)
    for (int m = 0; m < cols; ++m) {
      EXPECT_EQ(AB(n, m), A(n, m) - B(n, m));
      EXPECT_EQ(W1(n, m), A(n, m) - mid(B(n, m)));
      EXPECT_EQ(W2(n, m), mid(A(n, m)) - B(n, m));
      EXPECT_EQ(W3(n, m), A(n, m) - A_ia(n, m));
      EXPECT_EQ(W4(n, m), A_ia(n, m) - B(n, m));
      EXPECT_EQ(W5(n, m), A(n, m) - mid(A_ia(n, m)));
      EXPECT_EQ(W6(n, m), mid(A_ia(n, m)) - B(n, m));
    }

  if (rows == cols) {
    IACMatrix W7 = A - Herm;
    IACMatrix W8 = Herm - B;
    IACMatrix W9 = A - mid(Herm);
    IACMatrix W10 = mid(Herm) - B;
    IACMatrix W11 = A - Sym;
    IACMatrix W12 = Sym - B;
    IACMatrix W13 = A - mid(Sym);
    IACMatrix W14 = mid(Sym) - B;
    IACMatrix W15 = A - Antisym;
    IACMatrix W16 = Antisym - B;
    IACMatrix W17 = A - mid(Antisym);
    IACMatrix W18 = mid(Antisym) - B;
    for (int n = 0; n < rows; ++n)
      for (int m = 0; m < cols; ++m) {
        EXPECT_EQ(W7(n, m), A(n, m) - Herm(n, m));
        EXPECT_EQ(W8(n, m), Herm(n, m) - B(n, m));
        EXPECT_EQ(W9(n, m), A(n, m) - mid(Herm(n, m)));
        EXPECT_EQ(W10(n, m), mid(Herm(n, m)) - B(n, m));
        EXPECT_EQ(W11(n, m), A(n, m) - Sym(n, m));
        EXPECT_EQ(W12(n, m), Sym(n, m) - B(n, m));
        EXPECT_EQ(W13(n, m), A(n, m) - mid(Sym(n, m)));
        EXPECT_EQ(W14(n, m), mid(Sym(n, m)) - B(n, m));
        EXPECT_EQ(W15(n, m), A(n, m) - Antisym(n, m));
        EXPECT_EQ(W16(n, m), Antisym(n, m) - B(n, m));
        EXPECT_EQ(W17(n, m), A(n, m) - mid(Antisym(n, m)));
        EXPECT_EQ(W18(n, m), mid(Antisym(n, m)) - B(n, m));
      }
  }
}

TEST_P(IACMatrixTest, AdditiveInversTest) {
  IACMatrix W = -A;
  for (int n = 0; n < rows; ++n)
    for (int m = 0; m < cols; ++m)
      EXPECT_EQ(W(n, m), -A(n, m));
}

TEST_P(IACMatrixTest, ScalarMultiplicationTest) {
  IACMatrix W1 = value * A;
  IACMatrix W2 = A * value;
  IACMatrix W3 = value_ia * A;
  IACMatrix W4 = A * value_ia;
  IACMatrix W5 = value_c * A;
  IACMatrix W6 = A * value_c;
  IACMatrix W7 = value_d * A;
  IACMatrix W8 = A * value_d;
  IACMatrix W9 = value_i * A;
  IACMatrix W10 = A * value_i;
  for (int n = 0; n < rows; ++n)
    for (int m = 0; m < cols; ++m) {
      EXPECT_EQ(W1(n, m), value * A(n, m));
      EXPECT_EQ(W2(n, m), value * A(n, m));
      EXPECT_EQ(W3(n, m), value_ia * A(n, m));
      EXPECT_EQ(W4(n, m), value_ia * A(n, m));
      EXPECT_EQ(W5(n, m), value_c * A(n, m));
      EXPECT_EQ(W6(n, m), value_c * A(n, m));
      EXPECT_EQ(W7(n, m), value_d * A(n, m));
      EXPECT_EQ(W8(n, m), value_d * A(n, m));
      EXPECT_EQ(W9(n, m), value_i * A(n, m));
      EXPECT_EQ(W10(n, m), value_i * A(n, m));
    }
}

TEST_P(IACMatrixTest, ScalarDivisionTest) {
  IACMatrix W1 = A / value;
  IACMatrix W2 = A / value_ia;
  IACMatrix W3 = A / value_c;
  IACMatrix W4 = A / value_d;
  IACMatrix W5 = A / value_i;
  for (int n = 0; n < rows; ++n)
    for (int m = 0; m < cols; ++m) {
      EXPECT_EQ(W1(n, m), A(n, m) / value);
      EXPECT_EQ(W2(n, m), A(n, m) / value_ia);
      EXPECT_EQ(W3(n, m), A(n, m) / value_c);
      EXPECT_EQ(W4(n, m), A(n, m) / value_d);
      EXPECT_EQ(W5(n, m), A(n, m) / value_i);
    }
}

TEST_P(IACMatrixTest, MatrixVectorProductTest) {
  IACVector W0(rows);
  IACVector W1(rows);
  IACVector W2(rows);
  IACVector W3(rows);
  IACVector W4(rows);
  IACVector W5(rows);
  for (int n = 0; n < rows; ++n)
    for (int k = 0; k < cols; ++k) {
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

TEST_P(IACMatrixTest, MatrixMatrixProductTest) {
  IACMatrix transB = transpose(B);
  IACMatrix W0(rows, rows);
  IACMatrix W1(rows, rows);
  IACMatrix W2(rows, rows);
  IACMatrix W3(rows, rows);
  IACMatrix W4(rows, rows);
  for (int n = 0; n < rows; ++n)
    for (int m = 0; m < rows; ++m)
      for (int k = 0; k < cols; ++k) {
        W0[n][m] += A[n][k] * transB[k][m];
        W1[n][m] += A[n][k] * mid(transB[k][m]);
        W2[n][m] += mid(A[n][k]) * transB[k][m];
        W3[n][m] += A_ia[n][k] * transB[k][m];
        //        W4[n][m] += transpose(A_ia[n][k]) * transB[k][m];
      }
  EXPECT_EQ(W0, A * transpose(B));
  EXPECT_EQ(W1, A * mid(transpose(B)));
  EXPECT_EQ(W2, mid(A) * transpose(B));
  EXPECT_EQ(W3, A_ia * transpose(B));
  //  EXPECT_EQ(W4, mid(A_ia) * transpose(B));

  if (rows == cols) {
    IACMatrix W0(rows);
    IACMatrix W1(rows);
    IACMatrix W2(rows);
    IACMatrix W3(rows);
    IACMatrix V0(rows);
    IACMatrix V1(rows);
    IACMatrix V2(rows);
    IACMatrix V3(rows);
    IACMatrix X0(rows);
    IACMatrix X1(rows);
    IACMatrix X2(rows);
    IACMatrix X3(rows);
    for (int n = 0; n < rows; ++n)
      for (int m = 0; m < rows; ++m)
        for (int k = 0; k < cols; ++k) {
          W0[n][m] += A[n][k] * Sym(k, m);
          W1[n][m] += Sym(n, k) * B[k][m];
          W2[n][m] += A[n][k] * mid(Sym(k, m));
          W3[n][m] += IACInterval(mid(Sym(n, k))) * B[k][m];
          V0[n][m] += A[n][k] * Antisym(k, m);
          V1[n][m] += Antisym(n, k) * B[k][m];
          V2[n][m] += A[n][k] * mid(Antisym(k, m));
          V3[n][m] += mid(Antisym(n, k)) * B[k][m];
          X0[n][m] += A[n][k] * Herm(k, m);
          X1[n][m] += Herm(n, k) * B[k][m];
          X2[n][m] += A[n][k] * mid(Herm(k, m));
          X3[n][m] += mid(Herm(n, k)) * B[k][m];
        }
    EXPECT_EQ(W0, A * Sym);
    EXPECT_EQ(W1, Sym * B);
    EXPECT_EQ(W2, A * mid(Sym));
    EXPECT_EQ(W3, mid(Sym) * B);
    EXPECT_EQ(V0, A * Antisym);
    EXPECT_EQ(V1, Antisym * B);
    EXPECT_EQ(V2, A * mid(Antisym));
    EXPECT_EQ(V3, mid(Antisym) * B);
    EXPECT_EQ(X0, A * Herm);
    EXPECT_EQ(X1, Herm * B);
    EXPECT_EQ(X2, A * mid(Herm));
    EXPECT_EQ(X3, mid(Herm) * B);
  }
}

TEST_P(IACMatrixTest, ConjTest) {
  IACMatrix B = conj(A);
  for (int n = 0; n < rows; ++n)
    for (int m = 0; m < cols; ++m)
      EXPECT_EQ(B(n, m), conj(A(n, m)));
}

TEST_P(IACMatrixTest, RealTest) {
  IARMatrix B = real(A);
  for (int n = 0; n < rows; ++n)
    for (int m = 0; m < cols; ++m)
      EXPECT_EQ(B(n, m), real(A(n, m)));
}

TEST_P(IACMatrixTest, ImagTest) {
  IARMatrix B = imag(A);
  for (int n = 0; n < rows; ++n)
    for (int m = 0; m < cols; ++m)
      EXPECT_EQ(B(n, m), imag(A(n, m)));
}

TEST_P(IACMatrixTest, TransposeTest) {
  IACMatrix B = transpose(A);
  for (int n = 0; n < cols; ++n)
    for (int m = 0; m < rows; ++m)
      EXPECT_EQ(B(n, m), A(m, n));
}

TEST_P(IACMatrixTest, AdjointTest) {
  IACMatrix B = adjoint(A);
  for (int n = 0; n < cols; ++n)
    for (int m = 0; m < rows; ++m)
      EXPECT_EQ(B(n, m), conj(A(m, n)));
}

TEST_P(IACMatrixTest, DiagonalTest) {
  IACVector W = A.diag();
  for (int n = 0; n < std::min(rows, cols); ++n)
    EXPECT_EQ(W[n], A[n][n]);

  if (rows = cols) {
    IACMatrix W0(A);
    W0.diag(U);
    IACMatrix W1(A);
    W1.diag(mid(U));
    IACMatrix W2(A);
    W2.diag(U_ia);
    IACMatrix W3(A);
    W3.diag(mid(U_ia));
    for (int n = 0; n < U.size(); ++n)
      for (int m = 0; m < U.size(); ++m) {
        if (n == m) {
          EXPECT_EQ(W0(n, m), U[n]);
          EXPECT_EQ(W1(n, m), mid(U[n]));
          EXPECT_EQ(W2(n, m), IACInterval(U_ia[n]));
          EXPECT_EQ(W3(n, m), IACInterval(mid(U_ia[n])));
        } else {
          EXPECT_EQ(W0(n, m), IACInterval());
          EXPECT_EQ(W1(n, m), IACInterval());
          EXPECT_EQ(W2(n, m), IACInterval());
          EXPECT_EQ(W3(n, m), IACInterval());
        }
      }
  }
}

TEST_P(IACMatrixTest, EqualTest) {
  EXPECT_EQ(A, A);
  EXPECT_EQ(B, B);
  IACMatrix W1(Herm);
  EXPECT_EQ(W1, Herm);
  IACMatrix W2(Sym);
  EXPECT_EQ(W2, Sym);
  IACMatrix W3(Antisym);
  EXPECT_EQ(W3, Antisym);
  IACMatrix W4(A_ia);
  EXPECT_EQ(W4, A_ia);
}

TEST_P(IACMatrixTest, NonEqualTest) {
  EXPECT_NE(A, B);
  EXPECT_NE(B, A);
}

TEST_P(IACMatrixTest, SaveLoadTest) {
  Saver saver("SaveLoadTest");
  saver << A;
  saver.close();
  Loader loader("SaveLoadTest");
  IACMatrix Z;
  loader >> Z;
  loader.close();
  EXPECT_EQ(A, Z);
}

TEST_P(IACMatrixTest, MidTest) {
  CMatrix W = mid(A);
  for (int n = 0; n < rows; ++n)
    for (int m = 0; m < cols; ++m)
      EXPECT_EQ(W(n, m), mid(A(n, m)));
}

INSTANTIATE_TEST_CASE_P(
    IABasicAlgebraMatrixTest, IACMatrixTest,
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
