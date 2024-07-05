#include "gtest/gtest.h"

#include "TestBasicAlgebra.hpp"
#include "basicalgebra/RMatrix.hpp"

class RMatrixTest : public BasicAlgebraMatrixTest {
protected:
  RMatrix A, B;
  SymRMatrix Sym;
  AntisymRMatrix Antisym;
  RVector U;
  double value;
  int value_i;
public:
  RMatrixTest() : BasicAlgebraMatrixTest() {
    A.resize(rows, cols);
    B.resize(rows, cols);
    U.resize(cols);

    for (int n = 0; n < cols; ++n) {
      U[n] = RandomDouble();
    }
    for (int n = 0; n < rows; ++n) {
      for (int m = 0; m < cols; ++m) {
        A(RandomDouble(), n, m);
        B[n][m] = RandomDouble();
      }
    }

    if (rows == cols) {
      Sym.resize(rows);
      Antisym.resize(rows);
      for (int n = 0; n < rows; ++n) {
        for (int m = n; m < cols; ++m) {
          Sym(RandomDouble(), n, m);
          Antisym(RandomDouble(), n, m);
        }
      }
    }

    value = RandomNonZeroDouble();
    value_i = RandomNonZeroInt();
  }
};

TEST_P(RMatrixTest, ConstructorTest) {
  RMatrix W0(rows, cols);
  RMatrix W1(value, rows, cols);
  RMatrix W2(value_i, rows, cols);
  RMatrix W3(A);
  RMatrix W4(B);
  RMatrix W5(Sym);
  RMatrix W6(Antisym);
  for (int n = 0; n < rows; ++n)
    for (int m = 0; m < cols; ++m) {
      EXPECT_DOUBLE_EQ(W0(n, m), 0.0);
      EXPECT_DOUBLE_EQ(W1(n, m), value);
      EXPECT_DOUBLE_EQ(W2(n, m), value_i);
      EXPECT_DOUBLE_EQ(W3(n, m), A(n, m));
      EXPECT_DOUBLE_EQ(W4[n][m], B[n][m]);
      if (rows == cols) {
        EXPECT_DOUBLE_EQ(W5(n, m), Sym(n, m));
        EXPECT_DOUBLE_EQ(W6(n, m), Antisym(n, m));
      }
    }

  if (rows == cols) {
    RMatrix W(rows);
    for (int n = 0; n < rows; ++n)
      for (int m = 0; m < cols; ++m)
        EXPECT_DOUBLE_EQ(W[n][m], 0.0);
  }
}

TEST_P(RMatrixTest, AssignmentTest) {
  RMatrix Z = A;
  EXPECT_EQ(Z, A);
  if (rows == cols) {
    Z = Sym;
    EXPECT_EQ(Z, Sym);
    Z = Antisym;
    EXPECT_EQ(Z, Antisym);
  }
  RMatrix W(rows, cols);
  W = value;
  EXPECT_EQ(W, RMatrix(value, rows, cols));
  W = value_i;
  EXPECT_EQ(W, RMatrix(value_i, rows, cols));
}

TEST_P(RMatrixTest, SummationTest) {
  RMatrix AB = A + B;
  for (int n = 0; n < rows; ++n)
    for (int m = 0; m < cols; ++m) {
      EXPECT_DOUBLE_EQ(AB(n, m), A(n, m) + B(n, m));
    }

  if (rows == cols) {
    RMatrix W1 = A + Sym;
    RMatrix W2 = Sym + B;
    RMatrix W3 = A + Antisym;
    RMatrix W4 = Antisym + B;
    for (int n = 0; n < rows; ++n)
      for (int m = 0; m < cols; ++m) {
        EXPECT_DOUBLE_EQ(W1(n, m), A(n, m) + Sym(n, m));
        EXPECT_DOUBLE_EQ(W2(n, m), Sym(n, m) + B(n, m));
        EXPECT_DOUBLE_EQ(W3(n, m), A(n, m) + Antisym(n, m));
        EXPECT_DOUBLE_EQ(W4(n, m), Antisym(n, m) + B(n, m));
      }
  }
}

TEST_P(RMatrixTest, DifferenceTest) {
  RMatrix AB = A - B;
  for (int n = 0; n < rows; ++n)
    for (int m = 0; m < cols; ++m) {
      EXPECT_DOUBLE_EQ(AB(n, m), A(n, m) - B(n, m));
    }

  if (rows == cols) {
    RMatrix W1 = A - Sym;
    RMatrix W2 = Sym - B;
    RMatrix W3 = A - Antisym;
    RMatrix W4 = Antisym - B;
    for (int n = 0; n < rows; ++n)
      for (int m = 0; m < cols; ++m) {
        EXPECT_DOUBLE_EQ(W1(n, m), A(n, m) - Sym(n, m));
        EXPECT_DOUBLE_EQ(W2(n, m), Sym(n, m) - B(n, m));
        EXPECT_DOUBLE_EQ(W3(n, m), A(n, m) - Antisym(n, m));
        EXPECT_DOUBLE_EQ(W4(n, m), Antisym(n, m) - B(n, m));
      }
  }
}

TEST_P(RMatrixTest, AdditiveInversTest) {
  RMatrix W = -A;
  for (int n = 0; n < rows; ++n)
    for (int m = 0; m < cols; ++m)
      EXPECT_DOUBLE_EQ(W(n, m), -A(n, m));
}

TEST_P(RMatrixTest, ScalarMultiplicationTest) {
  RMatrix W1 = value * A;
  RMatrix W2 = A * value;
  RMatrix W3 = value_i * A;
  RMatrix W4 = A * value_i;
  for (int n = 0; n < rows; ++n)
    for (int m = 0; m < cols; ++m) {
      EXPECT_DOUBLE_EQ(W1(n, m), value * A(n, m));
      EXPECT_DOUBLE_EQ(W2(n, m), value * A(n, m));
      EXPECT_DOUBLE_EQ(W3(n, m), value_i * A(n, m));
      EXPECT_DOUBLE_EQ(W4(n, m), value_i * A(n, m));
    }
}

TEST_P(RMatrixTest, ScalarDivisionTest) {
  RMatrix W1 = A / value;
  RMatrix W2 = A / value_i;
  for (int n = 0; n < rows; ++n)
    for (int m = 0; m < cols; ++m) {
      EXPECT_DOUBLE_EQ(W1(n, m), A(n, m) / value);
      EXPECT_DOUBLE_EQ(W2(n, m), A(n, m) / value_i);
    }
}

TEST_P(RMatrixTest, MatrixVectorProductTest) {
  mpp_ba::SetTolerance(1e-13);
  RVector W0(rows);
  for (int n = 0; n < rows; ++n)
    for (int k = 0; k < cols; ++k) {
      W0[n] += A(n, k) * U[k];
    }
  EXPECT_EQ(W0, A * U);
}

TEST_P(RMatrixTest, MatrixMatrixProductTest) {
  RMatrix transB = transpose(B);
  RMatrix W0(rows, rows);
  for (int n = 0; n < rows; ++n)
    for (int m = 0; m < rows; ++m)
      for (int k = 0; k < cols; ++k) {
        W0[n][m] += A[n][k] * transB[k][m];
      }
  EXPECT_EQ(W0, A * transpose(B));

  if (rows == cols) {
    RMatrix W0(rows);
    RMatrix W1(rows);
    RMatrix V0(rows);
    RMatrix V1(rows);
    for (int n = 0; n < rows; ++n)
      for (int m = 0; m < rows; ++m)
        for (int k = 0; k < cols; ++k) {
          W0[n][m] += A[n][k] * Sym(k, m);
          W1[n][m] += Sym(n, k) * B[k][m];
          V0[n][m] += A[n][k] * Antisym(k, m);
          V1[n][m] += Antisym(n, k) * B[k][m];
        }
    EXPECT_EQ(W0, A * Sym);
    EXPECT_EQ(W1, Sym * B);
    EXPECT_EQ(V0, A * Antisym);
    EXPECT_EQ(V1, Antisym * B);
  }
}

TEST_P(RMatrixTest, TransposeTest) {
  RMatrix B = transpose(A);
  for (int n = 0; n < cols; ++n)
    for (int m = 0; m < rows; ++m)
      EXPECT_EQ(B[n][m], A[m][n]);
}

TEST_P(RMatrixTest, DiagonalTest) {
  RVector W = A.diag();
  for (int n = 0; n < std::min(rows, cols); ++n)
    EXPECT_DOUBLE_EQ(W[n], A[n][n]);

  if (rows = cols) {
    RMatrix W0(A);
    W0.diag(U);
    for (int n = 0; n < U.size(); ++n)
      for (int m = 0; m < U.size(); ++m) {
        if (n == m) {
          EXPECT_DOUBLE_EQ(W0(n, m), U[n]);
        } else {
          EXPECT_DOUBLE_EQ(W0(n, m), 0.0);
        }
      }
  }
}

TEST_P(RMatrixTest, EqualTest) {
  EXPECT_EQ(A, A);
  EXPECT_EQ(B, B);
  RMatrix W1(Sym);
  EXPECT_EQ(W1, Sym);
  RMatrix W2(Antisym);
  EXPECT_EQ(W2, Antisym);
}

TEST_P(RMatrixTest, NonEqualTest) {
  EXPECT_NE(A, B);
  EXPECT_NE(B, A);
}

TEST_P(RMatrixTest, SaveLoadTest) {
  Saver saver("SaveLoadTest");
  saver << A;
  saver.close();
  Loader loader("SaveLoadTest");
  RMatrix Z;
  loader >> Z;
  loader.close();
  EXPECT_EQ(A, Z);
}

#include "SparseRMatrix.hpp"

TEST_P(RMatrixTest, SparseTest) {
  int rows = RandomInt(500, 1000);
  int cols = RandomInt(500, 1000);
  RMatrix R(rows, cols);
  for (int i = 0; i < rows; ++i)
    for (int j = 0; j < cols; ++j) {
      if (double(rand()) / RAND_MAX < 0.5) R[i][j] = RandomDouble();
    }
  SparseRMatrix S(R);
  RMatrix T = S.ToRMatrix();
  EXPECT_EQ(R, T);

  for (int i = 0; i < rows; ++i)
    for (int j = 0; j < cols; ++j)
      EXPECT_EQ(R[i][j], S(i, j));
}

INSTANTIATE_TEST_CASE_P(
    BasicAlgebraMatrixTest, RMatrixTest,
    Values(BasicAlgebraMatrixTestParameter{1, 1}, BasicAlgebraMatrixTestParameter{1, 2},
           BasicAlgebraMatrixTestParameter{2, 1}, BasicAlgebraMatrixTestParameter{1, 3},
           BasicAlgebraMatrixTestParameter{2, 2}, BasicAlgebraMatrixTestParameter{3, 1},
           BasicAlgebraMatrixTestParameter{1, 4}, BasicAlgebraMatrixTestParameter{2, 3},
           BasicAlgebraMatrixTestParameter{3, 2}, BasicAlgebraMatrixTestParameter{4, 1},
           BasicAlgebraMatrixTestParameter{1, 5}, BasicAlgebraMatrixTestParameter{2, 4},
           BasicAlgebraMatrixTestParameter{3, 3}, BasicAlgebraMatrixTestParameter{4, 2},
           BasicAlgebraMatrixTestParameter{5, 1}));

int main(int argc, char **argv) {
  MppTest mppTest = MppTestBuilder(argc, argv).WithPPM();
  return mppTest.RUN_ALL_MPP_TESTS();
}
