#include "gtest/gtest.h"

#include "TestBasicAlgebra.hpp"
#include "basicalgebra/CMatrix.hpp"

class CMatrixTest : public BasicAlgebraMatrixTest {
protected:
  CMatrix A, B;
  RMatrix A_d;
  HermCMatrix Herm;
  SymRMatrix Sym;
  AntisymRMatrix Antisym;
  CVector U;
  RVector U_d;
  std::complex<double> value;
  double value_d;
  int value_i;
public:
  CMatrixTest() : BasicAlgebraMatrixTest() {
    mpp_ba::SetTolerance(5e-14);
    A.resize(rows, cols);
    B.resize(rows, cols);
    A_d.resize(rows, cols);
    U.resize(cols);
    U_d.resize(cols);

    for (int n = 0; n < cols; ++n) {
      U[n] = RandomComplex();
      U_d[n] = RandomDouble();
    }
    for (int n = 0; n < rows; ++n) {
      for (int m = 0; m < cols; ++m) {
        A(RandomComplex(), n, m);
        B[n][m] = RandomComplex();
        A_d[n][m] = RandomDouble();
      }
    }

    if (rows == cols) {
      Herm.resize(rows);
      Sym.resize(rows);
      Antisym.resize(rows);
      for (int n = 0; n < rows; ++n) {
        Herm(RandomDouble(), n, n);
        Sym(RandomDouble(), n, n);
        for (int m = n + 1; m < cols; ++m) {
          Herm(RandomComplex(), n, m);
          Sym(RandomDouble(), n, m);
          Antisym(RandomDouble(), n, m);
        }
      }
    }

    value = RandomNonZeroComplex();
    value_d = RandomNonZeroDouble();
    value_i = RandomNonZeroInt();
  }
};

TEST_P(CMatrixTest, ConstructorTest) {
  CMatrix W0(rows, cols);
  CMatrix W1(value, rows, cols);
  CMatrix W2(value_d, rows, cols);
  CMatrix W3(value_i, rows, cols);
  CMatrix W4(A);
  CMatrix W5(A_d);
  CMatrix W6(B);
  CMatrix W7(Herm);
  CMatrix W8(Sym);
  CMatrix W9(Antisym);
  for (int n = 0; n < rows; ++n)
    for (int m = 0; m < cols; ++m) {
      EXPECT_COMPLEX_EQ(W0(n, m), std::complex<double>(0.0));
      EXPECT_COMPLEX_EQ(W1(n, m), value);
      EXPECT_COMPLEX_EQ(W2(n, m), std::complex<double>(value_d));
      EXPECT_COMPLEX_EQ(W3(n, m), std::complex<double>(value_i));
      EXPECT_COMPLEX_EQ(W4(n, m), A(n, m));
      EXPECT_COMPLEX_EQ(W5(n, m), std::complex<double>(A_d(n, m)));
      EXPECT_COMPLEX_EQ(W6[n][m], B[n][m]);
      if (rows == cols) {
        EXPECT_COMPLEX_EQ(W7(n, m), Herm(n, m));
        EXPECT_COMPLEX_EQ(W8(n, m), std::complex<double>(Sym(n, m)));
        EXPECT_COMPLEX_EQ(W9(n, m), std::complex<double>(Antisym(n, m)));
      }
    }

  if (rows == cols) {
    CMatrix W(rows);
    for (int n = 0; n < rows; ++n)
      for (int m = 0; m < cols; ++m)
        EXPECT_COMPLEX_EQ(W[n][m], std::complex<double>(0.0));
  }
}

TEST_P(CMatrixTest, AssignmentTest) {
  CMatrix Z0 = A;
  EXPECT_EQ(Z0, A);
  Z0 = A_d;
  EXPECT_EQ(Z0, A_d);
  if (rows == cols) {
    Z0 = Herm;
    EXPECT_EQ(Z0, Herm);
    Z0 = Sym;
    EXPECT_EQ(Z0, Sym);
    Z0 = Antisym;
    EXPECT_EQ(Z0, Antisym);
  }
  CMatrix W(rows, cols);
  W = value;
  EXPECT_EQ(W, CMatrix(value, rows, cols));
  W = value_d;
  EXPECT_EQ(W, CMatrix(value_d, rows, cols));
  W = value_i;
  EXPECT_EQ(W, CMatrix(value_i, rows, cols));
}

TEST_P(CMatrixTest, SummationTest) {
  CMatrix AB = A + B;
  CMatrix W1 = A + A_d;
  CMatrix W2 = A_d + B;
  for (int n = 0; n < rows; ++n)
    for (int m = 0; m < cols; ++m) {
      EXPECT_COMPLEX_EQ(AB(n, m), A(n, m) + B(n, m));
      EXPECT_COMPLEX_EQ(W1(n, m), A(n, m) + A_d(n, m));
      EXPECT_COMPLEX_EQ(W2(n, m), A_d(n, m) + B(n, m));
    }

  if (rows == cols) {
    CMatrix W3 = A + Herm;
    CMatrix W4 = Herm + B;
    CMatrix W5 = A + Sym;
    CMatrix W6 = Sym + B;
    CMatrix W7 = A + Antisym;
    CMatrix W8 = Antisym + B;
    for (int n = 0; n < rows; ++n)
      for (int m = 0; m < cols; ++m) {
        EXPECT_COMPLEX_EQ(W3(n, m), A(n, m) + Herm(n, m));
        EXPECT_COMPLEX_EQ(W4(n, m), Herm(n, m) + B(n, m));
        EXPECT_COMPLEX_EQ(W5(n, m), A(n, m) + Sym(n, m));
        EXPECT_COMPLEX_EQ(W6(n, m), Sym(n, m) + B(n, m));
        EXPECT_COMPLEX_EQ(W7(n, m), A(n, m) + Antisym(n, m));
        EXPECT_COMPLEX_EQ(W8(n, m), Antisym(n, m) + B(n, m));
      }
  }
}

TEST_P(CMatrixTest, DifferenceTest) {
  CMatrix AB = A - B;
  CMatrix W1 = A - A_d;
  CMatrix W2 = A_d - B;
  for (int n = 0; n < rows; ++n)
    for (int m = 0; m < cols; ++m) {
      EXPECT_COMPLEX_EQ(AB(n, m), A(n, m) - B(n, m));
      EXPECT_COMPLEX_EQ(W1(n, m), A(n, m) - A_d(n, m));
      EXPECT_COMPLEX_EQ(W2(n, m), A_d(n, m) - B(n, m));
    }

  if (rows == cols) {
    CMatrix W3 = A - Herm;
    CMatrix W4 = Herm - B;
    CMatrix W5 = A - Sym;
    CMatrix W6 = Sym - B;
    CMatrix W7 = A - Antisym;
    CMatrix W8 = Antisym - B;
    for (int n = 0; n < rows; ++n)
      for (int m = 0; m < cols; ++m) {
        EXPECT_COMPLEX_EQ(W3(n, m), A(n, m) - Herm(n, m));
        EXPECT_COMPLEX_EQ(W4(n, m), Herm(n, m) - B(n, m));
        EXPECT_COMPLEX_EQ(W5(n, m), A(n, m) - Sym(n, m));
        EXPECT_COMPLEX_EQ(W6(n, m), Sym(n, m) - B(n, m));
        EXPECT_COMPLEX_EQ(W7(n, m), A(n, m) - Antisym(n, m));
        EXPECT_COMPLEX_EQ(W8(n, m), Antisym(n, m) - B(n, m));
      }
  }
}

TEST_P(CMatrixTest, AdditiveInversTest) {
  CMatrix W = -A;
  for (int n = 0; n < rows; ++n)
    for (int m = 0; m < cols; ++m)
      EXPECT_COMPLEX_EQ(W(n, m), -A(n, m));
}

TEST_P(CMatrixTest, ScalarMultiplicationTest) {
  CMatrix W1 = value * A;
  CMatrix W2 = A * value;
  CMatrix W3 = value_d * A;
  CMatrix W4 = A * value_d;
  CMatrix W5 = value_i * A;
  CMatrix W6 = A * value_i;
  for (int n = 0; n < rows; ++n)
    for (int m = 0; m < cols; ++m) {
      EXPECT_COMPLEX_EQ(W1(n, m), value * A(n, m));
      EXPECT_COMPLEX_EQ(W2(n, m), value * A(n, m));
      EXPECT_COMPLEX_EQ(W3(n, m), value_d * A(n, m));
      EXPECT_COMPLEX_EQ(W4(n, m), value_d * A(n, m));
      EXPECT_COMPLEX_EQ(W5(n, m), value_i * A(n, m));
      EXPECT_COMPLEX_EQ(W6(n, m), value_i * A(n, m));
    }
}

TEST_P(CMatrixTest, ScalarDivisionTest) {
  CMatrix W1 = A / value;
  CMatrix W2 = A / value_d;
  CMatrix W3 = A / value_i;
  for (int n = 0; n < rows; ++n)
    for (int m = 0; m < cols; ++m) {
      EXPECT_COMPLEX_EQ(W1(n, m), A(n, m) / value);
      EXPECT_COMPLEX_EQ(W2(n, m), A(n, m) / value_d);
      EXPECT_COMPLEX_EQ(W3(n, m), A(n, m) / value_i);
    }
}

TEST_P(CMatrixTest, MatrixVectorProductTest) {
  CVector W0(rows);
  CVector W1(rows);
  for (int n = 0; n < rows; ++n)
    for (int k = 0; k < cols; ++k) {
      W0[n] += A(n, k) * U[k];
      W1[n] += A(n, k) * U_d[k];
    }
  EXPECT_EQ(W0, A * U);
  EXPECT_EQ(W1, A * U_d);
}

TEST_P(CMatrixTest, MatrixMatrixProductTest) {
  CMatrix transB = transpose(B);
  CMatrix W0(rows, rows);
  CMatrix W1(rows, rows);
  for (int n = 0; n < rows; ++n)
    for (int m = 0; m < rows; ++m)
      for (int k = 0; k < cols; ++k) {
        W0[n][m] += A[n][k] * transB[k][m];
        W1[n][m] += A_d[n][k] * transB[k][m];
      }
  EXPECT_EQ(W0, A * transpose(B));
  EXPECT_EQ(W1, A_d * transpose(B));

  if (rows == cols) {
    CMatrix W0(rows);
    CMatrix W1(rows);
    CMatrix V0(rows);
    CMatrix V1(rows);
    CMatrix X0(rows);
    CMatrix X1(rows);
    for (int n = 0; n < rows; ++n)
      for (int m = 0; m < rows; ++m)
        for (int k = 0; k < cols; ++k) {
          W0[n][m] += A[n][k] * Sym(k, m);
          W1[n][m] += Sym(n, k) * B[k][m];
          V0[n][m] += A[n][k] * Antisym(k, m);
          V1[n][m] += Antisym(n, k) * B[k][m];
          X0[n][m] += A[n][k] * Herm(k, m);
          X1[n][m] += Herm(n, k) * B[k][m];
        }
    EXPECT_EQ(W0, A * Sym);
    EXPECT_EQ(W1, Sym * B);
    EXPECT_EQ(V0, A * Antisym);
    EXPECT_EQ(V1, Antisym * B);
    EXPECT_EQ(X0, A * Herm);
    EXPECT_EQ(X1, Herm * B);
  }
}

TEST_P(CMatrixTest, ConjTest) {
  CMatrix B = conj(A);
  for (int n = 0; n < rows; ++n)
    for (int m = 0; m < cols; ++m)
      EXPECT_COMPLEX_EQ(B(n, m), conj(A(n, m)));
}

TEST_P(CMatrixTest, RealTest) {
  RMatrix B = real(A);
  for (int n = 0; n < rows; ++n)
    for (int m = 0; m < cols; ++m)
      EXPECT_DOUBLE_EQ(B(n, m), real(A(n, m)));
}

TEST_P(CMatrixTest, ImagTest) {
  RMatrix B = imag(A);
  for (int n = 0; n < rows; ++n)
    for (int m = 0; m < cols; ++m)
      EXPECT_DOUBLE_EQ(B(n, m), imag(A(n, m)));
}

TEST_P(CMatrixTest, TransposeTest) {
  CMatrix B = transpose(A);
  for (int n = 0; n < cols; ++n)
    for (int m = 0; m < rows; ++m)
      EXPECT_COMPLEX_EQ(B(n, m), A(m, n));
}

TEST_P(CMatrixTest, AdjointTest) {
  CMatrix B = adjoint(A);
  for (int n = 0; n < cols; ++n)
    for (int m = 0; m < rows; ++m)
      EXPECT_COMPLEX_EQ(B(n, m), conj(A(m, n)));
}

TEST_P(CMatrixTest, DiagonalTest) {
  CVector W = A.diag();
  for (int n = 0; n < std::min(rows, cols); ++n)
    EXPECT_COMPLEX_EQ(W[n], A[n][n]);

  if (rows = cols) {
    CMatrix W0(A);
    W0.diag(U);
    CMatrix W1(A);
    W1.diag(U_d);
    for (int n = 0; n < U.size(); ++n)
      for (int m = 0; m < U.size(); ++m) {
        if (n == m) {
          EXPECT_COMPLEX_EQ(W0(n, m), U[n]);
          EXPECT_COMPLEX_EQ(W1(n, m), std::complex<double>(U_d[n]));
        } else {
          EXPECT_COMPLEX_EQ(W0(n, m), std::complex<double>(0.0));
          EXPECT_COMPLEX_EQ(W1(n, m), std::complex<double>(0.0));
        }
      }
  }
}

TEST_P(CMatrixTest, EqualTest) {
  EXPECT_EQ(A, A);
  EXPECT_EQ(B, B);
  CMatrix W1(Herm);
  EXPECT_EQ(W1, Herm);
  CMatrix W2(Sym);
  EXPECT_EQ(W2, Sym);
  CMatrix W3(Antisym);
  EXPECT_EQ(W3, Antisym);
  CMatrix W4(A_d);
  EXPECT_EQ(W4, A_d);
}

TEST_P(CMatrixTest, NonEqualTest) {
  EXPECT_NE(A, B);
  EXPECT_NE(B, A);
}

TEST_P(CMatrixTest, SaveLoadTest) {
  Saver saver("SaveLoadTest");
  saver << A;
  saver.close();
  Loader loader("SaveLoadTest");
  CMatrix Z;
  loader >> Z;
  loader.close();
  EXPECT_EQ(A, Z);
}

INSTANTIATE_TEST_CASE_P(
    IABasicAlgebraMatrixTest, CMatrixTest,
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

// #include "gtest/gtest.h"
//
// #include "basicalgebra/TestBasicAlgebra.hpp"
// #include "basicalgebra/CMatrix.hpp"
//
//
// class BasicAlgebraCMatrixTest : public BasicAlgebraMatrixTest {
// protected:
//     vector<double> valuesU_re, valuesU_im;
//     vector<vector<double>> valuesA_re, valuesB_re, valuesC_re, valuesS1, valuesS2,
//         valuesA1, valuesA2, valuesR1, valuesR2, valuesH1_re, valuesH2_re;
//     vector<vector<double>> valuesA_im, valuesB_im, valuesC_im, valuesH1_im, valuesH2_im;
////    vector<vector<int>> valuesA_int_re, valuesA_int_im;
//    double value_re, value_im;
//    int value_int;
//    CVector U;
//    RVector U_re;
//    CMatrix A, B, C;
//    RMatrix A_re;
////    CMatrixT<std::complex<int>> A_int;
//    HermCMatrix H1, H2;
//    RMatrix R1, R2;
//    SymRMatrix S1, S2;
//    AntisymRMatrix A1, A2;
// public:
//    BasicAlgebraCMatrixTest() : BasicAlgebraMatrixTest() {
//        valuesU_re.resize(cols);
//        valuesU_im.resize(cols);
//        valuesA_re.resize(rows);
//        valuesA_im.resize(rows);
//        valuesB_re.resize(rows);
//        valuesB_im.resize(rows);
//        valuesC_re.resize(cols);
//        valuesC_im.resize(cols);
////        valuesA_int_re.resize(rows);
////        valuesA_int_im.resize(rows);
//        valuesS1.resize(rows);
//        valuesS2.resize(cols);
//        valuesA1.resize(rows);
//        valuesA2.resize(cols);
//        valuesR1.resize(rows);
//        valuesR2.resize(cols);
//        valuesH1_re.resize(rows);
//        valuesH1_im.resize(rows);
//        valuesH2_re.resize(cols);
//        valuesH2_im.resize(cols);
//
//        U.resize(cols);
//        U_re.resize(cols);
//        A.resize(rows, cols);
//        B.resize(rows, cols);
//        C.resize(cols, rows);
//        A_re.resize(rows, cols);
////        A_int.resize(rows, cols);
//        H1.resize(rows);
//        H2.resize(cols);
//        R1.resize(rows);
//        R2.resize(cols);
//        S1.resize(rows);
//        S2.resize(cols);
//        A1.resize(rows);
//        A2.resize(cols);
//
//        for (int n = 0; n < cols; ++n) {
//            valuesU_re[n] = double(rand()) / RAND_MAX;
//            valuesU_im[n] = double(rand()) / RAND_MAX;
//            U[n] = std::complex<double>(valuesU_re[n], valuesU_im[n]);
//            U_re[n] = valuesU_re[n];
//            valuesC_re[n].resize(rows);
//            valuesC_im[n].resize(rows);
//            for (int m = 0; m < rows; ++m) {
//                valuesC_re[n][m] = double(rand()) / RAND_MAX;
//                valuesC_im[n][m] = double(rand()) / RAND_MAX;
//                C[n][m] = std::complex<double>(valuesC_re[n][m], valuesC_im[n][m]);
//            }
//            valuesH2_re[n].resize(cols);
//            valuesH2_im[n].resize(cols);
//            valuesR2[n].resize(cols);
//            valuesS2[n].resize(cols);
//            valuesA2[n].resize(cols);
//        }
//        for (int n = 0; n < rows; ++n) {
//            valuesA_re[n].resize(cols);
//            valuesA_im[n].resize(cols);
//            valuesB_re[n].resize(cols);
//            valuesB_im[n].resize(cols);
////            valuesA_int_re[n].resize(cols);
////            valuesA_int_im[n].resize(cols);
//            for (int m = 0; m < cols; ++m) {
//                valuesA_re[n][m] = double(rand()) / RAND_MAX;
//                valuesA_im[n][m] = double(rand()) / RAND_MAX;
//                valuesB_re[n][m] = double(rand()) / RAND_MAX;
//                valuesB_im[n][m] = double(rand()) / RAND_MAX;
//                A[n][m] = std::complex<double>(valuesA_re[n][m], valuesA_im[n][m]);
//                B(std::complex<double>(valuesB_re[n][m], valuesB_im[n][m]), n, m);
//                A_re[n][m] = valuesA_re[n][m];
////                valuesA_int_re[n][m] = rand();
////                valuesA_int_im[n][m] = rand();
////                A_int[n][m] =
////                    std::complex<int>(valuesA_int_re[n][m], valuesA_int_im[n][m]);
//            }
//            valuesH1_re[n].resize(rows);
//            valuesH1_im[n].resize(rows);
//            valuesR1[n].resize(rows);
//            valuesS1[n].resize(rows);
//            valuesA1[n].resize(rows);
//        }
//        for (int n = 0; n < rows; ++n) {
//            for (int m = 0; m < n; ++m) {
//                valuesH1_re[n][m] = double(rand()) / RAND_MAX;
//                valuesH1_re[m][n] = valuesH1_re[n][m];
//                valuesH1_im[n][m] = double(rand()) / RAND_MAX;
//                valuesH1_im[m][n] = -valuesH1_im[n][m];
//                valuesR1[n][m] = double(rand()) / RAND_MAX;
//                valuesS1[n][m] = double(rand()) / RAND_MAX;
//                valuesS1[m][n] = valuesS1[n][m];
//                valuesA1[n][m] = double(rand()) / RAND_MAX;
//                valuesA1[m][n] = -valuesA1[n][m];
//                H1(std::complex<double>(valuesH1_re[n][m], valuesH1_im[n][m]), n, m);
//                R1[n][m] = valuesR1[n][m];
//                S1(valuesS1[n][m], n, m);
//                A1(valuesA1[n][m], n, m);
//            }
//            valuesH1_re[n][n] = double(rand()) / RAND_MAX;
//            valuesH1_im[n][n] = 0.0;
//            valuesR1[n][n] = double(rand()) / RAND_MAX;
//            valuesS1[n][n] = double(rand()) / RAND_MAX;
//            valuesA1[n][n] = 0.0;
//            H1(std::complex<double>(valuesH1_re[n][n], valuesH1_im[n][n]), n, n);
//            R1[n][n] = valuesR1[n][n];
//            S1(valuesS1[n][n], n, n);
//            for (int m = n + 1; m < rows; ++m) {
//                valuesR1[n][m] = double(rand()) / RAND_MAX;
//                R1[n][m] = valuesR1[n][m];
//            }
//        }
//        for (int n = 0; n < cols; ++n) {
//            for (int m = 0; m < n; ++m) {
//                valuesH2_re[n][m] = double(rand()) / RAND_MAX;
//                valuesH2_re[m][n] = valuesH2_re[n][m];
//                valuesH2_im[n][m] = double(rand()) / RAND_MAX;
//                valuesH2_im[m][n] = -valuesH2_im[n][m];
//                valuesR2[n][m] = double(rand()) / RAND_MAX;
//                valuesS2[n][m] = double(rand()) / RAND_MAX;
//                valuesS2[m][n] = valuesS2[n][m];
//                valuesA2[n][m] = double(rand()) / RAND_MAX;
//                valuesA2[m][n] = -valuesA2[n][m];
//                H2(std::complex<double>(valuesH2_re[n][m], valuesH2_im[n][m]), n, m);
//                R2[n][m] = valuesR2[n][m];
//                S2(valuesS2[n][m], n, m);
//                A2(valuesA2[n][m], n, m);
//            }
//            valuesH2_re[n][n] = double(rand()) / RAND_MAX;
//            valuesH2_im[n][n] = 0.0;
//            valuesR2[n][n] = double(rand()) / RAND_MAX;
//            valuesS2[n][n] = double(rand()) / RAND_MAX;
//            valuesA2[n][n] = 0.0;
//            H2(std::complex<double>(valuesH2_re[n][n], valuesH2_im[n][n]), n, n);
//            R2[n][n] = valuesR2[n][n];
//            S2(valuesS2[n][n], n, n);
//            for (int m = n + 1; m < cols; ++m) {
//                valuesR2[n][m] = double(rand()) / RAND_MAX;
//                R2[n][m] = valuesR2[n][m];
//            }
//        }
//        value_re = double(rand()) / RAND_MAX;
//        value_im = double(rand()) / RAND_MAX;
//        value_int = rand();
//    }
//
//    std::complex<double> expectedProductAU(int n) {
//        double sum_re = 0.0, sum_im = 0.0;
//        for (int k = 0; k < cols; ++k) {
//            sum_re += valuesA_re[n][k] * valuesU_re[k] - valuesA_im[n][k] * valuesU_im[k];
//            sum_im += valuesA_re[n][k] * valuesU_im[k] + valuesA_im[n][k] * valuesU_re[k];
//        }
//        return std::complex<double>(sum_re, sum_im);
//    }
//
//    std::complex<double> expectedProductAU_re(int n) {
//        double sum_re = 0.0, sum_im = 0.0;
//        for (int k = 0; k < cols; ++k) {
//            sum_re += valuesA_re[n][k] * valuesU_re[k];
//            sum_im += valuesA_im[n][k] * valuesU_re[k];
//        }
//        return std::complex<double>(sum_re, sum_im);
//    }
//
//    std::complex<double> expectedProductAC(int n, int m) {
//        double sum_re = 0.0, sum_im = 0.0;
//        for (int k = 0; k < cols; ++k) {
//            sum_re +=
//                valuesA_re[n][k] * valuesC_re[k][m] - valuesA_im[n][k] * valuesC_im[k][m];
//            sum_im +=
//                valuesA_re[n][k] * valuesC_im[k][m] + valuesA_im[n][k] * valuesC_re[k][m];
//        }
//        return std::complex<double>(sum_re, sum_im);
//    }
//
//    std::complex<double> expectedProductH1H1(int n, int m) {
//        double sum_re = 0.0, sum_im = 0.0;
//        for (int k = 0; k < rows; ++k) {
//            sum_re += valuesH1_re[n][k] * valuesH1_re[k][m]
//                - valuesH1_im[n][k] * valuesH1_im[k][m];
//            sum_im += valuesH1_re[n][k] * valuesH1_im[k][m]
//                + valuesH1_im[n][k] * valuesH1_re[k][m];
//        }
//        return std::complex<double>(sum_re, sum_im);
//    }
//
//    std::complex<double> expectedProductH1A(int n, int m) {
//        double sum_re = 0.0, sum_im = 0.0;
//        for (int k = 0; k < rows; ++k) {
//            sum_re += valuesH1_re[n][k] * valuesA_re[k][m]
//                - valuesH1_im[n][k] * valuesA_im[k][m];
//            sum_im += valuesH1_re[n][k] * valuesA_im[k][m]
//                + valuesH1_im[n][k] * valuesA_re[k][m];
//        }
//        return std::complex<double>(sum_re, sum_im);
//    }
//
//    std::complex<double> expectedProductAH2(int n, int m) {
//        double sum_re = 0.0, sum_im = 0.0;
//        for (int k = 0; k < cols; ++k) {
//            sum_re += valuesA_re[n][k] * valuesH2_re[k][m]
//                - valuesA_im[n][k] * valuesH2_im[k][m];
//            sum_im += valuesA_re[n][k] * valuesH2_im[k][m]
//                + valuesA_im[n][k] * valuesH2_re[k][m];
//        }
//        return std::complex<double>(sum_re, sum_im);
//    }
//
//    std::complex<double> expectedProductR1A(int n, int m) {
//        double sum_re = 0.0, sum_im = 0.0;
//        for (int k = 0; k < rows; ++k) {
//            sum_re += valuesR1[n][k] * valuesA_re[k][m];
//            sum_im += valuesR1[n][k] * valuesA_im[k][m];
//        }
//        return std::complex<double>(sum_re, sum_im);
//    }
//
//    std::complex<double> expectedProductAR2(int n, int m) {
//        double sum_re = 0.0, sum_im = 0.0;
//        for (int k = 0; k < cols; ++k) {
//            sum_re += valuesA_re[n][k] * valuesR2[k][m];
//            sum_im += valuesA_im[n][k] * valuesR2[k][m];
//        }
//        return std::complex<double>(sum_re, sum_im);
//    }
//
//    std::complex<double> expectedProductS1A(int n, int m) {
//        double sum_re = 0.0, sum_im = 0.0;
//        for (int k = 0; k < rows; ++k) {
//            sum_re += valuesS1[n][k] * valuesA_re[k][m];
//            sum_im += valuesS1[n][k] * valuesA_im[k][m];
//        }
//        return std::complex<double>(sum_re, sum_im);
//    }
//
//    std::complex<double> expectedProductAS2(int n, int m) {
//        double sum_re = 0.0, sum_im = 0.0;
//        for (int k = 0; k < cols; ++k) {
//            sum_re += valuesA_re[n][k] * valuesS2[k][m];
//            sum_im += valuesA_im[n][k] * valuesS2[k][m];
//        }
//        return std::complex<double>(sum_re, sum_im);
//    }
//
//    std::complex<double> expectedProductA1A(int n, int m) {
//        double sum_re = 0.0, sum_im = 0.0;
//        for (int k = 0; k < rows; ++k) {
//            sum_re += valuesA1[n][k] * valuesA_re[k][m];
//            sum_im += valuesA1[n][k] * valuesA_im[k][m];
//        }
//        return std::complex<double>(sum_re, sum_im);
//    }
//
//    std::complex<double> expectedProductAA2(int n, int m) {
//        double sum_re = 0.0, sum_im = 0.0;
//        for (int k = 0; k < cols; ++k) {
//            sum_re += valuesA_re[n][k] * valuesA2[k][m];
//            sum_im += valuesA_im[n][k] * valuesA2[k][m];
//        }
//        return std::complex<double>(sum_re, sum_im);
//    }
//};
//
// TEST_P(BasicAlgebraCMatrixTest, ConstructorTest) {
//    CMatrix C0(rows, cols);
//    CMatrix C1(std::complex<double>(value_re, value_im), rows, cols);
//    CMatrix C2(value_re, rows, cols);
//    CMatrix C3(value_int, rows, cols);
//    CMatrix C4(A);
//    CMatrix C5(A_re);
//    CMatrix C6(H1);
//    CMatrix C7(R1);
//    CMatrix C8(S1);
//    CMatrix C9(A1);
//    CMatrix D1(U);
//    CMatrix D2(U_re);
//    for (int n = 0; n < rows; ++n) {
//        for (int m = 0; m < cols; ++m) {
//            EXPECT_DOUBLE_EQ(std::real(A[n][m]), valuesA_re[n][m]);
//            EXPECT_DOUBLE_EQ(std::imag(A[n][m]), valuesA_im[n][m]);
//            EXPECT_DOUBLE_EQ(std::real(B(n, m)), valuesB_re[n][m]);
//            EXPECT_DOUBLE_EQ(std::imag(B(n, m)), valuesB_im[n][m]);
//            EXPECT_DOUBLE_EQ(std::real(C0(n, m)), 0.0);
//            EXPECT_DOUBLE_EQ(std::imag(C0(n, m)), 0.0);
//            EXPECT_DOUBLE_EQ(std::real(C1(n, m)), value_re);
//            EXPECT_DOUBLE_EQ(std::imag(C1(n, m)), value_im);
//            EXPECT_DOUBLE_EQ(std::real(C2(n, m)), value_re);
//            EXPECT_DOUBLE_EQ(std::imag(C2(n, m)), 0.0);
//            EXPECT_DOUBLE_EQ(std::real(C3(n, m)), double(value_int));
//            EXPECT_DOUBLE_EQ(std::imag(C3(n, m)), 0.0);
//            EXPECT_DOUBLE_EQ(std::real(C4[n][m]), valuesA_re[n][m]);
//            EXPECT_DOUBLE_EQ(std::imag(C4[n][m]), valuesA_im[n][m]);
//            EXPECT_DOUBLE_EQ(std::real(C5[n][m]), valuesA_re[n][m]);
//            EXPECT_DOUBLE_EQ(std::imag(C5[n][m]), 0.0);
//        }
//        for (int m = 0; m < rows; ++m) {
//            EXPECT_DOUBLE_EQ(std::real(C6[n][m]), valuesH1_re[n][m]);
//            EXPECT_DOUBLE_EQ(std::imag(C6[n][m]), valuesH1_im[n][m]);
//            EXPECT_DOUBLE_EQ(std::real(C7[n][m]), valuesR1[n][m]);
//            EXPECT_DOUBLE_EQ(std::imag(C7[n][m]), 0.0);
//            EXPECT_DOUBLE_EQ(std::real(C8[n][m]), valuesS1[n][m]);
//            EXPECT_DOUBLE_EQ(std::imag(C8[n][m]), 0.0);
//            EXPECT_DOUBLE_EQ(std::real(C9[n][m]), valuesA1[n][m]);
//            EXPECT_DOUBLE_EQ(std::imag(C9[n][m]), 0.0);
//        }
//    }
//    for (int n = 0; n < cols; ++n)
//        for (int m = 0; m < cols; ++m) {
//            if (n == m) {
//                EXPECT_DOUBLE_EQ(std::real(D1[n][m]), valuesU_re[n]);
//                EXPECT_DOUBLE_EQ(std::imag(D1[n][m]), valuesU_im[n]);
//                EXPECT_DOUBLE_EQ(std::real(D2[n][m]), valuesU_re[n]);
//                EXPECT_DOUBLE_EQ(std::imag(D2[n][m]), 0.0);
//            } else {
//                EXPECT_DOUBLE_EQ(std::real(D1[n][m]), 0.0);
//                EXPECT_DOUBLE_EQ(std::imag(D1[n][m]), 0.0);
//                EXPECT_DOUBLE_EQ(std::real(D2[n][m]), 0.0);
//                EXPECT_DOUBLE_EQ(std::imag(D2[n][m]), 0.0);
//            }
//        }
//}
//
// TEST_P(BasicAlgebraCMatrixTest, AssignmentTest) {
//    CMatrix C1(rows, cols);
//    C1 = std::complex<double>(value_re, value_im);
//    CMatrix C2(rows, cols);
//    C2 = value_re;
//    CMatrix C3(rows, cols);
//    C3 = value_int;
//    CMatrix C4;
//    C4 = A;
//    CMatrix C5;
//    C5 = A_re;
//    CMatrix C6;
//    C6 = H1;
//    CMatrix C7;
//    C7 = R1;
//    CMatrix C8;
//    C8 = S1;
//    CMatrix C9;
//    C9 = A1;
//    for (int n = 0; n < rows; ++n) {
//        for (int m = 0; m < cols; ++m) {
//            EXPECT_DOUBLE_EQ(std::real(C1(n, m)), value_re);
//            EXPECT_DOUBLE_EQ(std::imag(C1(n, m)), value_im);
//            EXPECT_DOUBLE_EQ(std::real(C2(n, m)), value_re);
//            EXPECT_DOUBLE_EQ(std::imag(C2(n, m)), 0.0);
//            EXPECT_DOUBLE_EQ(std::real(C3(n, m)), double(value_int));
//            EXPECT_DOUBLE_EQ(std::imag(C3(n, m)), 0.0);
//            EXPECT_DOUBLE_EQ(std::real(C4[n][m]), valuesA_re[n][m]);
//            EXPECT_DOUBLE_EQ(std::imag(C4[n][m]), valuesA_im[n][m]);
//            EXPECT_DOUBLE_EQ(std::real(C5[n][m]), valuesA_re[n][m]);
//            EXPECT_DOUBLE_EQ(std::imag(C5[n][m]), 0.0);
//        }
//        for (int m = 0; m < rows; ++m) {
//            EXPECT_DOUBLE_EQ(std::real(C6[n][m]), valuesH1_re[n][m]);
//            EXPECT_DOUBLE_EQ(std::imag(C6[n][m]), valuesH1_im[n][m]);
//            EXPECT_DOUBLE_EQ(std::real(C7[n][m]), valuesR1[n][m]);
//            EXPECT_DOUBLE_EQ(std::imag(C7[n][m]), 0.0);
//            EXPECT_DOUBLE_EQ(std::real(C8[n][m]), valuesS1[n][m]);
//            EXPECT_DOUBLE_EQ(std::imag(C8[n][m]), 0.0);
//            EXPECT_DOUBLE_EQ(std::real(C9[n][m]), valuesA1[n][m]);
//            EXPECT_DOUBLE_EQ(std::imag(C9[n][m]), 0.0);
//        }
//    }
//}
//
// TEST_P(BasicAlgebraCMatrixTest, SummationTest) {
//    CMatrix AB = A + B;
//    CMatrix BA_re = B + A_re;
//    CMatrix A_reB = A_re + B;
//    for (int n = 0; n < rows; ++n)
//        for (int m = 0; m < cols; ++m) {
//            EXPECT_DOUBLE_EQ(std::real(AB[n][m]),
//                             valuesA_re[n][m] + valuesB_re[n][m]);
//            EXPECT_DOUBLE_EQ(std::imag(AB[n][m]),
//                             valuesA_im[n][m] + valuesB_im[n][m]);
//            EXPECT_DOUBLE_EQ(std::real(BA_re[n][m]),
//                             valuesB_re[n][m] + valuesA_re[n][m]);
//            EXPECT_DOUBLE_EQ(std::imag(BA_re[n][m]), valuesB_im[n][m]);
//            EXPECT_DOUBLE_EQ(std::real(A_reB[n][m]),
//                             valuesA_re[n][m] + valuesB_re[n][m]);
//            EXPECT_DOUBLE_EQ(std::imag(A_reB[n][m]), valuesB_im[n][m]);
//        }
//    if (rows == cols) {
//        CMatrix AH1 = A + H1;
//        CMatrix H1A = H1 + A;
//        CMatrix AS1 = A + S1;
//        CMatrix S1A = S1 + A;
//        CMatrix AA1 = A + A1;
//        CMatrix A1A = A1 + A;
//        for (int n = 0; n < rows; ++n)
//            for (int m = 0; m < cols; ++m) {
//                EXPECT_DOUBLE_EQ(std::real(AH1[n][m]),
//                                 valuesA_re[n][m] + valuesH1_re[n][m]);
//                EXPECT_DOUBLE_EQ(std::imag(AH1[n][m]),
//                                 valuesA_im[n][m] + valuesH1_im[n][m]);
//                EXPECT_DOUBLE_EQ(std::real(H1A[n][m]),
//                                 valuesH1_re[n][m] + valuesA_re[n][m]);
//                EXPECT_DOUBLE_EQ(std::imag(H1A[n][m]),
//                                 valuesH1_im[n][m] + valuesA_im[n][m]);
//                EXPECT_DOUBLE_EQ(std::real(AS1[n][m]),
//                                 valuesA_re[n][m] + valuesS1[n][m]);
//                EXPECT_DOUBLE_EQ(std::imag(AS1[n][m]), valuesA_im[n][m]);
//                EXPECT_DOUBLE_EQ(std::real(S1A[n][m]),
//                                 valuesS1[n][m] + valuesA_re[n][m]);
//                EXPECT_DOUBLE_EQ(std::imag(S1A[n][m]), valuesA_im[n][m]);
//                EXPECT_DOUBLE_EQ(std::real(AA1[n][m]),
//                                 valuesA_re[n][m] + valuesA1[n][m]);
//                EXPECT_DOUBLE_EQ(std::imag(AA1[n][m]), valuesA_im[n][m]);
//                EXPECT_DOUBLE_EQ(std::real(A1A[n][m]),
//                                 valuesA1[n][m] + valuesA_re[n][m]);
//                EXPECT_DOUBLE_EQ(std::imag(A1A[n][m]), valuesA_im[n][m]);
//            }
//    }
//}
//
// TEST_P(BasicAlgebraCMatrixTest, DifferenceTest) {
//    CMatrix AB = A - B;
//    CMatrix BA_re = B - A_re;
//    CMatrix A_reB = A_re - B;
//    for (int n = 0; n < rows; ++n)
//        for (int m = 0; m < cols; ++m) {
//            EXPECT_DOUBLE_EQ(std::real(AB[n][m]),
//                             valuesA_re[n][m] - valuesB_re[n][m]);
//            EXPECT_DOUBLE_EQ(std::imag(AB[n][m]),
//                             valuesA_im[n][m] - valuesB_im[n][m]);
//            EXPECT_DOUBLE_EQ(std::real(BA_re[n][m]),
//                             valuesB_re[n][m] - valuesA_re[n][m]);
//            EXPECT_DOUBLE_EQ(std::imag(BA_re[n][m]), valuesB_im[n][m]);
//            EXPECT_DOUBLE_EQ(std::real(A_reB[n][m]),
//                             valuesA_re[n][m] - valuesB_re[n][m]);
//            EXPECT_DOUBLE_EQ(std::imag(A_reB[n][m]), -valuesB_im[n][m]);
//        }
//    if (rows == cols) {
//        CMatrix AH1 = A - H1;
//        CMatrix H1A = H1 - A;
//        CMatrix AS1 = A - S1;
//        CMatrix S1A = S1 - A;
//        CMatrix AA1 = A - A1;
//        CMatrix A1A = A1 - A;
//        for (int n = 0; n < rows; ++n)
//            for (int m = 0; m < cols; ++m) {
//                EXPECT_DOUBLE_EQ(std::real(AH1[n][m]),
//                                 valuesA_re[n][m] - valuesH1_re[n][m]);
//                EXPECT_DOUBLE_EQ(std::imag(AH1[n][m]),
//                                 valuesA_im[n][m] - valuesH1_im[n][m]);
//                EXPECT_DOUBLE_EQ(std::real(H1A[n][m]),
//                                 valuesH1_re[n][m] - valuesA_re[n][m]);
//                EXPECT_DOUBLE_EQ(std::imag(H1A[n][m]),
//                                 valuesH1_im[n][m] - valuesA_im[n][m]);
//                EXPECT_DOUBLE_EQ(std::real(AS1[n][m]),
//                                 valuesA_re[n][m] - valuesS1[n][m]);
//                EXPECT_DOUBLE_EQ(std::imag(AS1[n][m]), valuesA_im[n][m]);
//                EXPECT_DOUBLE_EQ(std::real(S1A[n][m]),
//                                 valuesS1[n][m] - valuesA_re[n][m]);
//                EXPECT_DOUBLE_EQ(std::imag(S1A[n][m]), -valuesA_im[n][m]);
//                EXPECT_DOUBLE_EQ(std::real(AA1[n][m]),
//                                 valuesA_re[n][m] - valuesA1[n][m]);
//                EXPECT_DOUBLE_EQ(std::imag(AA1[n][m]), valuesA_im[n][m]);
//                EXPECT_DOUBLE_EQ(std::real(A1A[n][m]),
//                                 valuesA1[n][m] - valuesA_re[n][m]);
//                EXPECT_DOUBLE_EQ(std::imag(A1A[n][m]), -valuesA_im[n][m]);
//            }
//    }
//}
//
// TEST_P(BasicAlgebraCMatrixTest, AdditiveInversTest) {
//    CMatrix B = -A;
//    for (int n = 0; n < rows; ++n)
//        for (int m = 0; m < cols; ++m) {
//            EXPECT_DOUBLE_EQ(std::real(B[n][m]), -valuesA_re[n][m]);
//            EXPECT_DOUBLE_EQ(std::imag(B[n][m]), -valuesA_im[n][m]);
//        }
//}
//
// TEST_P(BasicAlgebraCMatrixTest, ScalarMultiplicationTest) {
//    CMatrix aA = std::complex<double>(value_re, value_im) * A;
//    CMatrix a_reA = value_re * A;
//    CMatrix a_intA = value_int * A;
//    CMatrix Aa = A * std::complex<double>(value_re, value_im);
//    CMatrix Aa_re = A * value_re;
//    CMatrix Aa_int = A * value_int;
//    for (int n = 0; n < rows; ++n)
//        for (int m = 0; m < cols; ++m) {
//            EXPECT_DOUBLE_EQ(std::real(aA[n][m]),
//                             value_re *valuesA_re[n][m]
//                                 -value_im * valuesA_im[n][m]);
//            EXPECT_DOUBLE_EQ(std::imag(aA[n][m]),
//                             value_re *valuesA_im[n][m]
//                                 +value_im * valuesA_re[n][m]);
//            EXPECT_DOUBLE_EQ(std::real(a_reA[n][m]), value_re * valuesA_re[n][m]);
//            EXPECT_DOUBLE_EQ(std::imag(a_reA[n][m]), value_re * valuesA_im[n][m]);
//            EXPECT_DOUBLE_EQ(std::real(a_intA[n][m]), value_int * valuesA_re[n][m]);
//            EXPECT_DOUBLE_EQ(std::imag(a_intA[n][m]), value_int * valuesA_im[n][m]);
//            EXPECT_DOUBLE_EQ(std::real(Aa[n][m]),
//                             value_re *valuesA_re[n][m]
//                                 -value_im * valuesA_im[n][m]);
//            EXPECT_DOUBLE_EQ(std::imag(Aa[n][m]),
//                             value_re *valuesA_im[n][m]
//                                 +value_im * valuesA_re[n][m]);
//            EXPECT_DOUBLE_EQ(std::real(Aa_re[n][m]), value_re * valuesA_re[n][m]);
//            EXPECT_DOUBLE_EQ(std::imag(Aa_re[n][m]), value_re * valuesA_im[n][m]);
//            EXPECT_DOUBLE_EQ(std::real(Aa_int[n][m]), value_int * valuesA_re[n][m]);
//            EXPECT_DOUBLE_EQ(std::imag(Aa_int[n][m]), value_int * valuesA_im[n][m]);
//        }
//}
//
// TEST_P(BasicAlgebraCMatrixTest, ScalarDivisionTest) {
//    CMatrix Aa = A / std::complex<double>(value_re, value_im);
//    CMatrix Aa_re = A / value_re;
//    CMatrix Aa_int = A / value_int;
//    for (int n = 0; n < rows; ++n)
//        for (int m = 0; m < cols; ++m) {
//            EXPECT_DOUBLE_EQ(std::real(Aa[n][m]),
//                             (value_re * valuesA_re[n][m]
//                                 + value_im * valuesA_im[n][m])
//                                 / (value_re * value_re + value_im * value_im));
//            EXPECT_DOUBLE_EQ(std::imag(Aa[n][m]),
//                             (value_re * valuesA_im[n][m]
//                                 - value_im * valuesA_re[n][m])
//                                 / (value_re * value_re + value_im * value_im));
//            EXPECT_DOUBLE_EQ(std::real(Aa_re[n][m]), valuesA_re[n][m] / value_re);
//            EXPECT_DOUBLE_EQ(std::imag(Aa_re[n][m]), valuesA_im[n][m] / value_re);
//            EXPECT_DOUBLE_EQ(std::real(Aa_int[n][m]), valuesA_re[n][m] / value_int);
//            EXPECT_DOUBLE_EQ(std::imag(Aa_int[n][m]), valuesA_im[n][m] / value_int);
//        }
//}
//
// TEST_P(BasicAlgebraCMatrixTest, MatrixVectorProductTest) {
//    CVector AU = A * U;
//    CVector AU_re = A * U_re;
//    for (int n = 0; n < rows; ++n) {
//        EXPECT_NEAR(std::real(AU[n]), std::real(expectedProductAU(n)), tol);
//        EXPECT_NEAR(std::imag(AU[n]), std::imag(expectedProductAU(n)), tol);
//        EXPECT_NEAR(std::real(AU_re[n]), std::real(expectedProductAU_re(n)), tol);
//        EXPECT_NEAR(std::imag(AU_re[n]), std::imag(expectedProductAU_re(n)), tol);
//    }
//}
//
// TEST_P(BasicAlgebraCMatrixTest, MatrixMatrixProductTest) {
//    CMatrix AC = A * C;
//    CMatrix H1H1 = H1 * H1;
//    CMatrix AH2 = A * H2;
//    CMatrix H1A = H1 * A;
//    CMatrix AR2 = A * R2;
//    CMatrix R1A = R1 * A;
//    CMatrix AS2 = A * S2;
//    CMatrix S1A = S1 * A;
//    CMatrix AA2 = A * A2;
//    CMatrix A1A = A1 * A;
//
//    for (int n = 0; n < AC.rows(); ++n)
//        for (int m = 0; m < AC.cols(); ++m) {
//            EXPECT_NEAR(std::real(AC[n][m]), std::real(expectedProductAC(n, m)), tol);
//            EXPECT_NEAR(std::imag(AC[n][m]), std::imag(expectedProductAC(n, m)), tol);
//        }
//    for (int n = 0; n < H1H1.rows(); ++n)
//        for (int m = 0; m < H1H1.cols(); ++m) {
//            EXPECT_NEAR(std::real(H1H1[n][m]),
//                        std::real(expectedProductH1H1(n, m)),
//                        tol);
//            EXPECT_NEAR(std::imag(H1H1[n][m]),
//                        std::imag(expectedProductH1H1(n, m)),
//                        tol);
//        }
//
//    for (int n = 0; n < rows; ++n)
//        for (int m = 0; m < cols; ++m) {
//            EXPECT_NEAR(std::real(H1A[n][m]),
//                        std::real(expectedProductH1A(n, m)),
//                        tol);
//            EXPECT_NEAR(std::imag(H1A[n][m]),
//                        std::imag(expectedProductH1A(n, m)),
//                        tol);
//            EXPECT_NEAR(std::real(AH2[n][m]),
//                        std::real(expectedProductAH2(n, m)),
//                        tol);
//            EXPECT_NEAR(std::imag(AH2[n][m]),
//                        std::imag(expectedProductAH2(n, m)),
//                        tol);
//            EXPECT_NEAR(std::real(R1A[n][m]),
//                        std::real(expectedProductR1A(n, m)),
//                        tol);
//            EXPECT_NEAR(std::imag(R1A[n][m]),
//                        std::imag(expectedProductR1A(n, m)),
//                        tol);
//            EXPECT_NEAR(std::real(AR2[n][m]),
//                        std::real(expectedProductAR2(n, m)),
//                        tol);
//            EXPECT_NEAR(std::imag(AR2[n][m]),
//                        std::imag(expectedProductAR2(n, m)),
//                        tol);
//            EXPECT_NEAR(std::real(S1A[n][m]),
//                        std::real(expectedProductS1A(n, m)),
//                        tol);
//            EXPECT_NEAR(std::imag(S1A[n][m]),
//                        std::imag(expectedProductS1A(n, m)),
//                        tol);
//            EXPECT_NEAR(std::real(AS2[n][m]),
//                        std::real(expectedProductAS2(n, m)),
//                        tol);
//            EXPECT_NEAR(std::imag(AS2[n][m]),
//                        std::imag(expectedProductAS2(n, m)),
//                        tol);
//            EXPECT_NEAR(std::real(A1A[n][m]),
//                        std::real(expectedProductA1A(n, m)),
//                        tol);
//            EXPECT_NEAR(std::imag(A1A[n][m]),
//                        std::imag(expectedProductA1A(n, m)),
//                        tol);
//            EXPECT_NEAR(std::real(AA2[n][m]),
//                        std::real(expectedProductAA2(n, m)),
//                        tol);
//            EXPECT_NEAR(std::imag(AA2[n][m]),
//                        std::imag(expectedProductAA2(n, m)),
//                        tol);
//        }
//}
//
// TEST_P(BasicAlgebraCMatrixTest, TransposeTest) {
//    CMatrix B = transpose(A);
//    for (int n = 0; n < B.rows(); ++n)
//        for (int m = 0; m < B.cols(); ++m) {
//            EXPECT_DOUBLE_EQ(std::real(B[n][m]), valuesA_re[m][n]);
//            EXPECT_DOUBLE_EQ(std::imag(B[n][m]), valuesA_im[m][n]);
//        }
//}
//
// TEST_P(BasicAlgebraCMatrixTest, ConjTest) {
//    CMatrix B = conj(A);
//    for (int n = 0; n < B.rows(); ++n)
//        for (int m = 0; m < B.cols(); ++m) {
//            EXPECT_DOUBLE_EQ(std::real(B[n][m]), valuesA_re[n][m]);
//            EXPECT_DOUBLE_EQ(std::imag(B[n][m]), -valuesA_im[n][m]);
//        }
//}
//
// TEST_P(BasicAlgebraCMatrixTest, AdjointTest) {
//    CMatrix B = adjoint(A);
//    for (int n = 0; n < B.rows(); ++n)
//        for (int m = 0; m < B.cols(); ++m) {
//            EXPECT_DOUBLE_EQ(std::real(B[n][m]), valuesA_re[m][n]);
//            EXPECT_DOUBLE_EQ(std::imag(B[n][m]), -valuesA_im[m][n]);
//        }
//}
//
// TEST_P(BasicAlgebraCMatrixTest, RealPartTest) {
//    RMatrix B = real(A);
//    for (int n = 0; n < B.rows(); ++n)
//        for (int m = 0; m < B.cols(); ++m)
//            EXPECT_DOUBLE_EQ(B[n][m], valuesA_re[n][m]);
//}
//
// TEST_P(BasicAlgebraCMatrixTest, ImaginaryPartTest) {
//    RMatrix B = imag(A);
//    for (int n = 0; n < B.rows(); ++n)
//        for (int m = 0; m < B.cols(); ++m)
//            EXPECT_DOUBLE_EQ(B[n][m], valuesA_im[n][m]);
//}
//
// TEST_P(BasicAlgebraCMatrixTest, DiagonalTest) {
//    CMatrix D(A);
//    D.diag(U);
//    CVector d = A.diag();
//
//    for (int n = 0; n < cols; ++n)
//        for (int m = 0; m < cols; ++m)
//            if (n == m) {
//                EXPECT_DOUBLE_EQ(std::real(D[n][m]), valuesU_re[n]);
//                EXPECT_DOUBLE_EQ(std::imag(D[n][m]), valuesU_im[n]);
//            } else {
//                EXPECT_DOUBLE_EQ(std::real(D[n][m]), 0.0);
//                EXPECT_DOUBLE_EQ(std::imag(D[n][m]), 0.0);
//            }
//    for (int n = 0; n < std::min(rows, cols); ++n) {
//        EXPECT_DOUBLE_EQ(std::real(d[n]), valuesA_re[n][n]);
//        EXPECT_DOUBLE_EQ(std::imag(d[n]), valuesA_im[n][n]);
//    }
//}
//
// TEST_P(BasicAlgebraCMatrixTest, EqualTest) {
//    EXPECT_EQ(A, A);
//    EXPECT_EQ(B, B);
//    CMatrix h(H1);
//    EXPECT_EQ(h, H1);
//    EXPECT_EQ(H1, h);
//    CMatrix r(R1);
//    EXPECT_EQ(r, R1);
//    EXPECT_EQ(R1, r);
//    CMatrix s(S1);
//    EXPECT_EQ(s, S1);
//    EXPECT_EQ(S1, s);
//    CMatrix a(A1);
//    EXPECT_EQ(a, A1);
//    EXPECT_EQ(A1, a);
//}
//
// TEST_P(BasicAlgebraCMatrixTest, NonEqualTest) {
//    EXPECT_NE(A, B);
//    EXPECT_NE(B, A);
//    EXPECT_NE(A, H1);
//    EXPECT_NE(H1, A);
//    EXPECT_NE(A, R1);
//    EXPECT_NE(R1, A);
//    EXPECT_NE(A, S1);
//    EXPECT_NE(S1, A);
//    EXPECT_NE(A, A1);
//    EXPECT_NE(A1, A);
//}
//
// TEST_P(BasicAlgebraCMatrixTest, SaveLoadTest) {
//    Saver saver("SaveLoadTest");
//    saver << A;
//    saver.close();
//    Loader loader("SaveLoadTest");
//    CMatrix C;
//    loader >> C;
//    loader.close();
//    for (int n = 0; n < rows; ++n)
//        for (int m = 0; m < cols; ++m) {
//            EXPECT_DOUBLE_EQ(std::real(C[n][m]), valuesA_re[n][m]);
//            EXPECT_DOUBLE_EQ(std::imag(C[n][m]), valuesA_im[n][m]);
//        }
//}
//
// INSTANTIATE_TEST_CASE_P(BasicAlgebraMatrixTest, BasicAlgebraCMatrixTest, Values(
//    BasicAlgebraMatrixTestParameter{1, 1},
//    BasicAlgebraMatrixTestParameter{1, 2},
//    BasicAlgebraMatrixTestParameter{1, 3},
//    BasicAlgebraMatrixTestParameter{2, 1},
//    BasicAlgebraMatrixTestParameter{2, 2},
//    BasicAlgebraMatrixTestParameter{2, 3},
//    BasicAlgebraMatrixTestParameter{3, 1},
//    BasicAlgebraMatrixTestParameter{3, 2},
//    BasicAlgebraMatrixTestParameter{3, 3}
//));
//
// int main(int argc, char **argv) {
//    MppTest mppTest = MppTestBuilder(argc, argv).WithPPM();
//    return mppTest.RUN_ALL_MPP_TESTS();
//}
