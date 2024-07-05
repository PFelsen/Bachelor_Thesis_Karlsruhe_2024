#include "gtest/gtest.h"

#include "TestBasicAlgebra.hpp"
#include "basicalgebra/CVector.hpp"

class CVectorTest : public BasicAlgebraVectorTest {
protected:
  CVector U, V;
  RVector U_r;
  std::complex<double> value;
  double value_d;
  int value_i;
public:
  CVectorTest() : BasicAlgebraVectorTest() {
    U.resize(dim);
    V.resize(dim);
    U_r.resize(dim);

    for (int n = 0; n < dim; ++n) {
      U[n] = RandomComplex();
      V[n] = RandomComplex();
      U_r[n] = RandomDouble();
    }
    value = RandomNonZeroComplex();
    value_d = RandomNonZeroDouble();
    value_i = RandomNonZeroInt();
  }
};

TEST_P(CVectorTest, ConstructorTest) {
  CVector W0(dim);
  CVector W1(value, dim);
  CVector W2(value_d, dim);
  CVector W3(value_i, dim);
  CVector W4(U);
  CVector W5(U_r);

  for (int n = 0; n < dim; ++n) {
    EXPECT_COMPLEX_EQ(W0[n], std::complex<double>(0.0));
    EXPECT_COMPLEX_EQ(W1[n], value);
    EXPECT_COMPLEX_EQ(W2[n], std::complex<double>(value_d));
    EXPECT_COMPLEX_EQ(W3[n], std::complex<double>(value_i));
    EXPECT_COMPLEX_EQ(W4[n], U[n]);
    EXPECT_COMPLEX_EQ(W5[n], std::complex<double>(U_r[n]));
  }
}

TEST_P(CVectorTest, AssignmentTest) {
  CVector Z0 = U;
  EXPECT_EQ(Z0, U);
  CVector Z1;
  Z1 = U_r;
  EXPECT_EQ(Z1, CVector(U_r));
  CVector W(dim);
  W = value;
  EXPECT_EQ(W, CVector(value, dim));
  W = value_d;
  EXPECT_EQ(W, CVector(value_d, dim));
  W = value_i;
  EXPECT_EQ(W, CVector(value_i, dim));
}

TEST_P(CVectorTest, SummationTest) {
  CVector UV = U + V;
  CVector W1 = U + U_r;
  CVector W2 = U_r + V;
  CVector W3(U);
  W3 += value;
  CVector W4(U);
  W4 += value_d;
  CVector W5(U);
  W5 += value_i;
  for (int n = 0; n < dim; ++n) {
    EXPECT_COMPLEX_EQ(UV[n], U[n] + V[n]);
    EXPECT_COMPLEX_EQ(W1[n], U[n] + U_r[n]);
    EXPECT_COMPLEX_EQ(W2[n], U_r[n] + V[n]);
    EXPECT_COMPLEX_EQ(W3[n], U[n] + value);
    EXPECT_COMPLEX_EQ(W4[n], U[n] + std::complex<double>(value_d));
    EXPECT_COMPLEX_EQ(W5[n], U[n] + std::complex<double>(value_i));
  }
}

TEST_P(CVectorTest, DifferenceTest) {
  CVector UV = U - V;
  CVector W1 = U - U_r;
  CVector W2 = U_r - V;
  CVector W3(U);
  W3 -= value;
  CVector W4(U);
  W4 -= value_d;
  CVector W5(U);
  W5 -= value_i;
  for (int n = 0; n < dim; ++n) {
    EXPECT_COMPLEX_EQ(UV[n], U[n] - V[n]);
    EXPECT_COMPLEX_EQ(W1[n], U[n] - U_r[n]);
    EXPECT_COMPLEX_EQ(W2[n], U_r[n] - V[n]);
    EXPECT_COMPLEX_EQ(W3[n], U[n] - value);
    EXPECT_COMPLEX_EQ(W4[n], U[n] - std::complex<double>(value_d));
    EXPECT_COMPLEX_EQ(W5[n], U[n] - std::complex<double>(value_i));
  }
}

TEST_P(CVectorTest, AdditiveInversTest) {
  CVector V = -U;
  for (int n = 0; n < dim; ++n)
    EXPECT_EQ(V[n], -U[n]);
}

TEST_P(CVectorTest, ScalarMultiplicationTest) {
  CVector W1 = value * U;
  CVector W2 = U * value;
  CVector W3 = value_d * U;
  CVector W4 = U * value_d;
  CVector W5 = value_i * U;
  CVector W6 = U * value_i;
  for (int n = 0; n < dim; ++n) {
    EXPECT_COMPLEX_EQ(W1[n], value * U[n]);
    EXPECT_COMPLEX_EQ(W2[n], value * U[n]);
    EXPECT_COMPLEX_EQ(W3[n], value_d * U[n]);
    EXPECT_COMPLEX_EQ(W4[n], value_d * U[n]);
    EXPECT_COMPLEX_EQ(W5[n], value_i * U[n]);
    EXPECT_COMPLEX_EQ(W6[n], value_i * U[n]);
  }
}

TEST_P(CVectorTest, ScalarDivisionTest) {
  CVector W1 = U / value;
  CVector W2 = U / value_d;
  CVector W3 = U / value_i;
  for (int n = 0; n < dim; ++n) {
    EXPECT_COMPLEX_EQ(W1[n], U[n] / value);
    EXPECT_COMPLEX_EQ(W2[n], U[n] / value_d);
    EXPECT_COMPLEX_EQ(W3[n], U[n] / value_i);
  }
}

TEST_P(CVectorTest, ScalarProductTest) {
  std::complex<double> a0{}, a1{}, a2{};
  for (int n = 0; n < dim; ++n) {
    a0 += U[n] * conj(V[n]);
    a1 += U[n] * U_r[n];
    a2 += U_r[n] * conj(V[n]);
  }
  EXPECT_COMPLEX_EQ(a0, U * V);
  EXPECT_COMPLEX_EQ(a1, U * U_r);
  EXPECT_COMPLEX_EQ(a2, U_r * V);
}

TEST_P(CVectorTest, NormSqrTest) {
  double sum = 0.0;
  for (int n = 0; n < dim; ++n)
    sum += std::norm(U[n]);
  EXPECT_DOUBLE_EQ(normSqr(U), sum);
}

TEST_P(CVectorTest, NormTest) {
  double sum = 0.0;
  for (int n = 0; n < dim; ++n)
    sum += std::norm(U[n]);
  EXPECT_DOUBLE_EQ(norm(U), sqrt(sum));
}

TEST_P(CVectorTest, EqualTest) {
  EXPECT_EQ(U, U);
  EXPECT_EQ(V, V);
}

TEST_P(CVectorTest, NonEqualTest) {
  EXPECT_NE(U, V);
  EXPECT_NE(V, U);
}

TEST_P(CVectorTest, SaveLoadTest) {
  Saver saver("SaveLoadTest");
  saver << U;
  saver.close();
  Loader loader("SaveLoadTest");
  CVector Z;
  loader >> Z;
  loader.close();
  EXPECT_EQ(U, Z);
}

TEST_P(CVectorTest, ConjTest) {
  CVector W = conj(U);
  for (int n = 0; n < dim; ++n)
    EXPECT_EQ(W[n], conj(U[n]));
}

TEST_P(CVectorTest, RealTest) {
  RVector W = real(U);
  for (int n = 0; n < dim; ++n)
    EXPECT_EQ(W[n], real(U[n]));
}

TEST_P(CVectorTest, ImagTest) {
  RVector W = imag(U);
  for (int n = 0; n < dim; ++n)
    EXPECT_EQ(W[n], imag(U[n]));
}

INSTANTIATE_TEST_CASE_P(BasicAlgebraTest, CVectorTest,
                        Values(BasicAlgebraVectorTestParameter{1},
                               BasicAlgebraVectorTestParameter{2},
                               BasicAlgebraVectorTestParameter{3},
                               BasicAlgebraVectorTestParameter{4},
                               BasicAlgebraVectorTestParameter{5}));

int main(int argc, char **argv) {
  MppTest mppTest = MppTestBuilder(argc, argv).WithPPM();
  return mppTest.RUN_ALL_MPP_TESTS();
}

// #include "gtest/gtest.h"
//
// #include "basicalgebra/TestBasicAlgebra.hpp"
// #include "basicalgebra/CVector.hpp"
//
//
// class BasicAlgebraCVectorTest : public BasicAlgebraVectorTest {
// protected:
//     vector<double> valuesU_re, valuesV_re, valuesU_im, valuesV_im;
//     double value_re, value_im;
//     int value_int;
//     CVector U, V, U_re_c;
//     RVector U_re, V_re;
//     double expectedNormSqrU;
//     std::complex<double> expectedScalarProductUV, expectedScalarProductU_reV,
//         expectedScalarProductUV_re;
// public:
//     BasicAlgebraCVectorTest() {
//         valuesU_re.resize(dim);
//         valuesU_im.resize(dim);
//         valuesV_re.resize(dim);
//         valuesV_im.resize(dim);
//         U.resize(dim);
//         V.resize(dim);
//         U_re.resize(dim);
//         V_re.resize(dim);
//         U_re_c.resize(dim);
//         expectedNormSqrU = 0.0;
//         expectedScalarProductUV = 0.0;
//         for (int n = 0; n < dim; ++n) {
//             valuesU_re[n] = double(rand()) / RAND_MAX;
//             valuesU_im[n] = double(rand()) / RAND_MAX;
//             valuesV_re[n] = double(rand()) / RAND_MAX;
//             valuesV_im[n] = double(rand()) / RAND_MAX;
//             U[n] = std::complex<double>(valuesU_re[n], valuesU_im[n]);
//             V[n] = std::complex<double>(valuesV_re[n], valuesV_im[n]);
//             U_re[n] = valuesU_re[n];
//             U_re_c[n] = valuesU_re[n];
//             V_re[n] = valuesV_re[n];
//             expectedNormSqrU +=
//                 valuesU_re[n] * valuesU_re[n] + valuesU_im[n] * valuesU_im[n];
//             expectedScalarProductUV += std::complex<double>(
//                 valuesU_re[n] * valuesV_re[n] + valuesU_im[n] * valuesV_im[n],
//                 valuesU_im[n] * valuesV_re[n] - valuesU_re[n] * valuesV_im[n]);
//             expectedScalarProductU_reV +=
//                 std::complex<double>(valuesU_re[n] * valuesV_re[n],
//                                      -valuesU_re[n] * valuesV_im[n]);
//             expectedScalarProductUV_re +=
//                 std::complex<double>(valuesU_re[n] * valuesV_re[n],
//                                      valuesU_im[n] * valuesV_re[n]);
//         }
//         value_re = double(rand()) / RAND_MAX;
//         value_im = double(rand()) / RAND_MAX;
//         value_int = rand();
//     }
// };
//
// TEST_P(BasicAlgebraCVectorTest, ConstructorTest) {
//     CVector W1(std::complex<double>(value_re, value_im), dim);
//     CVector W2(value_int, dim);
//     CVector W3(value_re, dim);
//     CVector W4(U);
//     CVector W5(U_re);
//     CVector W6(dim);
//     for (int n = 0; n < dim; ++n) {
//         EXPECT_DOUBLE_EQ(std::real(U[n]), valuesU_re[n]);
//         EXPECT_DOUBLE_EQ(std::imag(U[n]), valuesU_im[n]);
//         EXPECT_DOUBLE_EQ(std::real(V[n]), valuesV_re[n]);
//         EXPECT_DOUBLE_EQ(std::imag(V[n]), valuesV_im[n]);
//         EXPECT_DOUBLE_EQ(std::real(W1[n]), value_re);
//         EXPECT_DOUBLE_EQ(std::imag(W1[n]), value_im);
//         EXPECT_DOUBLE_EQ(std::real(W2[n]), double(value_int));
//         EXPECT_DOUBLE_EQ(std::imag(W2[n]), 0.0);
//         EXPECT_DOUBLE_EQ(std::real(W3[n]), value_re);
//         EXPECT_DOUBLE_EQ(std::imag(W3[n]), 0.0);
//         EXPECT_DOUBLE_EQ(std::real(W4[n]), valuesU_re[n]);
//         EXPECT_DOUBLE_EQ(std::imag(W4[n]), valuesU_im[n]);
//         EXPECT_DOUBLE_EQ(std::real(W5[n]), valuesU_re[n]);
//         EXPECT_DOUBLE_EQ(std::imag(W5[n]), 0.0);
//         EXPECT_DOUBLE_EQ(std::real(W6[n]), 0.0);
//         EXPECT_DOUBLE_EQ(std::imag(W6[n]), 0.0);
//     }
// }
//
// TEST_P(BasicAlgebraCVectorTest, AssignmentTest) {
//     CVector W1(dim);
//     W1 = std::complex<double>(value_re, value_im);
//     CVector W2(dim);
//     W2 = value_int;
//     CVector W3(dim);
//     W3 = value_re;
//     CVector W4;
//     W4 = U;
//     CVector W5;
//     W5 = U_re;
//     for (int n = 0; n < dim; ++n) {
//         EXPECT_DOUBLE_EQ(std::real(W1[n]), value_re);
//         EXPECT_DOUBLE_EQ(std::imag(W1[n]), value_im);
//         EXPECT_DOUBLE_EQ(std::real(W2[n]), double(value_int));
//         EXPECT_DOUBLE_EQ(std::imag(W2[n]), 0.0);
//         EXPECT_DOUBLE_EQ(std::real(W3[n]), value_re);
//         EXPECT_DOUBLE_EQ(std::imag(W3[n]), 0.0);
//         EXPECT_DOUBLE_EQ(std::real(W4[n]), valuesU_re[n]);
//         EXPECT_DOUBLE_EQ(std::imag(W4[n]), valuesU_im[n]);
//         EXPECT_DOUBLE_EQ(std::real(W5[n]), valuesU_re[n]);
//         EXPECT_DOUBLE_EQ(std::imag(W5[n]), 0.0);
//     }
// }
//
// TEST_P(BasicAlgebraCVectorTest, SummationTest) {
//     CVector UV = U + V;
//     CVector U_reV = U_re + V;
//     CVector UV_re = U + V_re;
//     CVector Uv(U);
//     Uv += std::complex<double>(value_re, value_im);
//     CVector Uv_re(U);
//     Uv_re += value_re;
//     CVector W1(U);
//     W1 += V_re;
//     for (int n = 0; n < dim; ++n) {
//         EXPECT_DOUBLE_EQ(std::real(UV[n]), valuesU_re[n] + valuesV_re[n]);
//         EXPECT_DOUBLE_EQ(std::imag(UV[n]), valuesU_im[n] + valuesV_im[n]);
//         EXPECT_DOUBLE_EQ(std::real(U_reV[n]), valuesU_re[n] + valuesV_re[n]);
//         EXPECT_DOUBLE_EQ(std::imag(U_reV[n]), valuesV_im[n]);
//         EXPECT_DOUBLE_EQ(std::real(UV_re[n]), valuesU_re[n] + valuesV_re[n]);
//         EXPECT_DOUBLE_EQ(std::imag(UV_re[n]), valuesU_im[n]);
//         EXPECT_DOUBLE_EQ(std::real(Uv[n]), valuesU_re[n] + value_re);
//         EXPECT_DOUBLE_EQ(std::imag(Uv[n]), valuesU_im[n] + value_im);
//         EXPECT_DOUBLE_EQ(std::real(Uv_re[n]), valuesU_re[n] + value_re);
//         EXPECT_DOUBLE_EQ(std::imag(Uv_re[n]), valuesU_im[n]);
//         EXPECT_DOUBLE_EQ(std::real(W1[n]), valuesU_re[n] + valuesV_re[n]);
//         EXPECT_DOUBLE_EQ(std::imag(W1[n]), valuesU_im[n]);
//     }
// }
//
// TEST_P(BasicAlgebraCVectorTest, DiffenceTest) {
//     CVector UV = U - V;
//     CVector U_reV = U_re - V;
//     CVector UV_re = U - V_re;
//     CVector Uv(U);
//     Uv -= std::complex<double>(value_re, value_im);
//     CVector Uv_re(U);
//     Uv_re -= value_re;
//     CVector W1(U);
//     W1 -= V_re;
//     for (int n = 0; n < dim; ++n) {
//         EXPECT_DOUBLE_EQ(std::real(UV[n]), valuesU_re[n] - valuesV_re[n]);
//         EXPECT_DOUBLE_EQ(std::imag(UV[n]), valuesU_im[n] - valuesV_im[n]);
//         EXPECT_DOUBLE_EQ(std::real(U_reV[n]), valuesU_re[n] - valuesV_re[n]);
//         EXPECT_DOUBLE_EQ(std::imag(U_reV[n]), -valuesV_im[n]);
//         EXPECT_DOUBLE_EQ(std::real(UV_re[n]), valuesU_re[n] - valuesV_re[n]);
//         EXPECT_DOUBLE_EQ(std::imag(UV_re[n]), valuesU_im[n]);
//         EXPECT_DOUBLE_EQ(std::real(Uv[n]), valuesU_re[n] - value_re);
//         EXPECT_DOUBLE_EQ(std::imag(Uv[n]), valuesU_im[n] - value_im);
//         EXPECT_DOUBLE_EQ(std::real(Uv_re[n]), valuesU_re[n] - value_re);
//         EXPECT_DOUBLE_EQ(std::imag(Uv_re[n]), valuesU_im[n]);
//         EXPECT_DOUBLE_EQ(std::real(W1[n]), valuesU_re[n] - valuesV_re[n]);
//         EXPECT_DOUBLE_EQ(std::imag(W1[n]), valuesU_im[n]);
//     }
// }
//
// TEST_P(BasicAlgebraCVectorTest, AdditiveInversTest) {
//     CVector V = -U;
//     for (int n = 0; n < dim; ++n) {
//         EXPECT_DOUBLE_EQ(std::real(V[n]), -valuesU_re[n]);
//         EXPECT_DOUBLE_EQ(std::imag(V[n]), -valuesU_im[n]);
//     }
// }
//
// TEST_P(BasicAlgebraCVectorTest, ScalarMultiplicationTest) {
//     CVector aU = std::complex<double>(value_re, value_im) * U;
//     CVector Ua = U * std::complex<double>(value_re, value_im);
//     CVector bU = value_re * U;
//     CVector Ub = U * value_re;
//     CVector cU = value_int * U;
//     CVector Uc = U * value_int;
//     for (int n = 0; n < dim; ++n) {
//         EXPECT_DOUBLE_EQ(std::real(aU[n]),
//                          value_re * valuesU_re[n] - value_im * valuesU_im[n]);
//         EXPECT_DOUBLE_EQ(std::imag(aU[n]),
//                          value_re * valuesU_im[n] + value_im * valuesU_re[n]);
//         EXPECT_DOUBLE_EQ(std::real(Ua[n]),
//                          value_re * valuesU_re[n] - value_im * valuesU_im[n]);
//         EXPECT_DOUBLE_EQ(std::imag(Ua[n]),
//                          value_re * valuesU_im[n] + value_im * valuesU_re[n]);
//         EXPECT_DOUBLE_EQ(std::real(bU[n]), value_re * valuesU_re[n]);
//         EXPECT_DOUBLE_EQ(std::imag(bU[n]), value_re * valuesU_im[n]);
//         EXPECT_DOUBLE_EQ(std::real(Ub[n]), valuesU_re[n] * value_re);
//         EXPECT_DOUBLE_EQ(std::imag(Ub[n]), valuesU_im[n] * value_re);
//         EXPECT_DOUBLE_EQ(std::real(cU[n]), value_int * valuesU_re[n]);
//         EXPECT_DOUBLE_EQ(std::imag(cU[n]), value_int * valuesU_im[n]);
//         EXPECT_DOUBLE_EQ(std::real(Uc[n]), valuesU_re[n] * value_int);
//         EXPECT_DOUBLE_EQ(std::imag(Uc[n]), valuesU_im[n] * value_int);
//     }
// }
//
// TEST_P(BasicAlgebraCVectorTest, ScalarDivisionTest) {
//     CVector Ua = U / std::complex<double>(value_re, value_im);
//     CVector Ub = U / value_re;
//     CVector Uc = U / value_int;
//     for (int n = 0; n < dim; ++n) {
//         EXPECT_DOUBLE_EQ(std::real(Ua[n]),
//                          (value_re * valuesU_re[n] + value_im * valuesU_im[n])
//                              / (value_re * value_re + value_im * value_im));
//         EXPECT_DOUBLE_EQ(std::imag(Ua[n]),
//                          (value_re * valuesU_im[n] - value_im * valuesU_re[n])
//                              / (value_re * value_re + value_im * value_im));
//         EXPECT_DOUBLE_EQ(std::real(Ub[n]), valuesU_re[n] / value_re);
//         EXPECT_DOUBLE_EQ(std::imag(Ub[n]), valuesU_im[n] / value_re);
//         EXPECT_DOUBLE_EQ(std::real(Uc[n]), valuesU_re[n] / value_int);
//         EXPECT_DOUBLE_EQ(std::imag(Uc[n]), valuesU_im[n] / value_int);
//     }
// }
//
// TEST_P(BasicAlgebraCVectorTest, ScalarProductTest) {
//     std::complex<double> UV = U * V;
//     std::complex<double> U_reV = U_re * V;
//     std::complex<double> UV_re = U * V_re;
//     EXPECT_DOUBLE_EQ(std::real(UV), std::real(expectedScalarProductUV));
//     EXPECT_DOUBLE_EQ(std::imag(UV), std::imag(expectedScalarProductUV));
//     EXPECT_DOUBLE_EQ(std::real(U_reV), std::real(expectedScalarProductU_reV));
//     EXPECT_DOUBLE_EQ(std::imag(U_reV), std::imag(expectedScalarProductU_reV));
//     EXPECT_DOUBLE_EQ(std::real(UV_re), std::real(expectedScalarProductUV_re));
//     EXPECT_DOUBLE_EQ(std::imag(UV_re), std::imag(expectedScalarProductUV_re));
// }
//
// TEST_P(BasicAlgebraCVectorTest, NormSqrTest) {
//     double n = normSqr(U);
//     EXPECT_DOUBLE_EQ(n, expectedNormSqrU);
// }
//
// TEST_P(BasicAlgebraCVectorTest, NormTest) {
//     double n = norm(U);
//     EXPECT_DOUBLE_EQ(n, sqrt(expectedNormSqrU));
// }
//
// TEST_P(BasicAlgebraCVectorTest, ConjTest) {
//     CVector W = conj(U);
//     for (int n = 0; n < dim; ++n) {
//         EXPECT_DOUBLE_EQ(std::real(W[n]), valuesU_re[n]);
//         EXPECT_DOUBLE_EQ(std::imag(W[n]), -valuesU_im[n]);
//     }
// }
//
// TEST_P(BasicAlgebraCVectorTest, RealPartTest) {
//     RVector W = real(U);
//     for (int n = 0; n < dim; ++n)
//         EXPECT_DOUBLE_EQ(W[n], valuesU_re[n]);
// }
//
// TEST_P(BasicAlgebraCVectorTest, ImaginaryPartTest) {
//     RVector W = imag(U);
//     for (int n = 0; n < dim; ++n)
//         EXPECT_DOUBLE_EQ(W[n], valuesU_im[n]);
// }
//
// TEST_P(BasicAlgebraCVectorTest, EqualTest) {
//     EXPECT_EQ(U, U);
//     EXPECT_EQ(V, V);
//     EXPECT_EQ(U_re_c, U_re);
//     EXPECT_EQ(U_re, U_re_c);
// }
//
// TEST_P(BasicAlgebraCVectorTest, NonEqualTest) {
//     EXPECT_NE(U, V);
//     EXPECT_NE(V, U);
//     EXPECT_NE(V_re, U_re_c);
//     EXPECT_NE(U_re_c, V_re);
// }
//
// TEST_P(BasicAlgebraCVectorTest, SaveLoadTest) {
//     Saver saver("SaveLoadTest");
//     saver << U;
//     saver.close();
//     Loader loader("SaveLoadTest");
//     CVector Z;
//     loader >> Z;
//     loader.close();
//     for (int n = 0; n < this->dim; ++n) {
//         EXPECT_DOUBLE_EQ(std::real(Z[n]), valuesU_re[n]);
//         EXPECT_DOUBLE_EQ(std::imag(Z[n]), valuesU_im[n]);
//     }
// }
//
// INSTANTIATE_TEST_CASE_P(BasicAlgebraVectorTest, BasicAlgebraCVectorTest, Values(
//     BasicAlgebraVectorTestParameter{1},
//     BasicAlgebraVectorTestParameter{2},
//     BasicAlgebraVectorTestParameter{3},
//     BasicAlgebraVectorTestParameter{4},
//     BasicAlgebraVectorTestParameter{5}
//));
//
// int main(int argc, char **argv) {
//     MppTest mppTest = MppTestBuilder(argc, argv).WithPPM();
//     return mppTest.RUN_ALL_MPP_TESTS();
// }
