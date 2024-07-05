#include "gtest/gtest.h"

#include "CVector.hpp"
#include "TestIABasicAlgebra.hpp"

class IACVectorTest : public IABasicAlgebraVectorTest {
protected:
  IACVector U, V;
  IARVector U_ia;
  IACInterval value;
  IAInterval value_ia;
  std::complex<double> value_c;
  double value_d;
  int value_i;
public:
  IACVectorTest() : IABasicAlgebraVectorTest() {
    U.resize(dim);
    V.resize(dim);
    U_ia.resize(dim);

    for (int n = 0; n < dim; ++n) {
      U[n] = RandomIACInterval();
      V[n] = RandomIACInterval();
      U_ia[n] = RandomIAInterval();
    }
    value = RandomNonZeroIACInterval();
    value_ia = RandomNonZeroIAInterval();
    value_c = RandomNonZeroComplex();
    value_d = RandomNonZeroDouble();
    value_i = RandomNonZeroInt();
  }
};

TEST_P(IACVectorTest, ConstructorTest) {
  IACVector W0(dim);
  IACVector W1(value, dim);
  IACVector W2(value_ia, dim);
  IACVector W3(value_c, dim);
  IACVector W4(value_d, dim);
  IACVector W5(value_i, dim);
  IACVector W6(U);
  IACVector W7(mid(U));
  IACVector W8(U_ia);
  IACVector W9(mid(U_ia));

  for (int n = 0; n < dim; ++n) {
    EXPECT_EQ(W0[n], IACInterval());
    EXPECT_EQ(W1[n], value);
    EXPECT_EQ(W2[n], IACInterval(value_ia));
    EXPECT_EQ(W3[n], IACInterval(value_c));
    EXPECT_EQ(W4[n], IACInterval(value_d));
    EXPECT_EQ(W5[n], IACInterval(value_i));
    EXPECT_EQ(W6[n], U[n]);
    EXPECT_EQ(W7[n], mid(U[n]));
    EXPECT_EQ(W8[n], U_ia[n]);
    EXPECT_EQ(W9[n], mid(U_ia[n]));
  }
}

TEST_P(IACVectorTest, AssignmentTest) {
  IACVector Z0 = U;
  EXPECT_EQ(Z0, U);
  IACVector Z1;
  Z1 = mid(U);
  EXPECT_EQ(Z1, IACVector(mid(U)));
  IACVector Z2;
  Z2 = U_ia;
  EXPECT_EQ(Z2, IACVector(U_ia));
  IACVector Z3;
  Z3 = mid(U_ia);
  EXPECT_EQ(Z3, IACVector(mid(U_ia)));
  IACVector W(dim);
  W = value;
  EXPECT_EQ(W, IACVector(value, dim));
  W = value_ia;
  EXPECT_EQ(W, IACVector(value_ia, dim));
  W = value_c;
  EXPECT_EQ(W, IACVector(value_c, dim));
  W = value_d;
  EXPECT_EQ(W, IACVector(value_d, dim));
  W = value_i;
  EXPECT_EQ(W, IACVector(value_i, dim));
}

TEST_P(IACVectorTest, SummationTest) {
  IACVector UV = U + V;
  IACVector W1 = U + mid(V);
  IACVector W2 = mid(U) + V;
  IACVector W3 = U + U_ia;
  IACVector W4 = U_ia + V;
  IACVector W5 = U + mid(U_ia);
  IACVector W6 = mid(U_ia) + V;
  IACVector W7(U);
  W7 += value;
  IACVector W8(U);
  W8 += value_ia;
  IACVector W9(U);
  W9 += value_c;
  IACVector W10(U);
  W10 += value_d;
  IACVector W11(U);
  W11 += value_i;
  for (int n = 0; n < dim; ++n) {
    EXPECT_EQ(UV[n], U[n] + V[n]);
    EXPECT_EQ(W1[n], U[n] + mid(V[n]));
    EXPECT_EQ(W2[n], mid(U[n]) + V[n]);
    EXPECT_EQ(W3[n], U[n] + U_ia[n]);
    EXPECT_EQ(W4[n], U_ia[n] + V[n]);
    EXPECT_EQ(W5[n], U[n] + mid(U_ia[n]));
    EXPECT_EQ(W6[n], mid(U_ia[n]) + V[n]);
    EXPECT_EQ(W7[n], U[n] + value);
    EXPECT_EQ(W8[n], U[n] + IACInterval(value_ia));
    EXPECT_EQ(W9[n], U[n] + IACInterval(value_c));
    EXPECT_EQ(W10[n], U[n] + IACInterval(value_d));
    EXPECT_EQ(W11[n], U[n] + IACInterval(value_i));
  }
}

TEST_P(IACVectorTest, DifferenceTest) {
  IACVector UV = U - V;
  IACVector W1 = U - mid(V);
  IACVector W2 = mid(U) - V;
  IACVector W3 = U - U_ia;
  IACVector W4 = U_ia - V;
  IACVector W5 = U - mid(U_ia);
  IACVector W6 = mid(U_ia) - V;
  IACVector W7(U);
  W7 -= value;
  IACVector W8(U);
  W8 -= value_ia;
  IACVector W9(U);
  W9 -= value_c;
  IACVector W10(U);
  W10 -= value_d;
  IACVector W11(U);
  W11 -= value_i;
  for (int n = 0; n < dim; ++n) {
    EXPECT_EQ(UV[n], U[n] - V[n]);
    EXPECT_EQ(W1[n], U[n] - mid(V[n]));
    EXPECT_EQ(W2[n], mid(U[n]) - V[n]);
    EXPECT_EQ(W3[n], U[n] - U_ia[n]);
    EXPECT_EQ(W4[n], U_ia[n] - V[n]);
    EXPECT_EQ(W5[n], U[n] - mid(U_ia[n]));
    EXPECT_EQ(W6[n], mid(U_ia[n]) - V[n]);
    EXPECT_EQ(W7[n], U[n] - value);
    EXPECT_EQ(W8[n], U[n] - IACInterval(value_ia));
    EXPECT_EQ(W9[n], U[n] - IACInterval(value_c));
    EXPECT_EQ(W10[n], U[n] - IACInterval(value_d));
    EXPECT_EQ(W11[n], U[n] - IACInterval(value_i));
  }
}

TEST_P(IACVectorTest, AdditiveInversTest) {
  IACVector V = -U;
  for (int n = 0; n < dim; ++n)
    EXPECT_EQ(V[n], -U[n]);
}

TEST_P(IACVectorTest, ScalarMultiplicationTest) {
  IACVector W1 = value * U;
  IACVector W2 = U * value;
  IACVector W3 = value_ia * U;
  IACVector W4 = U * value_ia;
  IACVector W5 = value_c * U;
  IACVector W6 = U * value_c;
  IACVector W7 = value_d * U;
  IACVector W8 = U * value_d;
  IACVector W9 = value_i * U;
  IACVector W10 = U * value_i;
  for (int n = 0; n < dim; ++n) {
    EXPECT_EQ(W1[n], value * U[n]);
    EXPECT_EQ(W2[n], value * U[n]);
    EXPECT_EQ(W3[n], value_ia * U[n]);
    EXPECT_EQ(W4[n], value_ia * U[n]);
    EXPECT_EQ(W5[n], value_c * U[n]);
    EXPECT_EQ(W6[n], value_c * U[n]);
    EXPECT_EQ(W7[n], value_d * U[n]);
    EXPECT_EQ(W8[n], value_d * U[n]);
    EXPECT_EQ(W9[n], value_i * U[n]);
    EXPECT_EQ(W10[n], value_i * U[n]);
  }
}

TEST_P(IACVectorTest, ScalarDivisionTest) {
  IACVector W1 = U / value;
  IACVector W2 = U / value_ia;
  IACVector W3 = U / value_c;
  IACVector W4 = U / value_d;
  IACVector W5 = U / value_i;
  for (int n = 0; n < dim; ++n) {
    EXPECT_EQ(W1[n], U[n] / value);
    EXPECT_EQ(W2[n], U[n] / value_ia);
    EXPECT_EQ(W3[n], U[n] / value_c);
    EXPECT_EQ(W4[n], U[n] / value_d);
    EXPECT_EQ(W5[n], U[n] / value_i);
  }
}

TEST_P(IACVectorTest, ScalarProductTest) {
  IACInterval a0, a1, a2, a3, a4, a5, a6;
  for (int n = 0; n < dim; ++n) {
    a0 += U[n] * conj(V[n]);
    a1 += U[n] * conj(mid(V[n]));
    a2 += mid(U[n]) * conj(V[n]);
    a3 += U[n] * U_ia[n];
    a4 += U_ia[n] * conj(V[n]);
    a5 += U[n] * mid(U_ia[n]);
    a6 += mid(U_ia[n]) * conj(V[n]);
  }
  EXPECT_EQ(a0, U * V);
  EXPECT_EQ(a1, U * mid(V));
  EXPECT_EQ(a2, mid(U) * V);
  EXPECT_EQ(a3, U * U_ia);
  EXPECT_EQ(a4, U_ia * V);
  EXPECT_EQ(a5, U * mid(U_ia));
  EXPECT_EQ(a6, mid(U_ia) * V);
}

TEST_P(IACVectorTest, NormSqrTest) {
  IAInterval sum;
  for (int n = 0; n < dim; ++n)
    sum += absSqr(U[n]);
  EXPECT_EQ(normSqr(U), sum);
}

TEST_P(IACVectorTest, NormTest) {
  IAInterval sum;
  for (int n = 0; n < dim; ++n)
    sum += absSqr(U[n]);
  EXPECT_EQ(norm(U), sqrt(sum));
}

TEST_P(IACVectorTest, EqualTest) {
  EXPECT_EQ(U, U);
  EXPECT_EQ(V, V);
}

TEST_P(IACVectorTest, NonEqualTest) {
  EXPECT_NE(U, V);
  EXPECT_NE(V, U);
}

TEST_P(IACVectorTest, SaveLoadTest) {
  Saver saver("SaveLoadTest");
  saver << U;
  saver.close();
  Loader loader("SaveLoadTest");
  IACVector Z;
  loader >> Z;
  loader.close();
  EXPECT_EQ(U, Z);
}

TEST_P(IACVectorTest, ConjTest) {
  IACVector W = conj(U);
  for (int n = 0; n < dim; ++n)
    EXPECT_EQ(W[n], conj(U[n]));
}

TEST_P(IACVectorTest, RealTest) {
  IARVector W = real(U);
  for (int n = 0; n < dim; ++n)
    EXPECT_EQ(W[n], real(U[n]));
}

TEST_P(IACVectorTest, ImagTest) {
  IARVector W = imag(U);
  for (int n = 0; n < dim; ++n)
    EXPECT_EQ(W[n], imag(U[n]));
}

TEST_P(IACVectorTest, MidTest) {
  CVector W = mid(U);
  for (int n = 0; n < dim; ++n)
    EXPECT_EQ(W[n], mid(U[n]));
}

INSTANTIATE_TEST_CASE_P(IABasicAlgebraTest, IACVectorTest,
                        Values(IABasicAlgebraVectorTestParameter{1},
                               IABasicAlgebraVectorTestParameter{2},
                               IABasicAlgebraVectorTestParameter{3},
                               IABasicAlgebraVectorTestParameter{4},
                               IABasicAlgebraVectorTestParameter{5}));

int main(int argc, char **argv) {
  MppTest mppTest = MppTestBuilder(argc, argv).WithPPM();
  return mppTest.RUN_ALL_MPP_TESTS();
}
