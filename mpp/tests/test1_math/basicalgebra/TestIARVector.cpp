#include "gtest/gtest.h"

#include "RVector.hpp"
#include "TestIABasicAlgebra.hpp"

class IARVectorTest : public IABasicAlgebraVectorTest {
protected:
  IARVector U, V;
  IAInterval value;
  double value_d;
  int value_i;
public:
  IARVectorTest() : IABasicAlgebraVectorTest() {
    U.resize(dim);
    V.resize(dim);

    for (int n = 0; n < dim; ++n) {
      U[n] = RandomIAInterval();
      V[n] = RandomIAInterval();
    }
    value = RandomNonZeroIAInterval();
    value_d = RandomNonZeroDouble();
    value_i = RandomNonZeroInt();
  }
};

TEST_P(IARVectorTest, ConstructorTest) {
  IARVector W0(dim);
  IARVector W1(value, dim);
  IARVector W2(value_d, dim);
  IARVector W3(value_i, dim);
  IARVector W4(U);
  IARVector W5(mid(U));

  for (int n = 0; n < dim; ++n) {
    EXPECT_EQ(W0[n], IAInterval());
    EXPECT_EQ(W1[n], value);
    EXPECT_EQ(W2[n], IAInterval(value_d));
    EXPECT_EQ(W3[n], IAInterval(value_i));
    EXPECT_EQ(W4[n], U[n]);
    EXPECT_EQ(W5[n], mid(U[n]));
  }
}

TEST_P(IARVectorTest, AssignmentTest) {
  IARVector Z0 = U;
  EXPECT_EQ(Z0, U);
  IARVector Z1;
  Z1 = mid(U);
  EXPECT_EQ(Z1, IARVector(mid(U)));
  IARVector W(dim);
  W = value;
  EXPECT_EQ(W, IARVector(value, dim));
  W = value_d;
  EXPECT_EQ(W, IARVector(value_d, dim));
  W = value_i;
  EXPECT_EQ(W, IARVector(value_i, dim));
}

TEST_P(IARVectorTest, SummationTest) {
  IARVector UV = U + V;
  IARVector W1 = U + mid(V);
  IARVector W2 = mid(U) + V;
  IARVector W3(U);
  W3 += value;
  IARVector W4(U);
  W4 += value_d;
  IARVector W5(U);
  W5 += value_i;
  for (int n = 0; n < dim; ++n) {
    EXPECT_EQ(UV[n], U[n] + V[n]);
    EXPECT_EQ(W1[n], U[n] + mid(V[n]));
    EXPECT_EQ(W2[n], mid(U[n]) + V[n]);
    EXPECT_EQ(W3[n], U[n] + value);
    EXPECT_EQ(W4[n], U[n] + IAInterval(value_d));
    EXPECT_EQ(W5[n], U[n] + IAInterval(value_i));
  }
}

TEST_P(IARVectorTest, DifferenceTest) {
  IARVector UV = U - V;
  IARVector W1 = U - mid(V);
  IARVector W2 = mid(U) - V;
  IARVector W3(U);
  W3 -= value;
  IARVector W4(U);
  W4 -= value_d;
  IARVector W5(U);
  W5 -= value_i;
  for (int n = 0; n < dim; ++n) {
    EXPECT_EQ(UV[n], U[n] - V[n]);
    EXPECT_EQ(W1[n], U[n] - mid(V[n]));
    EXPECT_EQ(W2[n], mid(U[n]) - V[n]);
    EXPECT_EQ(W3[n], U[n] - value);
    EXPECT_EQ(W4[n], U[n] - IAInterval(value_d));
    EXPECT_EQ(W5[n], U[n] - IAInterval(value_i));
  }
}

TEST_P(IARVectorTest, AdditiveInversTest) {
  IARVector V = -U;
  for (int n = 0; n < dim; ++n)
    EXPECT_EQ(V[n], -U[n]);
}

TEST_P(IARVectorTest, ScalarMultiplicationTest) {
  IARVector W1 = value * U;
  IARVector W2 = U * value;
  IARVector W3 = value_d * U;
  IARVector W4 = U * value_d;
  IARVector W5 = value_i * U;
  IARVector W6 = U * value_i;
  for (int n = 0; n < dim; ++n) {
    EXPECT_EQ(W1[n], value * U[n]);
    EXPECT_EQ(W2[n], value * U[n]);
    EXPECT_EQ(W3[n], value_d * U[n]);
    EXPECT_EQ(W4[n], value_d * U[n]);
    EXPECT_EQ(W5[n], value_i * U[n]);
    EXPECT_EQ(W6[n], value_i * U[n]);
  }
}

TEST_P(IARVectorTest, ScalarDivisionTest) {
  IARVector W1 = U / value;
  IARVector W2 = U / value_d;
  IARVector W3 = U / value_i;
  for (int n = 0; n < dim; ++n) {
    EXPECT_EQ(W1[n], U[n] / value);
    EXPECT_EQ(W2[n], U[n] / value_d);
    EXPECT_EQ(W3[n], U[n] / value_i);
  }
}

TEST_P(IARVectorTest, ScalarProductTest) {
  IAInterval a0, a1, a2;
  for (int n = 0; n < dim; ++n) {
    a0 += U[n] * V[n];
    a1 += U[n] * mid(V[n]);
    a2 += mid(U[n]) * V[n];
  }
  EXPECT_EQ(a0, U * V);
  EXPECT_EQ(a1, U * mid(V));
  EXPECT_EQ(a2, mid(U) * V);
}

TEST_P(IARVectorTest, NormSqrTest) {
  IAInterval sum;
  for (int n = 0; n < dim; ++n)
    sum += sqr(U[n]);
  EXPECT_EQ(normSqr(U), sum);
}

TEST_P(IARVectorTest, NormTest) {
  IAInterval sum;
  for (int n = 0; n < dim; ++n)
    sum += sqr(U[n]);
  EXPECT_EQ(norm(U), sqrt(sum));
}

TEST_P(IARVectorTest, EqualTest) {
  EXPECT_EQ(U, U);
  EXPECT_EQ(V, V);
}

TEST_P(IARVectorTest, NonEqualTest) {
  EXPECT_NE(U, V);
  EXPECT_NE(V, U);
}

TEST_P(IARVectorTest, SaveLoadTest) {
  Saver saver("SaveLoadTest");
  saver << U;
  saver.close();
  Loader loader("SaveLoadTest");
  IARVector Z;
  loader >> Z;
  loader.close();
  EXPECT_EQ(U, Z);
}

TEST_P(IARVectorTest, InfTest) {
  RVector W = inf(U);
  for (int n = 0; n < dim; ++n)
    EXPECT_EQ(W[n], inf(U[n]));
}

TEST_P(IARVectorTest, SupTest) {
  RVector W = sup(U);
  for (int n = 0; n < dim; ++n)
    EXPECT_EQ(W[n], sup(U[n]));
}

TEST_P(IARVectorTest, MidTest) {
  RVector W = mid(U);
  for (int n = 0; n < dim; ++n)
    EXPECT_EQ(W[n], mid(U[n]));
}

INSTANTIATE_TEST_CASE_P(IABasicAlgebraTest, IARVectorTest,
                        Values(IABasicAlgebraVectorTestParameter{1},
                               IABasicAlgebraVectorTestParameter{2},
                               IABasicAlgebraVectorTestParameter{3},
                               IABasicAlgebraVectorTestParameter{4},
                               IABasicAlgebraVectorTestParameter{5}));

int main(int argc, char **argv) {
  MppTest mppTest = MppTestBuilder(argc, argv).WithPPM();
  return mppTest.RUN_ALL_MPP_TESTS();
}
