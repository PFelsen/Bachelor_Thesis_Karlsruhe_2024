#include "gtest/gtest.h"

#include "TestBasicAlgebra.hpp"
#include "basicalgebra/RVector.hpp"

class RVectorTest : public BasicAlgebraVectorTest {
protected:
  RVector U, V;
  double value;
  int value_i;
public:
  RVectorTest() : BasicAlgebraVectorTest() {
    U.resize(dim);
    V.resize(dim);

    for (int n = 0; n < dim; ++n) {
      U[n] = RandomDouble();
      V[n] = RandomDouble();
    }
    value = RandomNonZeroDouble();
    value_i = RandomNonZeroInt();
  }
};

TEST_P(RVectorTest, ConstructorTest) {
  RVector W0(dim);
  RVector W1(value, dim);
  RVector W2(value_i, dim);
  RVector W3(U);

  for (int n = 0; n < dim; ++n) {
    EXPECT_DOUBLE_EQ(W0[n], 0.0);
    EXPECT_DOUBLE_EQ(W1[n], value);
    EXPECT_DOUBLE_EQ(W2[n], value_i);
    EXPECT_DOUBLE_EQ(W3[n], U[n]);
  }
}

TEST_P(RVectorTest, AssignmentTest) {
  RVector Z = U;
  EXPECT_EQ(Z, U);
  RVector W(dim);
  W = value;
  EXPECT_EQ(W, RVector(value, dim));
  W = value_i;
  EXPECT_EQ(W, RVector(value_i, dim));
}

TEST_P(RVectorTest, SummationTest) {
  RVector UV = U + V;
  RVector W1(U);
  W1 += value;
  RVector W2(U);
  W2 += value_i;
  for (int n = 0; n < dim; ++n) {
    EXPECT_DOUBLE_EQ(UV[n], U[n] + V[n]);
    EXPECT_DOUBLE_EQ(W1[n], U[n] + value);
    EXPECT_DOUBLE_EQ(W2[n], U[n] + value_i);
  }
}

TEST_P(RVectorTest, DifferenceTest) {
  RVector UV = U - V;
  RVector W1(U);
  W1 -= value;
  RVector W2(U);
  W2 -= value_i;
  for (int n = 0; n < dim; ++n) {
    EXPECT_DOUBLE_EQ(UV[n], U[n] - V[n]);
    EXPECT_DOUBLE_EQ(W1[n], U[n] - value);
    EXPECT_DOUBLE_EQ(W2[n], U[n] - value_i);
  }
}

TEST_P(RVectorTest, AdditiveInversTest) {
  RVector V = -U;
  for (int n = 0; n < dim; ++n)
    EXPECT_DOUBLE_EQ(V[n], -U[n]);
}

TEST_P(RVectorTest, ScalarMultiplicationTest) {
  RVector W1 = value * U;
  RVector W2 = U * value;
  RVector W3 = value_i * U;
  RVector W4 = U * value_i;
  for (int n = 0; n < dim; ++n) {
    EXPECT_DOUBLE_EQ(W1[n], value * U[n]);
    EXPECT_DOUBLE_EQ(W2[n], value * U[n]);
    EXPECT_DOUBLE_EQ(W3[n], value_i * U[n]);
    EXPECT_DOUBLE_EQ(W4[n], value_i * U[n]);
  }
}

TEST_P(RVectorTest, ScalarDivisionTest) {
  RVector W1 = U / value;
  RVector W2 = U / value_i;
  for (int n = 0; n < dim; ++n) {
    EXPECT_DOUBLE_EQ(W1[n], U[n] / value);
    EXPECT_DOUBLE_EQ(W2[n], U[n] / value_i);
  }
}

TEST_P(RVectorTest, ScalarProductTest) {
  double a0 = 0.0;
  for (int n = 0; n < dim; ++n)
    a0 += U[n] * V[n];
  EXPECT_DOUBLE_EQ(a0, U * V);
}

TEST_P(RVectorTest, NormSqrTest) {
  double sum = 0.0;
  for (int n = 0; n < dim; ++n)
    sum += U[n] * U[n];
  EXPECT_DOUBLE_EQ(normSqr(U), sum);
}

TEST_P(RVectorTest, NormTest) {
  double sum = 0.0;
  for (int n = 0; n < dim; ++n)
    sum += U[n] * U[n];
  EXPECT_DOUBLE_EQ(norm(U), sqrt(sum));
}

TEST_P(RVectorTest, EqualTest) {
  EXPECT_EQ(U, U);
  EXPECT_EQ(V, V);
}

TEST_P(RVectorTest, NonEqualTest) {
  EXPECT_NE(U, V);
  EXPECT_NE(V, U);
}

TEST_P(RVectorTest, SaveLoadTest) {
  Saver saver("SaveLoadTest");
  saver << U;
  saver.close();
  Loader loader("SaveLoadTest");
  RVector Z;
  loader >> Z;
  loader.close();
  EXPECT_EQ(U, Z);
}

INSTANTIATE_TEST_CASE_P(BasicAlgebraTest, RVectorTest,
                        Values(BasicAlgebraVectorTestParameter{1},
                               BasicAlgebraVectorTestParameter{2},
                               BasicAlgebraVectorTestParameter{3},
                               BasicAlgebraVectorTestParameter{4},
                               BasicAlgebraVectorTestParameter{5}));

int main(int argc, char **argv) {
  MppTest mppTest = MppTestBuilder(argc, argv).WithPPM();
  return mppTest.RUN_ALL_MPP_TESTS();
}
