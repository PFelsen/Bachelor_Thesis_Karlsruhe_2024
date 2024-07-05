#include "gtest/gtest.h"

#include "SymRMatrix.hpp"
#include "TestIABasicAlgebra.hpp"

class IASymRMatrixTest : public IABasicAlgebraMatrixTest {
protected:
  IASymRMatrix A, B;
  IARVector U;
  IAInterval value;
  double value_d;
  int value_i;

  int dim;
public:
  IASymRMatrixTest() : IABasicAlgebraMatrixTest() {
    dim = cols = rows;
    A.resize(dim);
    B.resize(dim);
    U.resize(dim);

    for (int n = 0; n < dim; ++n) {
      U[n] = RandomIAInterval();
      for (int m = n; m < dim; ++m) {
        A(RandomIAInterval(), n, m);
        B(RandomIAInterval(), n, m);
      }
    }
    value = RandomNonZeroIAInterval();
    value_d = RandomNonZeroDouble();
    value_i = RandomNonZeroInt();
  }
};

TEST_P(IASymRMatrixTest, ConstructorTest) {
  IASymRMatrix W0(dim);
  IASymRMatrix W1(value, dim);
  IASymRMatrix W2(value_d, dim);
  IASymRMatrix W3(value_i, dim);
  IASymRMatrix W4(A);
  IASymRMatrix W5(mid(A));

  for (int n = 0; n < dim; ++n)
    for (int m = 0; m < dim; ++m) {
      EXPECT_EQ(W0(n, m), IAInterval());
      EXPECT_EQ(W1(n, m), value);
      EXPECT_EQ(W2(n, m), IAInterval(value_d));
      EXPECT_EQ(W3(n, m), IAInterval(value_i));
      EXPECT_EQ(W4(n, m), A(n, m));
      EXPECT_EQ(W5(n, m), mid(A(n, m)));
    }
}

TEST_P(IASymRMatrixTest, AssignmentTest) {
  IASymRMatrix Z0 = A;
  EXPECT_EQ(Z0, A);
  IASymRMatrix Z1;
  Z1 = mid(A);
  EXPECT_EQ(Z1, IASymRMatrix(mid(A)));
  IASymRMatrix W(dim);
  W = value;
  EXPECT_EQ(W, IASymRMatrix(value, dim));
  W = value_d;
  EXPECT_EQ(W, IASymRMatrix(value_d, dim));
  W = value_i;
  EXPECT_EQ(W, IASymRMatrix(value_i, dim));
}

TEST_P(IASymRMatrixTest, SummationTest) {
  IASymRMatrix AB = A + B;
  IASymRMatrix W1 = A + mid(B);
  IASymRMatrix W2 = mid(A) + B;
  for (int n = 0; n < dim; ++n)
    for (int m = 0; m < dim; ++m) {
      EXPECT_EQ(AB(n, m), A(n, m) + B(n, m));
      EXPECT_EQ(W1(n, m), A(n, m) + mid(B(n, m)));
      EXPECT_EQ(W2(n, m), mid(A(n, m)) + B(n, m));
    }
}

TEST_P(IASymRMatrixTest, DifferenceTest) {
  IASymRMatrix AB = A - B;
  IASymRMatrix W1 = A - mid(B);
  IASymRMatrix W2 = mid(A) - B;
  for (int n = 0; n < dim; ++n)
    for (int m = 0; m < dim; ++m) {
      EXPECT_EQ(AB(n, m), A(n, m) - B(n, m));
      EXPECT_EQ(W1(n, m), A(n, m) - mid(B(n, m)));
      EXPECT_EQ(W2(n, m), mid(A(n, m)) - B(n, m));
    }
}

TEST_P(IASymRMatrixTest, AdditiveInversTest) {
  IASymRMatrix W = -A;
  for (int n = 0; n < dim; ++n)
    for (int m = 0; m < dim; ++m)
      EXPECT_EQ(W(n, m), -A(n, m));
}

TEST_P(IASymRMatrixTest, ScalarMultiplicationTest) {
  IASymRMatrix W1 = value * A;
  IASymRMatrix W2 = A * value;
  IASymRMatrix W3 = value_d * A;
  IASymRMatrix W4 = A * value_d;
  IASymRMatrix W5 = value_i * A;
  IASymRMatrix W6 = A * value_i;
  for (int n = 0; n < dim; ++n)
    for (int m = 0; m < dim; ++m) {
      EXPECT_EQ(W1(n, m), value * A(n, m));
      EXPECT_EQ(W2(n, m), value * A(n, m));
      EXPECT_EQ(W3(n, m), value_d * A(n, m));
      EXPECT_EQ(W4(n, m), value_d * A(n, m));
      EXPECT_EQ(W5(n, m), value_i * A(n, m));
      EXPECT_EQ(W6(n, m), value_i * A(n, m));
    }
}

TEST_P(IASymRMatrixTest, ScalarDivisionTest) {
  IASymRMatrix W1 = A / value;
  IASymRMatrix W2 = A / value_d;
  IASymRMatrix W3 = A / value_i;
  for (int n = 0; n < dim; ++n)
    for (int m = 0; m < dim; ++m) {
      EXPECT_EQ(W1(n, m), A(n, m) / value);
      EXPECT_EQ(W2(n, m), A(n, m) / value_d);
      EXPECT_EQ(W3(n, m), A(n, m) / value_i);
    }
}

TEST_P(IASymRMatrixTest, MatrixVectorProductTest) {
  IARVector W0(dim);
  IARVector W1(dim);
  IARVector W2(dim);
  for (int n = 0; n < dim; ++n)
    for (int k = 0; k < dim; ++k) {
      W0[n] += A(n, k) * U[k];
      W1[n] += A(n, k) * mid(U[k]);
      W2[n] += mid(A(n, k)) * U[k];
    }
  EXPECT_EQ(W0, A * U);
  EXPECT_EQ(W1, A * mid(U));
  EXPECT_EQ(W2, mid(A) * U);
}

TEST_P(IASymRMatrixTest, EqualTest) {
  EXPECT_EQ(A, A);
  EXPECT_EQ(B, B);
}

TEST_P(IASymRMatrixTest, NonEqualTest) {
  EXPECT_NE(A, B);
  EXPECT_NE(B, A);
}

TEST_P(IASymRMatrixTest, SaveLoadTest) {
  Saver saver("SaveLoadTest");
  saver << A;
  saver.close();
  Loader loader("SaveLoadTest");
  IASymRMatrix Z;
  loader >> Z;
  loader.close();
  EXPECT_EQ(A, Z);
}

TEST_P(IASymRMatrixTest, InfTest) {
  SymRMatrix W = inf(A);
  for (int n = 0; n < dim; ++n)
    for (int m = 0; m < dim; ++m)
      EXPECT_EQ(W(n, m), inf(A(n, m)));
}

TEST_P(IASymRMatrixTest, SupTest) {
  SymRMatrix W = sup(A);
  for (int n = 0; n < dim; ++n)
    for (int m = 0; m < dim; ++m)
      EXPECT_EQ(W(n, m), sup(A(n, m)));
}

TEST_P(IASymRMatrixTest, MidTest) {
  SymRMatrix W = mid(A);
  for (int n = 0; n < dim; ++n)
    for (int m = 0; m < dim; ++m)
      EXPECT_EQ(W(n, m), mid(A(n, m)));
}

INSTANTIATE_TEST_CASE_P(IABasicAlgebraMatrixTest, IASymRMatrixTest,
                        Values(IABasicAlgebraMatrixTestParameter{1, 1},
                               IABasicAlgebraMatrixTestParameter{2, 2},
                               IABasicAlgebraMatrixTestParameter{3, 3},
                               IABasicAlgebraMatrixTestParameter{4, 4},
                               IABasicAlgebraMatrixTestParameter{5, 5}));

int main(int argc, char **argv) {
  MppTest mppTest = MppTestBuilder(argc, argv).WithPPM();
  return mppTest.RUN_ALL_MPP_TESTS();
}