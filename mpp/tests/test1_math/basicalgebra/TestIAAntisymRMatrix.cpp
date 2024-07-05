#include "gtest/gtest.h"

#include "RMatrix.hpp"
#include "TestIABasicAlgebra.hpp"

class IAAntisymRMatrixTest : public IABasicAlgebraMatrixTest {
protected:
  IAAntisymRMatrix A, B;
  IARVector U;
  IAInterval value;
  double value_d;
  int value_i;

  int dim;
public:
  IAAntisymRMatrixTest() : IABasicAlgebraMatrixTest() {
    dim = cols = rows;
    A.resize(dim);
    B.resize(dim);
    U.resize(dim);

    for (int n = 0; n < dim; ++n) {
      U[n] = RandomIAInterval();
      for (int m = n + 1; m < dim; ++m) {
        A(RandomIAInterval(), n, m);
        B(RandomIAInterval(), n, m);
      }
    }
    value = RandomNonZeroIAInterval();
    value_d = RandomNonZeroDouble();
    value_i = RandomNonZeroInt();
  }
};

TEST_P(IAAntisymRMatrixTest, ConstructorTest) {
  IAAntisymRMatrix W0(dim);
  IAAntisymRMatrix W1(A);
  IAAntisymRMatrix W2(mid(A));

  for (int n = 0; n < dim; ++n)
    for (int m = 0; m < dim; ++m) {
      EXPECT_EQ(W0(n, m), IAInterval());
      EXPECT_EQ(W1(n, m), A(n, m));
      EXPECT_EQ(W2(n, m), mid(A(n, m)));
    }
}

TEST_P(IAAntisymRMatrixTest, AssignmentTest) {
  IAAntisymRMatrix Z0 = A;
  EXPECT_EQ(Z0, A);
  IAAntisymRMatrix Z1;
  Z1 = mid(A);
  EXPECT_EQ(Z1, IAAntisymRMatrix(mid(A)));
}

TEST_P(IAAntisymRMatrixTest, SummationTest) {
  IAAntisymRMatrix AB = A + B;
  IAAntisymRMatrix W1 = A + mid(B);
  IAAntisymRMatrix W2 = mid(A) + B;
  for (int n = 0; n < dim; ++n)
    for (int m = 0; m < dim; ++m) {
      EXPECT_EQ(AB(n, m), A(n, m) + B(n, m));
      EXPECT_EQ(W1(n, m), A(n, m) + mid(B(n, m)));
      EXPECT_EQ(W2(n, m), mid(A(n, m)) + B(n, m));
    }
}

TEST_P(IAAntisymRMatrixTest, DifferenceTest) {
  IAAntisymRMatrix AB = A - B;
  IAAntisymRMatrix W1 = A - mid(B);
  IAAntisymRMatrix W2 = mid(A) - B;
  for (int n = 0; n < dim; ++n)
    for (int m = 0; m < dim; ++m) {
      EXPECT_EQ(AB(n, m), A(n, m) - B(n, m));
      EXPECT_EQ(W1(n, m), A(n, m) - mid(B(n, m)));
      EXPECT_EQ(W2(n, m), mid(A(n, m)) - B(n, m));
    }
}

TEST_P(IAAntisymRMatrixTest, AdditiveInversTest) {
  IAAntisymRMatrix W = -A;
  for (int n = 0; n < dim; ++n)
    for (int m = 0; m < dim; ++m)
      EXPECT_EQ(W(n, m), -A(n, m));
}

TEST_P(IAAntisymRMatrixTest, ScalarMultiplicationTest) {
  IAAntisymRMatrix W1 = value * A;
  IAAntisymRMatrix W2 = A * value;
  IAAntisymRMatrix W3 = value_d * A;
  IAAntisymRMatrix W4 = A * value_d;
  IAAntisymRMatrix W5 = value_i * A;
  IAAntisymRMatrix W6 = A * value_i;
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

TEST_P(IAAntisymRMatrixTest, ScalarDivisionTest) {
  IAAntisymRMatrix W1 = A / value;
  IAAntisymRMatrix W2 = A / value_d;
  IAAntisymRMatrix W3 = A / value_i;
  for (int n = 0; n < dim; ++n)
    for (int m = 0; m < dim; ++m) {
      EXPECT_EQ(W1(n, m), A(n, m) / value);
      EXPECT_EQ(W2(n, m), A(n, m) / value_d);
      EXPECT_EQ(W3(n, m), A(n, m) / value_i);
    }
}

TEST_P(IAAntisymRMatrixTest, MatrixVectorProductTest) {
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

TEST_P(IAAntisymRMatrixTest, TransposeTest) {
  IAAntisymRMatrix B = transpose(A);
  for (int n = 0; n < dim; ++n)
    for (int m = 0; m < dim; ++m)
      EXPECT_EQ(B(n, m), A(m, n));
}

TEST_P(IAAntisymRMatrixTest, EqualTest) {
  EXPECT_EQ(A, A);
  EXPECT_EQ(B, B);
}

TEST_P(IAAntisymRMatrixTest, NonEqualTest) {
  if (dim == 1) {
    EXPECT_EQ(A, B);
    EXPECT_EQ(B, A);
  } else {
    EXPECT_NE(A, B);
    EXPECT_NE(B, A);
  }
}

TEST_P(IAAntisymRMatrixTest, SaveLoadTest) {
  Saver saver("SaveLoadTest");
  saver << A;
  saver.close();
  Loader loader("SaveLoadTest");
  IAAntisymRMatrix Z;
  loader >> Z;
  loader.close();
  EXPECT_EQ(A, Z);
}

TEST_P(IAAntisymRMatrixTest, InfTest) {
  RMatrix W = inf(A);
  for (int n = 0; n < dim; ++n)
    for (int m = 0; m < dim; ++m)
      EXPECT_EQ(W(n, m), inf(A(n, m)));
}

TEST_P(IAAntisymRMatrixTest, SupTest) {
  RMatrix W = sup(A);
  for (int n = 0; n < dim; ++n)
    for (int m = 0; m < dim; ++m)
      EXPECT_EQ(W(n, m), sup(A(n, m)));
}

TEST_P(IAAntisymRMatrixTest, MidTest) {
  AntisymRMatrix W = mid(A);
  for (int n = 0; n < dim; ++n)
    for (int m = 0; m < dim; ++m)
      EXPECT_EQ(W(n, m), mid(A(n, m)));
}

INSTANTIATE_TEST_CASE_P(IABasicAlgebraMatrixTest, IAAntisymRMatrixTest,
                        Values(IABasicAlgebraMatrixTestParameter{1, 1},
                               IABasicAlgebraMatrixTestParameter{2, 2},
                               IABasicAlgebraMatrixTestParameter{3, 3},
                               IABasicAlgebraMatrixTestParameter{4, 4},
                               IABasicAlgebraMatrixTestParameter{5, 5}));

int main(int argc, char **argv) {
  MppTest mppTest = MppTestBuilder(argc, argv).WithPPM();
  return mppTest.RUN_ALL_MPP_TESTS();
}