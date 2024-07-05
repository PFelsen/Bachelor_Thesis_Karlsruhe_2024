#include "gtest/gtest.h"

#include "TestBasicAlgebra.hpp"
#include "basicalgebra/AntisymRMatrix.hpp"

class AntisymRMatrixTest : public BasicAlgebraMatrixTest {
protected:
  AntisymRMatrix A, B;
  RVector U;
  double value;
  int value_i;

  int dim;
public:
  AntisymRMatrixTest() : BasicAlgebraMatrixTest() {
    dim = cols = rows;
    A.resize(dim);
    B.resize(dim);
    U.resize(dim);

    for (int n = 0; n < dim; ++n) {
      U[n] = RandomDouble();
      for (int m = n; m < dim; ++m) {
        A(RandomDouble(), n, m);
        B(RandomDouble(), n, m);
      }
    }
    value = RandomNonZeroDouble();
    value_i = RandomNonZeroInt();
  }
};

TEST_P(AntisymRMatrixTest, ConstructorTest) {
  AntisymRMatrix W0(dim);
  AntisymRMatrix W1(A);
  for (int n = 0; n < dim; ++n)
    for (int m = 0; m < dim; ++m) {
      EXPECT_DOUBLE_EQ(W0(n, m), 0.0);
      EXPECT_DOUBLE_EQ(W1(n, m), A(n, m));
    }
}

TEST_P(AntisymRMatrixTest, AssignmentTest) {
  AntisymRMatrix Z = A;
  EXPECT_EQ(Z, A);
}

TEST_P(AntisymRMatrixTest, SummationTest) {
  AntisymRMatrix AB = A + B;
  for (int n = 0; n < dim; ++n)
    for (int m = 0; m < dim; ++m)
      EXPECT_DOUBLE_EQ(AB(n, m), A(n, m) + B(n, m));
}

TEST_P(AntisymRMatrixTest, DifferenceTest) {
  AntisymRMatrix AB = A - B;
  for (int n = 0; n < dim; ++n)
    for (int m = 0; m < dim; ++m)
      EXPECT_DOUBLE_EQ(AB(n, m), A(n, m) - B(n, m));
}

TEST_P(AntisymRMatrixTest, AdditiveInversTest) {
  AntisymRMatrix W = -A;
  for (int n = 0; n < dim; ++n)
    for (int m = 0; m < dim; ++m)
      EXPECT_DOUBLE_EQ(W(n, m), -A(n, m));
}

TEST_P(AntisymRMatrixTest, ScalarMultiplicationTest) {
  AntisymRMatrix W1 = value * A;
  AntisymRMatrix W2 = A * value;
  AntisymRMatrix W3 = value_i * A;
  AntisymRMatrix W4 = A * value_i;
  for (int n = 0; n < dim; ++n)
    for (int m = 0; m < dim; ++m) {
      EXPECT_DOUBLE_EQ(W1(n, m), value * A(n, m));
      EXPECT_DOUBLE_EQ(W2(n, m), value * A(n, m));
      EXPECT_DOUBLE_EQ(W3(n, m), value_i * A(n, m));
      EXPECT_DOUBLE_EQ(W4(n, m), value_i * A(n, m));
    }
}

TEST_P(AntisymRMatrixTest, ScalarDivisionTest) {
  AntisymRMatrix W1 = A / value;
  AntisymRMatrix W2 = A / value_i;
  for (int n = 0; n < dim; ++n)
    for (int m = 0; m < dim; ++m) {
      EXPECT_DOUBLE_EQ(W1(n, m), A(n, m) / value);
      EXPECT_DOUBLE_EQ(W2(n, m), A(n, m) / value_i);
    }
}

TEST_P(AntisymRMatrixTest, MatrixVectorProductTest) {
  RVector W0(dim);
  for (int n = 0; n < dim; ++n)
    for (int k = 0; k < dim; ++k)
      W0[n] += A(n, k) * U[k];
  EXPECT_EQ(W0, A * U);
}

TEST_P(AntisymRMatrixTest, TransposeTest) {
  AntisymRMatrix B = transpose(A);
  for (int n = 0; n < dim; ++n)
    for (int m = 0; m < dim; ++m)
      EXPECT_EQ(B(n, m), A(m, n));
}

TEST_P(AntisymRMatrixTest, EqualTest) {
  EXPECT_EQ(A, A);
  EXPECT_EQ(B, B);
}

TEST_P(AntisymRMatrixTest, NonEqualTest) {
  if (dim == 1) {
    EXPECT_EQ(A, B);
    EXPECT_EQ(B, A);
  } else {
    EXPECT_NE(A, B);
    EXPECT_NE(B, A);
  }
}

TEST_P(AntisymRMatrixTest, SaveLoadTest) {
  Saver saver("SaveLoadTest");
  saver << A;
  saver.close();
  Loader loader("SaveLoadTest");
  AntisymRMatrix Z;
  loader >> Z;
  loader.close();
  EXPECT_EQ(A, Z);
}

INSTANTIATE_TEST_CASE_P(BasicAlgebraMatrixTest, AntisymRMatrixTest,
                        Values(BasicAlgebraMatrixTestParameter{1, 1},
                               BasicAlgebraMatrixTestParameter{2, 2},
                               BasicAlgebraMatrixTestParameter{3, 3},
                               BasicAlgebraMatrixTestParameter{4, 4},
                               BasicAlgebraMatrixTestParameter{5, 5}));

int main(int argc, char **argv) {
  MppTest mppTest = MppTestBuilder(argc, argv).WithPPM();
  return mppTest.RUN_ALL_MPP_TESTS();
}
