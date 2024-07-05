#include "gtest/gtest.h"

#include "TestBasicAlgebra.hpp"
#include "basicalgebra/SymRMatrix.hpp"

class SymRMatrixTest : public BasicAlgebraMatrixTest {
protected:
  SymRMatrix A, B;
  RVector U;
  double value;
  int value_i;

  int dim;
public:
  SymRMatrixTest() : BasicAlgebraMatrixTest() {
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

TEST_P(SymRMatrixTest, ConstructorTest) {
  SymRMatrix W0(dim);
  SymRMatrix W1(value, dim);
  SymRMatrix W2(value_i, dim);
  SymRMatrix W3(A);
  for (int n = 0; n < dim; ++n)
    for (int m = 0; m < dim; ++m) {
      EXPECT_DOUBLE_EQ(W0(n, m), 0.0);
      EXPECT_DOUBLE_EQ(W1(n, m), value);
      EXPECT_DOUBLE_EQ(W2(n, m), value_i);
      EXPECT_DOUBLE_EQ(W3(n, m), A(n, m));
    }
}

TEST_P(SymRMatrixTest, AssignmentTest) {
  SymRMatrix Z = A;
  EXPECT_EQ(Z, A);
  SymRMatrix W(dim);
  W = value;
  EXPECT_EQ(W, SymRMatrix(value, dim));
  W = value_i;
  EXPECT_EQ(W, SymRMatrix(value_i, dim));
}

TEST_P(SymRMatrixTest, SummationTest) {
  SymRMatrix AB = A + B;
  for (int n = 0; n < dim; ++n)
    for (int m = 0; m < dim; ++m)
      EXPECT_DOUBLE_EQ(AB(n, m), A(n, m) + B(n, m));
}

TEST_P(SymRMatrixTest, DifferenceTest) {
  SymRMatrix AB = A - B;
  for (int n = 0; n < dim; ++n)
    for (int m = 0; m < dim; ++m)
      EXPECT_DOUBLE_EQ(AB(n, m), A(n, m) - B(n, m));
}

TEST_P(SymRMatrixTest, AdditiveInversTest) {
  SymRMatrix W = -A;
  for (int n = 0; n < dim; ++n)
    for (int m = 0; m < dim; ++m)
      EXPECT_DOUBLE_EQ(W(n, m), -A(n, m));
}

TEST_P(SymRMatrixTest, ScalarMultiplicationTest) {
  SymRMatrix W1 = value * A;
  SymRMatrix W2 = A * value;
  SymRMatrix W3 = value_i * A;
  SymRMatrix W4 = A * value_i;
  for (int n = 0; n < dim; ++n)
    for (int m = 0; m < dim; ++m) {
      EXPECT_DOUBLE_EQ(W1(n, m), value * A(n, m));
      EXPECT_DOUBLE_EQ(W2(n, m), value * A(n, m));
      EXPECT_DOUBLE_EQ(W3(n, m), value_i * A(n, m));
      EXPECT_DOUBLE_EQ(W4(n, m), value_i * A(n, m));
    }
}

TEST_P(SymRMatrixTest, ScalarDivisionTest) {
  SymRMatrix W1 = A / value;
  SymRMatrix W2 = A / value_i;
  for (int n = 0; n < dim; ++n)
    for (int m = 0; m < dim; ++m) {
      EXPECT_DOUBLE_EQ(W1(n, m), A(n, m) / value);
      EXPECT_DOUBLE_EQ(W2(n, m), A(n, m) / value_i);
    }
}

TEST_P(SymRMatrixTest, MatrixVectorProductTest) {
  RVector W0(dim);
  for (int n = 0; n < dim; ++n)
    for (int k = 0; k < dim; ++k)
      W0[n] += A(n, k) * U[k];
  EXPECT_EQ(W0, A * U);
}

TEST_P(SymRMatrixTest, EqualTest) {
  EXPECT_EQ(A, A);
  EXPECT_EQ(B, B);
}

TEST_P(SymRMatrixTest, NonEqualTest) {
  EXPECT_NE(A, B);
  EXPECT_NE(B, A);
}

TEST_P(SymRMatrixTest, SaveLoadTest) {
  Saver saver("SaveLoadTest");
  saver << A;
  saver.close();
  Loader loader("SaveLoadTest");
  SymRMatrix Z;
  loader >> Z;
  loader.close();
  EXPECT_EQ(A, Z);
}

INSTANTIATE_TEST_CASE_P(BasicAlgebraMatrixTest, SymRMatrixTest,
                        Values(BasicAlgebraMatrixTestParameter{1, 1},
                               BasicAlgebraMatrixTestParameter{2, 2},
                               BasicAlgebraMatrixTestParameter{3, 3},
                               BasicAlgebraMatrixTestParameter{4, 4},
                               BasicAlgebraMatrixTestParameter{5, 5}));

int main(int argc, char **argv) {
  MppTest mppTest = MppTestBuilder(argc, argv).WithPPM();
  return mppTest.RUN_ALL_MPP_TESTS();
}
