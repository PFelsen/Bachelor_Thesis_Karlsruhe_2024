#include "SymTensor.hpp"
#include "TestEnvironment.hpp"

#include <vector>

template<int sDim>
class SymTensorTest : public ::testing::Test {
protected:
  using SYMTENSOR = SymTensorT<double, sDim>;
  using VECTORFIELD = VectorFieldT<double, sDim>;

  int dim = -1;
  SYMTENSOR A, B;
  VECTORFIELD v;
  double value;
  int value_i;

  std::vector<std::vector<double>> values;
  std::vector<SYMTENSOR> tensors;

  SymTensorTest() {
    mpp_geometry::SetTolerance(1e-15);
    dim = A.Dim();
    for (int i = 0; i < A.size(); ++i) {
      A[i] = RandomDouble();
      B[i] = RandomDouble();
    }
    for (int i = 0; i < dim; ++i)
      v[i] = RandomDouble();
    value = RandomNonZeroDouble();
    value_i = RandomNonZeroInt();


    values.resize(dim);
    tensors.resize(dim);
    for (int i = 0; i < dim; ++i)
      values[i].resize(dim);
    for (int i = 0; i < dim; ++i) {
      values[i][i] = RandomDouble();
      for (int j = 0; j < i; ++j)
        values[i][j] = values[j][i] = RandomDouble();
    }

    if constexpr (sDim >= 2) tensors[1] = SYMTENSOR(values[0][0], values[1][0], values[1][1]);
    if constexpr (sDim >= 3)
      tensors[2] = SYMTENSOR(values[0][0], values[1][0], values[1][1], values[2][0], values[2][1],
                             values[2][2]);
  }

  void checkConstructor() {
    SYMTENSOR z0;
    SYMTENSOR z1(A);
    SYMTENSOR z2(value);
    SYMTENSOR z3(value_i);
    for (int i = 0; i < dim; ++i)
      for (int j = 0; j < dim; ++j) {
        EXPECT_DOUBLE_EQ(z0(i, j), 0.0);
        EXPECT_DOUBLE_EQ(z1(i, j), A(i, j));
        EXPECT_DOUBLE_EQ(z2(i, j), value);
        EXPECT_DOUBLE_EQ(z3(i, j), double(value_i));
      }

    for (int i = 1; i < dim; ++i) {
      for (int k = 0; k < dim; ++k) {
        for (int l = 0; l < dim; ++l) {
          if (k <= i && l <= i) EXPECT_DOUBLE_EQ(tensors[i](k, l), values[k][l]);
          else if (k == l) EXPECT_DOUBLE_EQ(tensors[i](k, l), 1.0);
          else EXPECT_DOUBLE_EQ(tensors[i](k, l), 0.0);
        }
      }
    }
  }

  virtual void checkAssignment() {
    SYMTENSOR z1 = A;
    SYMTENSOR z2(A);
    z2 = value;
    SYMTENSOR z3(A);
    z3 = value_i;
    for (int i = 0; i < dim; ++i)
      for (int j = 0; j < dim; ++j) {
        EXPECT_DOUBLE_EQ(z1(i, j), A(i, j));
        EXPECT_DOUBLE_EQ(z2(i, j), value);
        EXPECT_DOUBLE_EQ(z3(i, j), double(value_i));
      }
  }

  virtual void checkSummation() {
    auto r0 = A + B;
    for (int i = 0; i < dim; ++i)
      for (int j = 0; j < dim; ++j) {
        EXPECT_DOUBLE_EQ(r0(i, j), A(i, j) + B(i, j));
      }
  }

  virtual void checkDifference() {
    auto r0 = A - B;
    for (int i = 0; i < dim; ++i)
      for (int j = 0; j < dim; ++j) {
        EXPECT_DOUBLE_EQ(r0(i, j), A(i, j) - B(i, j));
      }
  }

  virtual void checkAdditiveInverse() {
    auto r = -A;
    for (int i = 0; i < dim; ++i)
      for (int j = 0; j < dim; ++j)
        EXPECT_DOUBLE_EQ(r(i, j), -A(i, j));
  }

  virtual void checkScalarMultiplication() {
    auto r1 = A * value;
    auto r2 = value * A;
    auto r3 = A * value_i;
    auto r4 = value_i * A;
    for (int i = 0; i < dim; ++i)
      for (int j = 0; j < dim; ++j) {
        EXPECT_DOUBLE_EQ(r1(i, j), A(i, j) * value);
        EXPECT_DOUBLE_EQ(r2(i, j), A(i, j) * value);
        EXPECT_DOUBLE_EQ(r3(i, j), A(i, j) * value_i);
        EXPECT_DOUBLE_EQ(r4(i, j), A(i, j) * value_i);
      }
  }

  virtual void checkScalarDivision() {
    auto r1 = A / value;
    auto r2 = A / value_i;
    for (int i = 0; i < dim; ++i)
      for (int j = 0; j < dim; ++j) {
        EXPECT_DOUBLE_EQ(r1(i, j), A(i, j) / value);
        EXPECT_DOUBLE_EQ(r2(i, j), A(i, j) / value_i);
      }
  }

  virtual void checkMatrixVectorProduct() {
    VECTORFIELD r0;
    for (int i = 0; i < dim; ++i)
      for (int k = 0; k < dim; ++k) {
        r0[i] += A(i, k) * v[k];
      }
    EXPECT_EQ(r0, A * v);
  }

  virtual void checkDet() {
    if (A.Dim() == 1) {
      double d = A(0, 0);
      EXPECT_DOUBLE_EQ(d, det(A));
    } else if (A.Dim() == 2) {
      double d = A(0, 0) * A(1, 1) - A(1, 0) * A(0, 1);
      EXPECT_DOUBLE_EQ(d, det(A));
    } else if (A.Dim() == 3) {
      double d = A(0, 0) * (A(1, 1) * A(2, 2) - A(2, 1) * A(1, 2))
                 + A(1, 0) * (A(2, 1) * A(0, 2) - A(0, 1) * A(2, 2))
                 + A(2, 0) * (A(0, 1) * A(1, 2) - A(1, 1) * A(0, 2));
      EXPECT_DOUBLE_EQ(d, det(A));
    } else {
      EXPECT_ANY_THROW(det(A));
    }
  }

  virtual void checkTrace() {
    double t = 0.0;
    for (int i = 0; i < dim; ++i)
      t += A(i, i);
    EXPECT_DOUBLE_EQ(t, trace(A));
  }

  virtual void checkNorm() {
    double sum = 0.0;
    for (int i = 0; i < dim; ++i)
      for (int j = 0; j < dim; ++j)
        sum += A(i, j) * A(i, j);
    EXPECT_DOUBLE_EQ(norm(A), sqrt(sum));
  }

  virtual void checkEQ() {
    EXPECT_EQ(A, A);
    EXPECT_EQ(B, B);
  }

  virtual void checkNonEQ() {
    EXPECT_NE(A, B);
    EXPECT_NE(B, A);
  }

  virtual void checkSaveLoad() {
    Saver saver("SaveLoadTest");
    saver << A;
    saver.close();
    Loader loader("SaveLoadTest");
    SYMTENSOR C;
    loader >> C;
    loader.close();
    EXPECT_EQ(A, C);
  }
};

/// To avoid the same code in each test
#define SYMTENSOR_TESTS(tensorClass)                                                               \
                                                                                                   \
  TEST_F(tensorClass, ConstructorTest) { checkConstructor(); }                                     \
                                                                                                   \
  TEST_F(tensorClass, AssignmentTest) { checkAssignment(); }                                       \
                                                                                                   \
  TEST_F(tensorClass, SummationTest) { checkSummation(); }                                         \
                                                                                                   \
  TEST_F(tensorClass, DifferenceTest) { checkDifference(); }                                       \
                                                                                                   \
  TEST_F(tensorClass, AdditiveInverseTest) { checkAdditiveInverse(); }                             \
                                                                                                   \
  TEST_F(tensorClass, ScalarMultiplicationTest) { checkScalarMultiplication(); }                   \
                                                                                                   \
  TEST_F(tensorClass, ScalarDivisionTest) { checkScalarDivision(); }                               \
                                                                                                   \
  TEST_F(tensorClass, MatrixVectorProductTest) { checkMatrixVectorProduct(); }                     \
                                                                                                   \
  TEST_F(tensorClass, DeterminantTest) { checkDet(); }                                             \
                                                                                                   \
  TEST_F(tensorClass, TraceTest) { checkTrace(); }                                                 \
                                                                                                   \
  TEST_F(tensorClass, NormTest) { checkNorm(); }                                                   \
                                                                                                   \
  TEST_F(tensorClass, EqualTest) { checkEQ(); }                                                    \
                                                                                                   \
  TEST_F(tensorClass, NonEqualTest) { checkNonEQ(); }                                              \
                                                                                                   \
  TEST_F(tensorClass, SaveLoadTest) { checkSaveLoad(); }

using SymTensorTest_3 = SymTensorTest<3>;

SYMTENSOR_TESTS(SymTensorTest_3)

using SymTensorTest_2 = SymTensorTest<2>;

SYMTENSOR_TESTS(SymTensorTest_2)

using SymTensorTest_1 = SymTensorTest<1>;

SYMTENSOR_TESTS(SymTensorTest_1)

int main(int argc, char **argv) {
  MppTest mppTest = MppTestBuilder(argc, argv).WithPPM();
  return mppTest.RUN_ALL_MPP_TESTS();
}
