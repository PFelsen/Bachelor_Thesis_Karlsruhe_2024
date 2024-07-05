#include "SymTensor.hpp"
#include "TestEnvironment.hpp"

#include <vector>

template<int sDim>
class IASymTensorTest : public ::testing::Test {
protected:
  using IASYMTENSOR = SymTensorT<IAInterval, sDim>;
  using IAVECTORFIELD = VectorFieldT<IAInterval, sDim>;

  int dim = -1;
  IASYMTENSOR A, B;
  IAVECTORFIELD v;
  IAInterval value;
  double value_d;
  int value_i;

  std::vector<std::vector<IAInterval>> values;
  std::vector<IASYMTENSOR> tensors;

  IASymTensorTest() {
    mpp_geometry::SetTolerance(1e-15);
    dim = A.Dim();
    for (int i = 0; i < A.size(); ++i) {
      A[i] = RandomIAInterval();
      B[i] = RandomIAInterval();
    }
    for (int i = 0; i < dim; ++i)
      v[i] = RandomIAInterval();
    value = RandomNonZeroIAInterval();
    value_d = RandomNonZeroDouble();
    value_i = RandomNonZeroInt();


    values.resize(dim);
    tensors.resize(dim);
    for (int i = 0; i < dim; ++i)
      values[i].resize(dim);
    for (int i = 0; i < dim; ++i) {
      values[i][i] = RandomIAInterval();
      for (int j = 0; j < i; ++j)
        values[i][j] = values[j][i] = RandomIAInterval();
    }

    if constexpr (sDim >= 2) tensors[1] = IASYMTENSOR(values[0][0], values[1][0], values[1][1]);
    if constexpr (sDim >= 3)
      tensors[2] = IASYMTENSOR(values[0][0], values[1][0], values[1][1], values[2][0], values[2][1],
                               values[2][2]);
  }

  void checkConstructor() {
    IASYMTENSOR z0;
    IASYMTENSOR z1(A);
    IASYMTENSOR z2(mid(A));
    IASYMTENSOR z3(value);
    IASYMTENSOR z4(value_d);
    IASYMTENSOR z5(value_i);
    for (int i = 0; i < dim; ++i)
      for (int j = 0; j < dim; ++j) {
        EXPECT_EQ(z0(i, j), IAInterval());
        EXPECT_EQ(z1(i, j), A(i, j));
        EXPECT_EQ(z2(i, j), IAInterval(mid(A(i, j))));
        EXPECT_EQ(z3(i, j), value);
        EXPECT_EQ(z4(i, j), IAInterval(value_d));
        EXPECT_EQ(z5(i, j), IAInterval(value_i));
      }

    for (int i = 1; i < dim; ++i) {
      for (int k = 0; k < dim; ++k) {
        for (int l = 0; l < dim; ++l) {
          if (k <= i && l <= i) EXPECT_EQ(tensors[i](k, l), values[k][l]);
          else if (k == l) EXPECT_EQ(tensors[i](k, l), IAInterval(1.0));
          else EXPECT_EQ(tensors[i](k, l), IAInterval());
        }
      }
    }
  }

  virtual void checkAssignment() {
    IASYMTENSOR z1 = A;
    IASYMTENSOR z2 = mid(A);
    IASYMTENSOR z3(A);
    z3 = value;
    IASYMTENSOR z4(A);
    z4 = value_d;
    IASYMTENSOR z5(A);
    z5 = value_i;
    for (int i = 0; i < dim; ++i)
      for (int j = 0; j < dim; ++j) {
        EXPECT_EQ(z1(i, j), A(i, j));
        EXPECT_EQ(z2(i, j), IAInterval(mid(A(i, j))));
        EXPECT_EQ(z3(i, j), value);
        EXPECT_EQ(z4(i, j), IAInterval(value_d));
        EXPECT_EQ(z5(i, j), IAInterval(value_i));
      }
  }

  virtual void checkSummation() {
    auto r0 = A + B;
    auto r1 = mid(A) + B;
    auto r2 = A + mid(B);
    for (int i = 0; i < dim; ++i)
      for (int j = 0; j < dim; ++j) {
        EXPECT_EQ(r0(i, j), A(i, j) + B(i, j));
        EXPECT_EQ(r1(i, j), mid(A(i, j)) + B(i, j));
        EXPECT_EQ(r2(i, j), A(i, j) + mid(B(i, j)));
      }
  }

  virtual void checkDifference() {
    auto r0 = A - B;
    auto r1 = mid(A) - B;
    auto r2 = A - mid(B);
    for (int i = 0; i < dim; ++i)
      for (int j = 0; j < dim; ++j) {
        EXPECT_EQ(r0(i, j), A(i, j) - B(i, j));
        EXPECT_EQ(r1(i, j), mid(A(i, j)) - B(i, j));
        EXPECT_EQ(r2(i, j), A(i, j) - mid(B(i, j)));
      }
  }

  virtual void checkAdditiveInverse() {
    auto r = -A;
    for (int i = 0; i < dim; ++i)
      for (int j = 0; j < dim; ++j)
        EXPECT_EQ(r(i, j), -A(i, j));
  }

  virtual void checkScalarMultiplication() {
    auto r1 = A * value;
    auto r2 = value * A;
    auto r3 = A * value_d;
    auto r4 = value_d * A;
    auto r5 = A * value_i;
    auto r6 = value_i * A;
    auto r7 = value * mid(A);
    auto r8 = mid(A) * value;
    for (int i = 0; i < dim; ++i)
      for (int j = 0; j < dim; ++j) {
        EXPECT_EQ(r1(i, j), A(i, j) * value);
        EXPECT_EQ(r2(i, j), A(i, j) * value);
        EXPECT_EQ(r3(i, j), A(i, j) * value_d);
        EXPECT_EQ(r4(i, j), A(i, j) * value_d);
        EXPECT_EQ(r5(i, j), A(i, j) * value_i);
        EXPECT_EQ(r6(i, j), A(i, j) * value_i);
        EXPECT_EQ(r7(i, j), mid(A(i, j)) * value);
        EXPECT_EQ(r8(i, j), mid(A(i, j)) * value);
      }
  }

  virtual void checkScalarDivision() {
    auto r1 = A / value;
    auto r2 = A / value_d;
    auto r3 = A / value_i;
    auto r4 = mid(A) / value;
    for (int i = 0; i < dim; ++i)
      for (int j = 0; j < dim; ++j) {
        EXPECT_EQ(r1(i, j), A(i, j) / value);
        EXPECT_EQ(r2(i, j), A(i, j) / value_d);
        EXPECT_EQ(r3(i, j), A(i, j) / value_i);
        EXPECT_EQ(r4(i, j), mid(A(i, j)) / value);
      }
  }

  virtual void checkMatrixVectorProduct() {
    IAVECTORFIELD r0, r1, r2;
    for (int i = 0; i < dim; ++i)
      for (int k = 0; k < dim; ++k) {
        r0[i] += A(i, k) * v[k];
        r1[i] += A(i, k) * mid(v[k]);
        r2[i] += mid(A(i, k)) * v[k];
      }
    EXPECT_EQ(r0, A * v);
    EXPECT_EQ(r1, A * mid(v));
    EXPECT_EQ(r2, mid(A) * v);
  }

  virtual void checkDet() {
    if (A.Dim() == 1) {
      IAInterval d = A(0, 0);
      EXPECT_EQ(d, det(A));
    } else if (A.Dim() == 2) {
      IAInterval d = A(0, 0) * A(1, 1) - A(1, 0) * A(0, 1);
      EXPECT_EQ(d, det(A));
    } else if (A.Dim() == 3) {
      IAInterval d = A(0, 0) * (A(1, 1) * A(2, 2) - A(2, 1) * A(1, 2))
                     + A(1, 0) * (A(2, 1) * A(0, 2) - A(0, 1) * A(2, 2))
                     + A(2, 0) * (A(0, 1) * A(1, 2) - A(1, 1) * A(0, 2));
      EXPECT_EQ(d, det(A));
    } else {
      EXPECT_ANY_THROW(det(A));
    }
  }

  virtual void checkTrace() {
    IAInterval t;
    for (int i = 0; i < dim; ++i)
      t += A(i, i);
    EXPECT_EQ(t, trace(A));
  }

  virtual void checkNorm() {
    IAInterval sum;
    for (int i = 0; i < dim; ++i)
      for (int j = 0; j < dim; ++j)
        sum += sqr(A(i, j));
    EXPECT_EQ(norm(A), sqrt(sum));
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
    IASYMTENSOR C;
    loader >> C;
    loader.close();
    EXPECT_EQ(A, C);
  }

  virtual void checkMid() {
    auto r = mid(A);
    for (int i = 0; i < dim; ++i)
      for (int j = 0; j < dim; ++j)
        EXPECT_DOUBLE_EQ(r(i, j), mid(A(i, j)));
  }

  virtual void checkInf() {
    auto r = inf(A);
    for (int i = 0; i < dim; ++i)
      for (int j = 0; j < dim; ++j)
        EXPECT_DOUBLE_EQ(r(i, j), inf(A(i, j)));
  }

  virtual void checkSup() {
    auto r = sup(A);
    for (int i = 0; i < dim; ++i)
      for (int j = 0; j < dim; ++j)
        EXPECT_DOUBLE_EQ(r(i, j), sup(A(i, j)));
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
  TEST_F(tensorClass, SaveLoadTest) { checkSaveLoad(); }                                           \
                                                                                                   \
  TEST_F(tensorClass, MidTest) { checkMid(); }                                                     \
                                                                                                   \
  TEST_F(tensorClass, InfTest) { checkInf(); }                                                     \
                                                                                                   \
  TEST_F(tensorClass, SupTest) { checkSup(); }

using IASymTensorTest_3 = IASymTensorTest<3>;

SYMTENSOR_TESTS(IASymTensorTest_3)

using IASymTensorTest_2 = IASymTensorTest<2>;

SYMTENSOR_TESTS(IASymTensorTest_2)

using IASymTensorTest_1 = IASymTensorTest<1>;

SYMTENSOR_TESTS(IASymTensorTest_1)

int main(int argc, char **argv) {
  MppTest mppTest = MppTestBuilder(argc, argv).WithPPM();
  return mppTest.RUN_ALL_MPP_TESTS();
}
