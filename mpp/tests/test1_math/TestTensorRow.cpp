#include "Tensor.hpp"
#include "TestEnvironment.hpp"

#include <vector>

template<int dim>
class TensorRowTest : public ::testing::Test {
protected:
  using VECTORFIELD = VectorFieldT<double, dim>;
  using TENSOR = TensorT<double, dim>;
  using COMPONENT = VectorFieldComponentT<double>;
  using ROW = TensorRowT<double, dim>;

  double scalar;
  int scalar_i;

  VECTORFIELD U, V;
  TENSOR T;

  TensorRowTest() {
    mpp_geometry::SetTolerance(1e-15);
    scalar = RandomNonZeroDouble();
    scalar_i = RandomNonZeroInt();
    for (int i = 0; i < dim; ++i) {
      U[i] = RandomNonZeroDouble();
      V[i] = RandomNonZeroDouble();
      for (int j = 0; j < dim; ++j)
        T[i][j] = RandomNonZeroDouble();
    }
  }

  void checkConstructor() {
    ROW z;
    EXPECT_EQ(z.Value(), VECTORFIELD{});
    EXPECT_EQ(z.Component(), 0);

    for (int k = 0; k < dim; ++k) {
      ROW v(k, V);
      EXPECT_EQ(v.Value(), V);
      EXPECT_EQ(v.Component(), k);
      for (int i = 0; i < dim; ++i) {
        if (i == k) EXPECT_EQ(V, v[i]);
        else EXPECT_EQ(VECTORFIELD{}, v[i]);
      }

      ROW u(v);
      EXPECT_EQ(u.Value(), V);
      EXPECT_EQ(u.Component(), k);

      TENSOR w(v);
      for (int i = 0; i < dim; ++i) {
        if (i == k) EXPECT_EQ(V, w[i]);
        else EXPECT_EQ(VECTORFIELD{}, w[i]);
      }
    }
  }

  virtual void checkAssignment() {
    for (int k = 0; k < dim; ++k) {
      ROW v(k, V);
      ROW u = v;
      EXPECT_EQ(u.Value(), V);
      EXPECT_EQ(u.Component(), k);
    }
  }

  virtual void checkSummation() {
    for (int k = 0; k < dim; ++k) {
      ROW a(k, V);
      TENSOR A(a);
      EXPECT_EQ(a + T, A + T);
      EXPECT_EQ(T + a, T + A);
    }
  }

  virtual void checkDifference() {
    for (int k = 0; k < dim; ++k) {
      ROW a(k, V);
      TENSOR A(a);
      EXPECT_EQ(a - T, A - T);
      EXPECT_EQ(T - a, T - A);
    }
  }

  virtual void checkAdditiveInverse() {
    for (int k = 0; k < dim; ++k) {
      ROW a(k, V);
      TENSOR A(a);
      EXPECT_EQ(TENSOR(-a), -A);
    }
  }

  virtual void checkScalarMultiplication() {
    for (int k = 0; k < dim; ++k) {
      ROW v(k, V);

      auto r1 = v * scalar;
      auto r2 = scalar * v;
      auto r3 = v * scalar_i;
      auto r4 = scalar_i * v;

      EXPECT_EQ(r1.Value(), V * scalar);
      EXPECT_EQ(r2.Value(), V * scalar);
      EXPECT_EQ(r3.Value(), V * scalar_i);
      EXPECT_EQ(r4.Value(), V * scalar_i);
      EXPECT_EQ(r1.Component(), k);
      EXPECT_EQ(r2.Component(), k);
      EXPECT_EQ(r3.Component(), k);
      EXPECT_EQ(r4.Component(), k);
    }
  }

  virtual void checkScalarDivision() {
    for (int k = 0; k < dim; ++k) {
      ROW v(k, V);

      auto r1 = v / scalar;
      auto r3 = v / scalar_i;

      EXPECT_EQ(r1.Value(), V / scalar);
      EXPECT_EQ(r3.Value(), V / scalar_i);
      EXPECT_EQ(r1.Component(), k);
      EXPECT_EQ(r3.Component(), k);
    }
  }

  virtual void checkVectorTensorProduct() {
    for (int k = 0; k < dim; ++k) {
      ROW A(k, V);
      TENSOR B(A);

      EXPECT_EQ(VECTORFIELD(A * U), B * U);
      EXPECT_EQ(VECTORFIELD(A.multiply(U)), B * U);
      EXPECT_EQ(U * A, U * B);
      EXPECT_EQ(A.applyTransposed(U), U * B);

      for (int j = 0; j < dim; ++j) {
        COMPONENT v(j, scalar);

        EXPECT_EQ(VECTORFIELD(A * v), B * v);
        EXPECT_EQ(VECTORFIELD(A.multiply(v)), B * v);
        EXPECT_EQ(v * A, v * B);
        EXPECT_EQ(A.applyTransposed(v), v * B);
      }
    }
  }

  virtual void checkTensorTensorProduct() {
    for (int k = 0; k < dim; ++k) {
      ROW a(k, V);
      TENSOR A(a);
      for (int l = 0; l < dim; ++l) {
        ROW b(l, U);
        TENSOR B(b);

        EXPECT_EQ(TENSOR(a * b), A * B);
      }
    }
  }

  virtual void checkColRowProduct() {
    for (int k = 0; k < dim; ++k) {
      ROW a(k, V);
      TENSOR A(a);
      for (int l = 0; l < dim; ++l) {
        ROW b(l, U);
        TENSOR B(b);

        EXPECT_EQ(ColRowProduct(a, b), transpose(A) * B);
      }
    }
  }

  virtual void checkFrobenius() {
    for (int k = 0; k < dim; ++k) {
      ROW v(k, V);
      TENSOR a(v);
      for (int l = 0; l < dim; ++l) {
        ROW u(l, U);
        TENSOR b(u);
        EXPECT_DOUBLE_EQ(Frobenius(v, u), Frobenius(a, b));
      }
    }
  }

  virtual void checkNorm() {
    for (int k = 0; k < dim; ++k) {
      ROW v(k, V);
      TENSOR a(v);
      EXPECT_DOUBLE_EQ(norm(v), norm(a));
    }
  }

  virtual void checkTrace() {
    for (int k = 0; k < dim; ++k) {
      ROW v(k, V);
      TENSOR a(v);
      EXPECT_DOUBLE_EQ(trace(v), trace(a));
    }
  }

  virtual void checkEQ() {
    for (int k = 0; k < dim; ++k) {
      ROW v(k, V);
      EXPECT_EQ(v, v);
    }
  }

  virtual void checkNonEQ() {
    std::vector<std::vector<ROW>> A;
    for (int k = 0; k < dim; ++k) {
      std::vector<ROW> A_k{};
      for (int l = 0; l < dim; ++l) {
        A_k.push_back(ROW(k, V));
      }
      A.push_back(A_k);
    }
    for (int k = 0; k < dim; ++k) {
      for (int l = 0; l < dim; ++l) {
        if (k != l) EXPECT_NE(A[k], A[l]);
      }
    }
  }

  virtual void checkSaveLoad() {
    for (int k = 0; k < dim; ++k) {
      ROW v(k, V);
      Saver saver("SaveLoadTest");
      saver << v;
      saver.close();
      Loader loader("SaveLoadTest");
      ROW w;
      loader >> w;
      loader.close();
      EXPECT_EQ(v, w);
    }
  }
};

/// To avoid the same code in each test
#define TENSOR_ROW_TESTS(tensorClass)                                                              \
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
  TEST_F(tensorClass, VectorTensorProductTest) { checkVectorTensorProduct(); }                     \
                                                                                                   \
  TEST_F(tensorClass, TensorTensorProductTest) { checkTensorTensorProduct(); }                     \
                                                                                                   \
  TEST_F(tensorClass, ColRowProductTest) { checkColRowProduct(); }                                 \
                                                                                                   \
  TEST_F(tensorClass, FrobeniusTest) { checkFrobenius(); }                                         \
                                                                                                   \
  TEST_F(tensorClass, NormTest) { checkNorm(); }                                                   \
                                                                                                   \
  TEST_F(tensorClass, TraceTest) { checkTrace(); }                                                 \
                                                                                                   \
  TEST_F(tensorClass, EqualTest) { checkEQ(); }                                                    \
                                                                                                   \
  TEST_F(tensorClass, NonEqualTest) { checkNonEQ(); }                                              \
                                                                                                   \
  TEST_F(tensorClass, SaveLoadTest) { checkSaveLoad(); }

using TensorRowTest_3 = TensorRowTest<SpaceDimension>;

TENSOR_ROW_TESTS(TensorRowTest_3)

/*using TensorRowTest_2 = TensorRowTest<2>;

TENSOR_ROW_TESTS(TensorRowTest_2)

using TensorRowTest_1 = TensorRowTest<1>;

TENSOR_ROW_TESTS(TensorRowTest_1)*/

int main(int argc, char **argv) {
  MppTest mppTest = MppTestBuilder(argc, argv).WithPPM();
  return mppTest.RUN_ALL_MPP_TESTS();
}
