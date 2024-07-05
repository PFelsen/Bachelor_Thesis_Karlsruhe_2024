#include "Tensor.hpp"
#include "TestEnvironment.hpp"

#include <vector>

template<int dim>
class TensorComponentTest : public ::testing::Test {
protected:
  using VECTORFIELD = VectorFieldT<double, dim>;
  using TENSOR = TensorT<double, dim>;
  using VCOMPONENT = VectorFieldComponentT<double>;
  using TCOMPONENT = TensorComponentT<double>;
  using ROW = TensorRowT<double, dim>;

  double scalar;
  int scalar_i;

  VECTORFIELD U, V;
  TENSOR T;

  TensorComponentTest() {
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
    TCOMPONENT z;
    EXPECT_EQ(z.Value(), double{});
    EXPECT_EQ(z.Row(), 0);
    EXPECT_EQ(z.Column(), 0);

    for (int k = 0; k < dim; ++k) {
      for (int l = 0; l < dim; ++l) {
        TCOMPONENT v(k, l, scalar);
        EXPECT_EQ(v.Value(), scalar);
        EXPECT_EQ(v.Row(), k);
        EXPECT_EQ(v.Column(), l);
        for (int i = 0; i < dim; ++i) {
          if (i == k) {
            for (int j = 0; j < dim; ++j) {
              if (j == l) {
                EXPECT_EQ(v.Entry(i, j), scalar);
              } else {
                EXPECT_EQ(v.Entry(i, j), double{});
              }
            }
          }
        }


        TCOMPONENT u(v);
        EXPECT_EQ(u.Value(), scalar);
        EXPECT_EQ(u.Row(), k);
        EXPECT_EQ(u.Column(), l);

        TENSOR w(v);
        for (int i = 0; i < dim; ++i) {
          if (i == k) {
            for (int j = 0; j < dim; ++j) {
              if (j == l) {
                EXPECT_EQ(w[i][j], scalar);
              } else {
                EXPECT_EQ(w[i][j], double{});
              }
            }
          }
        }
      }
    }
  }

  virtual void checkAssignment() {
    for (int k = 0; k < dim; ++k) {
      for (int l = 0; l < dim; ++l) {
        TCOMPONENT v(k, l, scalar);
        TCOMPONENT u = v;
        EXPECT_EQ(u.Value(), scalar);
        EXPECT_EQ(u.Row(), k);
        EXPECT_EQ(u.Column(), l);
      }
    }
  }

  virtual void checkSummation() {
    for (int k = 0; k < dim; ++k) {
      for (int l = 0; l < dim; ++l) {
        TCOMPONENT a(k, l, scalar);
        TENSOR A(a);
        EXPECT_EQ(a + T, A + T);
        EXPECT_EQ(T + a, T + A);
      }
    }
  }

  virtual void checkDifference() {
    for (int k = 0; k < dim; ++k) {
      for (int l = 0; l < dim; ++l) {
        TCOMPONENT a(k, l, scalar);
        TENSOR A(a);
        EXPECT_EQ(a - T, A - T);
        EXPECT_EQ(T - a, T - A);
      }
    }
  }

  virtual void checkAdditiveInverse() {
    for (int k = 0; k < dim; ++k) {
      for (int l = 0; l < dim; ++l) {
        TCOMPONENT a(k, l, scalar);
        TENSOR A(a);
        EXPECT_EQ(TENSOR(-a), -A);
      }
    }
  }

  virtual void checkScalarMultiplication() {
    for (int k = 0; k < dim; ++k) {
      for (int l = 0; l < dim; ++l) {
        TCOMPONENT v(k, l, scalar);

        auto r1 = v * scalar;
        auto r2 = scalar * v;
        auto r3 = v * scalar_i;
        auto r4 = scalar_i * v;

        EXPECT_EQ(r1.Value(), scalar * scalar);
        EXPECT_EQ(r2.Value(), scalar * scalar);
        EXPECT_EQ(r3.Value(), scalar * scalar_i);
        EXPECT_EQ(r4.Value(), scalar * scalar_i);
        EXPECT_EQ(r1.Row(), k);
        EXPECT_EQ(r2.Row(), k);
        EXPECT_EQ(r3.Row(), k);
        EXPECT_EQ(r4.Row(), k);
        EXPECT_EQ(r1.Column(), l);
        EXPECT_EQ(r2.Column(), l);
        EXPECT_EQ(r3.Column(), l);
        EXPECT_EQ(r4.Column(), l);
      }
    }
  }

  virtual void checkScalarDivision() {
    for (int k = 0; k < dim; ++k) {
      for (int l = 0; l < dim; ++l) {
        TCOMPONENT v(k, l, scalar);

        auto r1 = v / scalar;
        auto r3 = v / scalar_i;

        EXPECT_EQ(r1.Value(), scalar / scalar);
        EXPECT_EQ(r3.Value(), scalar / scalar_i);
        EXPECT_EQ(r1.Row(), k);
        EXPECT_EQ(r3.Row(), k);
        EXPECT_EQ(r1.Column(), l);
        EXPECT_EQ(r3.Column(), l);
      }
    }
  }

  virtual void checkVectorTensorProduct() {
    for (int k = 0; k < dim; ++k) {
      for (int l = 0; l < dim; ++l) {
        TCOMPONENT A(k, l, scalar);
        TENSOR B(A);

        EXPECT_EQ(VECTORFIELD(A * U), B * U);
        EXPECT_EQ(VECTORFIELD(A.multiply(U)), B * U);
        EXPECT_EQ(VECTORFIELD(U * A), U * B);
        EXPECT_EQ(VECTORFIELD(A.applyTransposed(U)), U * B);

        for (int j = 0; j < dim; ++j) {
          VCOMPONENT v(j, scalar);

          EXPECT_EQ(VECTORFIELD(A * v), B * v);
          EXPECT_EQ(VECTORFIELD(A.multiply(v)), B * v);
          EXPECT_EQ(VECTORFIELD(v * A), v * B);
          EXPECT_EQ(VECTORFIELD(A.applyTransposed(v)), v * B);
        }
      }
    }
  }

  virtual void checkTensorTensorProduct() {
    for (int k = 0; k < dim; ++k) {
      for (int l = 0; l < dim; ++l) {
        TCOMPONENT a(k, l, scalar);
        TENSOR A(a);
        for (int i = 0; i < dim; ++i) {
          for (int j = 0; j < dim; ++j) {
            TCOMPONENT b(i, j, scalar);
            TENSOR B(b);
            EXPECT_EQ(TENSOR(a * b), A * B);
          }
        }
      }
    }
  }

  virtual void checkFrobenius() {
    for (int k = 0; k < dim; ++k) {
      for (int l = 0; l < dim; ++l) {
        TCOMPONENT a(k, l, scalar);
        TENSOR A(a);
        for (int i = 0; i < dim; ++i) {
          for (int j = 0; j < dim; ++j) {
            TCOMPONENT b(i, j, scalar);
            TENSOR B(b);
            EXPECT_DOUBLE_EQ(Frobenius(a, b), Frobenius(A, B));
          }
        }
      }
    }
  }

  virtual void checkNorm() {
    for (int k = 0; k < dim; ++k) {
      for (int l = 0; l < dim; ++l) {
        TCOMPONENT a(k, l, scalar);
        TENSOR A(a);
        EXPECT_DOUBLE_EQ(norm(a), norm(A));
      }
    }
  }

  virtual void checkTrace() {
    for (int k = 0; k < dim; ++k) {
      for (int l = 0; l < dim; ++l) {
        TCOMPONENT a(k, l, scalar);
        TENSOR A(a);
        EXPECT_DOUBLE_EQ(trace(a), trace(A));
      }
    }
  }

  virtual void checkEQ() {
    for (int k = 0; k < dim; ++k) {
      for (int l = 0; l < dim; ++l) {
        TCOMPONENT a(k, l, scalar);
        EXPECT_EQ(a, a);
      }
    }
  }

  virtual void checkNonEQ() {
    std::vector<TCOMPONENT> A;
    for (int k = 0; k < dim; ++k) {
      for (int l = 0; l < dim; ++l) {
        A.push_back(TCOMPONENT(k, l, scalar));
      }
    }
    for (int k = 0; k < dim; ++k) {
      for (int l = 0; l < dim; ++l) {
        if (k != l) EXPECT_NE(A[k], A[l]);
      }
    }
  }

  virtual void checkSaveLoad() {
    for (int k = 0; k < dim; ++k) {
      for (int l = 0; l < dim; ++l) {
        TCOMPONENT v(k, l, scalar);
        Saver saver("SaveLoadTest");
        saver << v;
        saver.close();
        Loader loader("SaveLoadTest");
        TCOMPONENT w;
        loader >> w;
        loader.close();
        EXPECT_EQ(v, w);
      }
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

using TensorComponentTest_3 = TensorComponentTest<SpaceDimension>;

TENSOR_ROW_TESTS(TensorComponentTest_3)

/*
using TensorComponentTest_2 = TensorComponentTest<2>;

TENSOR_ROW_TESTS(TensorComponentTest_2)

using TensorComponentTest_1 = TensorComponentTest<1>;

TENSOR_ROW_TESTS(TensorComponentTest_1)*/

int main(int argc, char **argv) {
  MppTest mppTest = MppTestBuilder(argc, argv).WithPPM().WithScreenLogging();
  return mppTest.RUN_ALL_MPP_TESTS();
}
