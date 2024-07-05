#include "Tensor.hpp"
#include "TestEnvironment.hpp"

#include <vector>

template<int dim>
class VectorFieldComponentTest : public ::testing::Test {
protected:
  using VECTORFIELD = VectorFieldT<double, dim>;
  using TENSOR = TensorT<double, dim>;
  using COMPONENT = VectorFieldComponentT<double>;

  double value;
  double scalar;
  int scalar_i;

  VECTORFIELD U;
  TENSOR T;

  VectorFieldComponentTest() {
    mpp_geometry::SetTolerance(1e-15);
    value = RandomNonZeroDouble();
    scalar = RandomNonZeroDouble();
    scalar_i = RandomNonZeroInt();
    for (int i = 0; i < dim; ++i) {
      U[i] = RandomNonZeroDouble();
      for (int j = 0; j < dim; ++j)
        T[i][j] = RandomNonZeroDouble();
    }
  }

  void checkConstructor() {
    COMPONENT z;
    EXPECT_DOUBLE_EQ(z.Value(), 0.0);
    EXPECT_EQ(z.Component(), 0);

    for (int k = 0; k < dim; ++k) {
      COMPONENT v(k, value);
      EXPECT_DOUBLE_EQ(v.Value(), value);
      EXPECT_EQ(v.Component(), k);
      for (int i = 0; i < dim; ++i) {
        if (i == k) EXPECT_DOUBLE_EQ(value, v[i]);
        else EXPECT_DOUBLE_EQ(0.0, v[i]);
      }

      COMPONENT u(v);
      EXPECT_DOUBLE_EQ(u.Value(), value);
      EXPECT_EQ(u.Component(), k);

      VECTORFIELD w(v);
      for (int i = 0; i < dim; ++i) {
        if (i == k) EXPECT_DOUBLE_EQ(value, w[i]);
        else EXPECT_DOUBLE_EQ(0.0, w[i]);
      }
    }
  }

  virtual void checkAssignment() {
    for (int k = 0; k < dim; ++k) {
      COMPONENT v(k, value);
      COMPONENT u = v;
      EXPECT_DOUBLE_EQ(u.Value(), value);
      EXPECT_EQ(u.Component(), k);
    }
  }

  virtual void checkSummation() {
    for (int k = 0; k < dim; ++k) {
      COMPONENT v(k, value);
      VECTORFIELD V(v);
      EXPECT_EQ(U + v, U + V);
      EXPECT_EQ(v + U, V + U);
    }
  }

  virtual void checkDifference() {
    for (int k = 0; k < dim; ++k) {
      COMPONENT v(k, value);
      VECTORFIELD V(v);
      EXPECT_EQ(U - v, U - V);
      EXPECT_EQ(v - U, V - U);
    }
  }

  virtual void checkAdditiveInverse() {
    for (int k = 0; k < dim; ++k) {
      COMPONENT v(k, value);
      VECTORFIELD V(v);
      EXPECT_EQ(VECTORFIELD(-v), -V);
    }
  }

  virtual void checkScalarMultiplication() {
    for (int k = 0; k < dim; ++k) {
      COMPONENT v(k, value);

      auto r1 = v * scalar;
      auto r2 = scalar * v;
      auto r3 = v * scalar_i;
      auto r4 = scalar_i * v;

      EXPECT_DOUBLE_EQ(r1.Value(), value * scalar);
      EXPECT_DOUBLE_EQ(r2.Value(), value * scalar);
      EXPECT_DOUBLE_EQ(r3.Value(), value * scalar_i);
      EXPECT_DOUBLE_EQ(r4.Value(), value * scalar_i);
      EXPECT_EQ(r1.Component(), k);
      EXPECT_EQ(r2.Component(), k);
      EXPECT_EQ(r3.Component(), k);
      EXPECT_EQ(r4.Component(), k);
    }
  }

  virtual void checkScalarDivision() {
    for (int k = 0; k < dim; ++k) {
      COMPONENT v(k, value);

      auto r1 = v / scalar;
      auto r3 = v / scalar_i;

      EXPECT_DOUBLE_EQ(r1.Value(), value / scalar);
      EXPECT_DOUBLE_EQ(r3.Value(), value / scalar_i);
      EXPECT_EQ(r1.Component(), k);
      EXPECT_EQ(r3.Component(), k);
    }
  }

  virtual void checkVectorTensorProduct() {
    for (int k = 0; k < dim; ++k) {
      COMPONENT v(k, value);
      VECTORFIELD w(v);
      mout << DOUT(T) << DOUT(w) << endl;

      auto res = T * v;

      mout << DOUT(res) << endl;

      EXPECT_EQ(T * v, T * w);
      EXPECT_EQ(T.multiply(v), T * w);
      EXPECT_EQ(v * T, w * T);
      EXPECT_EQ(T.applyTransposed(v), w * T);
    }
  }

  virtual void checkScalarProduct() {
    for (int k = 0; k < dim; ++k) {
      COMPONENT v(k, value);

      EXPECT_DOUBLE_EQ(v * U, value * U[k]);
      EXPECT_DOUBLE_EQ(U * v, value * U[k]);
      for (int l = 0; l < dim; ++l) {
        COMPONENT w(l, scalar);
        if (k == l) {
          EXPECT_DOUBLE_EQ(v * w, value * scalar);
        } else {
          EXPECT_DOUBLE_EQ(v * w, 0.0);
        }
      }
    }
  }

  virtual void checkNorm() {
    for (int k = 0; k < dim; ++k) {
      COMPONENT v(k, value);
      EXPECT_DOUBLE_EQ(norm(v), std::abs(value));
    }
  }

  virtual void checkEQ() {
    for (int k = 0; k < dim; ++k) {
      COMPONENT v(k, value);
      EXPECT_EQ(v, v);
    }
  }

  virtual void checkNonEQ() {
    std::vector<COMPONENT> A;
    for (int k = 0; k < dim; ++k) {
      A.push_back(COMPONENT(k, value));
    }
    for (int k = 0; k < dim; ++k) {
      for (int l = 0; l < dim; ++l) {
        if (k != l) EXPECT_NE(A[k], A[l]);
      }
    }
  }

  virtual void checkSaveLoad() {
    for (int k = 0; k < dim; ++k) {
      COMPONENT v(k, value);
      Saver saver("SaveLoadTest");
      saver << v;
      saver.close();
      Loader loader("SaveLoadTest");
      COMPONENT w;
      loader >> w;
      loader.close();
      EXPECT_EQ(v, w);
    }
  }
};

/// To avoid the same code in each test
#define VECTORFIELD_COMPONENT_TESTS(vectorfieldClass)                                              \
                                                                                                   \
  TEST_F(vectorfieldClass, ConstructorTest) { checkConstructor(); }                                \
                                                                                                   \
  TEST_F(vectorfieldClass, AssignmentTest) { checkAssignment(); }                                  \
                                                                                                   \
  TEST_F(vectorfieldClass, SummationTest) { checkSummation(); }                                    \
                                                                                                   \
  TEST_F(vectorfieldClass, DifferenceTest) { checkDifference(); }                                  \
                                                                                                   \
  TEST_F(vectorfieldClass, AdditiveInverseTest) { checkAdditiveInverse(); }                        \
                                                                                                   \
  TEST_F(vectorfieldClass, ScalarMultiplicationTest) { checkScalarMultiplication(); }              \
                                                                                                   \
  TEST_F(vectorfieldClass, ScalarDivisionTest) { checkScalarDivision(); }                          \
                                                                                                   \
  TEST_F(vectorfieldClass, VectorTensorProductTest) { checkVectorTensorProduct(); }                \
                                                                                                   \
  TEST_F(vectorfieldClass, ScalarProductTest) { checkScalarProduct(); }                            \
                                                                                                   \
  TEST_F(vectorfieldClass, NormTest) { checkNorm(); }                                              \
                                                                                                   \
  TEST_F(vectorfieldClass, EqualTest) { checkEQ(); }                                               \
                                                                                                   \
  TEST_F(vectorfieldClass, NonEqualTest) { checkNonEQ(); }                                         \
                                                                                                   \
  TEST_F(vectorfieldClass, SaveLoadTest) { checkSaveLoad(); }

using VectorFieldComponentTest_3 = VectorFieldComponentTest<SpaceDimension>;

VECTORFIELD_COMPONENT_TESTS(VectorFieldComponentTest_3)

/*
using VectorFieldComponentTest_2 = VectorFieldComponentTest<2>;

VECTORFIELD_COMPONENT_TESTS(VectorFieldComponentTest_2)

using VectorFieldComponentTest_1 = VectorFieldComponentTest<1>;

VECTORFIELD_COMPONENT_TESTS(VectorFieldComponentTest_1)
*/
int main(int argc, char **argv) {
  MppTest mppTest = MppTestBuilder(argc, argv).WithPPM().WithScreenLogging();
  return mppTest.RUN_ALL_MPP_TESTS();
}
