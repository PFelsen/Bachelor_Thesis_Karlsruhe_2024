#include "TestEnvironment.hpp"
#include "VectorField.hpp"

#include <vector>

template<int sDim, int tDim>
class VectorFieldTest : public ::testing::Test {
protected:
  using VECTORFIELD = VectorFieldT<double, sDim>;
  using POINT = PointT<double, sDim, tDim>;

  int dim = -1;
  VECTORFIELD p, q;
  POINT z;
  double value;
  int value_i;

  std::vector<double> values;
  std::vector<VECTORFIELD> vectorfields;

  VectorFieldTest() {
    mpp_geometry::SetTolerance(1e-15);
    dim = p.Dim();
    for (int i = 0; i < dim; ++i) {
      p[i] = RandomDouble();
      q[i] = RandomDouble();
    }
    for (int i = 0; i < z.SpaceDim() + z.TimeDim(); ++i)
      z[i] = RandomDouble();
    value = RandomNonZeroDouble();
    value_i = RandomNonZeroInt();

    VECTORFIELD tmp;
    values.resize(dim);
    vectorfields.resize(dim);
    for (int i = 0; i < dim; ++i)
      values[i] = RandomDouble();
    if constexpr (tmp.Dim() >= 2) vectorfields[1] = VECTORFIELD(values[0], values[1]);
    if constexpr (tmp.Dim() >= 3) vectorfields[2] = VECTORFIELD(values[0], values[1], values[2]);
  }

  void checkConstructor() {
    VECTORFIELD z0;
    VECTORFIELD z1(p);
    VECTORFIELD z4(value);
    VECTORFIELD z5(value_i);
    VECTORFIELD z6(z);
    for (int i = 0; i < dim; ++i) {
      EXPECT_DOUBLE_EQ(z0[i], 0.0);
      EXPECT_DOUBLE_EQ(z1[i], p[i]);
      EXPECT_DOUBLE_EQ(z4[i], value);
      EXPECT_DOUBLE_EQ(z5[i], double(value_i));
      EXPECT_DOUBLE_EQ(z6[i], z[i]);
    }

    for (int i = 1; i < dim; ++i) {
      for (int k = 0; k < i + 1; ++k)
        EXPECT_DOUBLE_EQ(vectorfields[i][k], values[k]);
      for (int k = i + 1; k < dim; ++k)
        EXPECT_DOUBLE_EQ(vectorfields[i][k], 0.0);
    }

    for (int i = 0; i < dim; ++i) {
      VECTORFIELD z(i, value);
      for (int k = 0; k < dim; ++k)
        if (i == k) EXPECT_DOUBLE_EQ(z[k], value);
        else EXPECT_DOUBLE_EQ(z[k], 0.0);
    }
  }

  virtual void checkAssignment() {
    VECTORFIELD z1 = p;
    VECTORFIELD z4(p);
    z4 = value;
    VECTORFIELD z5(p);
    z5 = value_i;
    VECTORFIELD z6 = z;
    for (int i = 0; i < dim; ++i) {
      EXPECT_DOUBLE_EQ(z1[i], p[i]);
      EXPECT_DOUBLE_EQ(z4[i], value);
      EXPECT_DOUBLE_EQ(z5[i], double(value_i));
      EXPECT_DOUBLE_EQ(z6[i], z[i]);
    }
  }

  virtual void checkSummation() {
    auto r0 = p + q;
    auto r3 = p + z;
    auto r4 = z + q;
    for (int i = 0; i < dim; ++i) {
      EXPECT_DOUBLE_EQ(r0[i], p[i] + q[i]);
      EXPECT_DOUBLE_EQ(r3[i], p[i] + z[i]);
      EXPECT_DOUBLE_EQ(r4[i], z[i] + q[i]);
    }
  }

  virtual void checkDifference() {
    auto r0 = p - q;
    auto r3 = p - z;
    auto r4 = z - q;
    for (int i = 0; i < dim; ++i) {
      EXPECT_DOUBLE_EQ(r0[i], p[i] - q[i]);
      EXPECT_DOUBLE_EQ(r3[i], p[i] - z[i]);
      EXPECT_DOUBLE_EQ(r4[i], z[i] - q[i]);
    }
  }

  virtual void checkAdditiveInverse() {
    auto z = -p;
    for (int i = 0; i < dim; ++i)
      EXPECT_DOUBLE_EQ(z[i], -p[i]);
  }

  virtual void checkScalarMultiplication() {
    auto r1 = p * value;
    auto r2 = value * p;
    auto r5 = p * value_i;
    auto r6 = value_i * p;
    for (int i = 0; i < dim; ++i) {
      EXPECT_DOUBLE_EQ(r1[i], p[i] * value);
      EXPECT_DOUBLE_EQ(r2[i], p[i] * value);
      EXPECT_DOUBLE_EQ(r5[i], p[i] * value_i);
      EXPECT_DOUBLE_EQ(r6[i], p[i] * value_i);
    }
  }

  virtual void checkScalarDivision() {
    auto r1 = p / value;
    auto r3 = p / value_i;
    for (int i = 0; i < dim; ++i) {
      EXPECT_DOUBLE_EQ(r1[i], p[i] / value);
      EXPECT_DOUBLE_EQ(r3[i], p[i] / value_i);
    }
  }

  virtual void checkScalarProduct() {
    double p0{}, p3{}, p4{};
    for (int i = 0; i < dim; ++i) {
      p0 += p[i] * q[i];
      p3 += p[i] * z[i];
      p4 += z[i] * q[i];
    }
    EXPECT_DOUBLE_EQ(p0, p * q);
    EXPECT_DOUBLE_EQ(p3, p * z);
    EXPECT_DOUBLE_EQ(p4, z * q);
  }

  virtual void checkNorm() {
    double sum = 0.0;
    for (int i = 0; i < dim; ++i)
      sum += p[i] * p[i];
    EXPECT_DOUBLE_EQ(norm(p), sqrt(sum));
  }

  virtual void checkCrossProduct() {
    EXPECT_ANY_THROW(p ^ q);
    EXPECT_ANY_THROW(curl(p, q));
    EXPECT_ANY_THROW(p ^ z);
    EXPECT_ANY_THROW(curl(p, z));
    EXPECT_ANY_THROW(z ^ q);
    EXPECT_ANY_THROW(curl(z, q));
  };

  virtual void checkVFtoPoint() {
    POINT r;
    VFtoPoint(p, r);
    for (int i = 0; i < r.SpaceDim(); ++i)
      EXPECT_DOUBLE_EQ(r[i], p[i]);
    if (r.TimeDim() != 0) EXPECT_DOUBLE_EQ(r[r.SpaceDim()], 0.0);
  }

  virtual void checkEQ() {
    EXPECT_EQ(p, p);
    EXPECT_EQ(q, q);
  }

  virtual void checkNonEQ() {
    EXPECT_NE(p, q);
    EXPECT_NE(q, p);
  }

  virtual void checkSaveLoad() {
    Saver saver("SaveLoadTest");
    saver << p;
    saver.close();
    Loader loader("SaveLoadTest");
    VECTORFIELD z;
    loader >> z;
    loader.close();
    EXPECT_EQ(p, z);
  }
};

/// To avoid the same code in each test
#define VECTORFIELD_TESTS(vectorfieldClass)                                                        \
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
  TEST_F(vectorfieldClass, ScalarProductTest) { checkScalarProduct(); }                            \
                                                                                                   \
  TEST_F(vectorfieldClass, NormTest) { checkNorm(); }                                              \
                                                                                                   \
  TEST_F(vectorfieldClass, CrossProductTest) { checkCrossProduct(); }                              \
                                                                                                   \
  TEST_F(vectorfieldClass, VectorFieldToPointTest) { checkVFtoPoint(); }                           \
                                                                                                   \
  TEST_F(vectorfieldClass, EqualTest) { checkEQ(); }                                               \
                                                                                                   \
  TEST_F(vectorfieldClass, NonEqualTest) { checkNonEQ(); }                                         \
                                                                                                   \
  TEST_F(vectorfieldClass, SaveLoadTest) { checkSaveLoad(); }

template<int tDim>
class VectorFieldTest_3_ : public VectorFieldTest<3, tDim> {
protected:
  void checkCrossProduct() override {
    VectorFieldT<double, 3> z0;
    z0[0] = this->p[1] * this->q[2] - this->p[2] * this->q[1];
    z0[1] = this->p[2] * this->q[0] - this->p[0] * this->q[2];
    z0[2] = this->p[0] * this->q[1] - this->p[1] * this->q[0];
    EXPECT_EQ(z0, this->p ^ this->q);
    EXPECT_EQ(z0, curl(this->p, this->q));
    VectorFieldT<double, 3> z1;
    z1[0] = this->p[1] * this->z[2] - this->p[2] * this->z[1];
    z1[1] = this->p[2] * this->z[0] - this->p[0] * this->z[2];
    z1[2] = this->p[0] * this->z[1] - this->p[1] * this->z[0];
    EXPECT_EQ(z1, this->p ^ this->z);
    EXPECT_EQ(z1, curl(this->p, this->z));
    VectorFieldT<double, 3> z2;
    z2[0] = this->z[1] * this->q[2] - this->z[2] * this->q[1];
    z2[1] = this->z[2] * this->q[0] - this->z[0] * this->q[2];
    z2[2] = this->z[0] * this->q[1] - this->z[1] * this->q[0];
    EXPECT_EQ(z2, this->z ^ this->q);
    EXPECT_EQ(z2, curl(this->z, this->q));
  }
};

using VectorFieldTest_3_1 = VectorFieldTest_3_<1>;

VECTORFIELD_TESTS(VectorFieldTest_3_1)

using VectorFieldTest_2_1 = VectorFieldTest<2, 1>;

VECTORFIELD_TESTS(VectorFieldTest_2_1)

using VectorFieldTest_1_1 = VectorFieldTest<1, 1>;

VECTORFIELD_TESTS(VectorFieldTest_1_1)

using VectorFieldTest_3_0 = VectorFieldTest_3_<0>;

VECTORFIELD_TESTS(VectorFieldTest_3_0)

using VectorFieldTest_2_0 = VectorFieldTest<2, 0>;

VECTORFIELD_TESTS(VectorFieldTest_2_0)

using VectorFieldTest_1_0 = VectorFieldTest<1, 0>;

VECTORFIELD_TESTS(VectorFieldTest_1_0)

int main(int argc, char **argv) {
  MppTest mppTest = MppTestBuilder(argc, argv).WithPPM();
  return mppTest.RUN_ALL_MPP_TESTS();
}
