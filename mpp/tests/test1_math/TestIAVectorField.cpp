#include "TestEnvironment.hpp"
#include "VectorField.hpp"

#include <vector>

template<int sDim, int tDim>
class IAVectorFieldTest : public ::testing::Test {
protected:
  using IAVECTORFIELD = VectorFieldT<IAInterval, sDim>;
  using IAPOINT = PointT<IAInterval, sDim, tDim>;

  int dim = -1;
  IAVECTORFIELD p, q;
  IAPOINT z;
  IAInterval value;
  double value_d;
  int value_i;

  std::vector<IAInterval> values;
  std::vector<IAVECTORFIELD> vectorfields;

  IAVectorFieldTest() {
    mpp_geometry::SetTolerance(1e-15);
    dim = p.Dim();
    for (int i = 0; i < dim; ++i) {
      p[i] = RandomIAInterval();
      q[i] = RandomIAInterval();
    }
    for (int i = 0; i < z.SpaceDim() + z.TimeDim(); ++i)
      z[i] = RandomIAInterval();
    value = RandomNonZeroIAInterval();
    value_d = RandomNonZeroDouble();
    value_i = RandomNonZeroInt();

    IAVECTORFIELD tmp;
    values.resize(dim);
    vectorfields.resize(dim);
    for (int i = 0; i < dim; ++i)
      values[i] = RandomIAInterval();
    if constexpr (tmp.Dim() >= 2) vectorfields[1] = IAVECTORFIELD(values[0], values[1]);
    if constexpr (tmp.Dim() >= 3) vectorfields[2] = IAVECTORFIELD(values[0], values[1], values[2]);
  }

  void checkConstructor() {
    IAVECTORFIELD z0;
    IAVECTORFIELD z1(p);
    IAVECTORFIELD z2(mid(p));
    IAVECTORFIELD z3(value);
    IAVECTORFIELD z4(value_d);
    IAVECTORFIELD z5(value_i);
    IAVECTORFIELD z6(z);
    for (int i = 0; i < dim; ++i) {
      EXPECT_EQ(z0[i], IAInterval());
      EXPECT_EQ(z1[i], p[i]);
      EXPECT_EQ(z2[i], IAInterval(mid(p[i])));
      EXPECT_EQ(z3[i], value);
      EXPECT_EQ(z4[i], IAInterval(value_d));
      EXPECT_EQ(z5[i], IAInterval(value_i));
      EXPECT_EQ(z6[i], z[i]);
    }

    for (int i = 1; i < dim; ++i) {
      for (int k = 0; k < i + 1; ++k)
        EXPECT_EQ(vectorfields[i][k], values[k]);
      for (int k = i + 1; k < dim; ++k)
        EXPECT_EQ(vectorfields[i][k], IAInterval());
    }

    for (int i = 0; i < dim; ++i) {
      IAVECTORFIELD z(i, value);
      for (int k = 0; k < dim; ++k)
        if (i == k) EXPECT_EQ(z[k], value);
        else EXPECT_EQ(z[k], IAInterval());
    }
  }

  virtual void checkAssignment() {
    IAVECTORFIELD z1 = p;
    IAVECTORFIELD z2 = mid(p);
    IAVECTORFIELD z3(p);
    z3 = value;
    IAVECTORFIELD z4(p);
    z4 = value_d;
    IAVECTORFIELD z5(p);
    z5 = value_i;
    IAVECTORFIELD z6 = z;
    for (int i = 0; i < dim; ++i) {
      EXPECT_EQ(z1[i], p[i]);
      EXPECT_EQ(z2[i], IAInterval(mid(p[i])));
      EXPECT_EQ(z3[i], value);
      EXPECT_EQ(z4[i], IAInterval(value_d));
      EXPECT_EQ(z5[i], IAInterval(value_i));
      EXPECT_EQ(z6[i], z[i]);
    }
  }

  virtual void checkSummation() {
    auto r0 = p + q;
    auto r1 = mid(p) + q;
    auto r2 = p + mid(q);
    auto r3 = p + z;
    auto r4 = z + q;
    for (int i = 0; i < dim; ++i) {
      EXPECT_EQ(r0[i], p[i] + q[i]);
      EXPECT_EQ(r1[i], mid(p[i]) + q[i]);
      EXPECT_EQ(r2[i], p[i] + mid(q[i]));
      EXPECT_EQ(r3[i], p[i] + z[i]);
      EXPECT_EQ(r4[i], z[i] + q[i]);
    }
  }

  virtual void checkDifference() {
    auto r0 = p - q;
    auto r1 = mid(p) - q;
    auto r2 = p - mid(q);
    auto r3 = p - z;
    auto r4 = z - q;
    for (int i = 0; i < dim; ++i) {
      EXPECT_EQ(r0[i], p[i] - q[i]);
      EXPECT_EQ(r1[i], mid(p[i]) - q[i]);
      EXPECT_EQ(r2[i], p[i] - mid(q[i]));
      EXPECT_EQ(r3[i], p[i] - z[i]);
      EXPECT_EQ(r4[i], z[i] - q[i]);
    }
  }

  virtual void checkAdditiveInverse() {
    auto z = -p;
    for (int i = 0; i < dim; ++i)
      EXPECT_EQ(z[i], -p[i]);
  }

  virtual void checkScalarMultiplication() {
    auto r1 = p * value;
    auto r2 = value * p;
    auto r3 = p * value_d;
    auto r4 = value_d * p;
    auto r5 = p * value_i;
    auto r6 = value_i * p;
    auto r7 = value * mid(p);
    auto r8 = mid(p) * value;
    for (int i = 0; i < dim; ++i) {
      EXPECT_EQ(r1[i], p[i] * value);
      EXPECT_EQ(r2[i], p[i] * value);
      EXPECT_EQ(r3[i], p[i] * value_d);
      EXPECT_EQ(r4[i], p[i] * value_d);
      EXPECT_EQ(r5[i], p[i] * value_i);
      EXPECT_EQ(r6[i], p[i] * value_i);
      EXPECT_EQ(r7[i], mid(p[i]) * value);
      EXPECT_EQ(r8[i], mid(p[i]) * value);
    }
  }

  virtual void checkScalarDivision() {
    auto r1 = p / value;
    auto r2 = p / value_d;
    auto r3 = p / value_i;
    auto r4 = mid(p) / value;
    for (int i = 0; i < dim; ++i) {
      EXPECT_EQ(r1[i], p[i] / value);
      EXPECT_EQ(r2[i], p[i] / value_d);
      EXPECT_EQ(r3[i], p[i] / value_i);
      EXPECT_EQ(r4[i], mid(p[i]) / value);
    }
  }

  virtual void checkScalarProduct() {
    IAInterval p0, p1, p2, p3, p4;
    for (int i = 0; i < dim; ++i) {
      p0 += p[i] * q[i];
      p1 += mid(p[i]) * q[i];
      p2 += p[i] * mid(q[i]);
      p3 += p[i] * z[i];
      p4 += z[i] * q[i];
    }
    EXPECT_EQ(p0, p * q);
    EXPECT_EQ(p1, mid(p) * q);
    EXPECT_EQ(p2, p * mid(q));
    EXPECT_EQ(p3, p * z);
    EXPECT_EQ(p4, z * q);
  }

  virtual void checkNorm() {
    IAInterval sum;
    for (int i = 0; i < dim; ++i)
      sum += sqr(p[i]);
    EXPECT_TRUE(norm(p) <= sqrt(sum));
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
    IAPOINT r;
    VFtoPoint(p, r);
    for (int i = 0; i < r.SpaceDim(); ++i)
      EXPECT_EQ(r[i], p[i]);
    if (r.TimeDim() != 0) EXPECT_EQ(r[r.SpaceDim()], IAInterval());
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
    IAVECTORFIELD z;
    loader >> z;
    loader.close();
    EXPECT_EQ(p, z);
  }

  virtual void checkMid() {
    auto r = mid(p);
    for (int i = 0; i < dim; ++i)
      EXPECT_DOUBLE_EQ(r[i], mid(p[i]));
  }

  virtual void checkInf() {
    auto r = inf(p);
    for (int i = 0; i < dim; ++i)
      EXPECT_DOUBLE_EQ(r[i], inf(p[i]));
  }

  virtual void checkSup() {
    auto r = sup(p);
    for (int i = 0; i < dim; ++i)
      EXPECT_DOUBLE_EQ(r[i], sup(p[i]));
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
  TEST_F(vectorfieldClass, SaveLoadTest) { checkSaveLoad(); }                                      \
                                                                                                   \
  TEST_F(vectorfieldClass, MidTest) { checkMid(); }                                                \
                                                                                                   \
  TEST_F(vectorfieldClass, InfTest) { checkInf(); }                                                \
                                                                                                   \
  TEST_F(vectorfieldClass, SupTest) { checkSup(); }

template<int tDim>
class IAVectorFieldTest_3_ : public IAVectorFieldTest<3, tDim> {
protected:
  void checkCrossProduct() override {
    VectorFieldT<IAInterval, 3> z0;
    z0[0] = this->p[1] * this->q[2] - this->p[2] * this->q[1];
    z0[1] = this->p[2] * this->q[0] - this->p[0] * this->q[2];
    z0[2] = this->p[0] * this->q[1] - this->p[1] * this->q[0];
    EXPECT_EQ(z0, this->p ^ this->q);
    EXPECT_EQ(z0, curl(this->p, this->q));
    VectorFieldT<IAInterval, 3> z1;
    z1[0] = this->p[1] * this->z[2] - this->p[2] * this->z[1];
    z1[1] = this->p[2] * this->z[0] - this->p[0] * this->z[2];
    z1[2] = this->p[0] * this->z[1] - this->p[1] * this->z[0];
    EXPECT_EQ(z1, this->p ^ this->z);
    EXPECT_EQ(z1, curl(this->p, this->z));
    VectorFieldT<IAInterval, 3> z2;
    z2[0] = this->z[1] * this->q[2] - this->z[2] * this->q[1];
    z2[1] = this->z[2] * this->q[0] - this->z[0] * this->q[2];
    z2[2] = this->z[0] * this->q[1] - this->z[1] * this->q[0];
    EXPECT_EQ(z2, this->z ^ this->q);
    EXPECT_EQ(z2, curl(this->z, this->q));
  }
};

using IAVectorFieldTest_3_1 = IAVectorFieldTest_3_<1>;

VECTORFIELD_TESTS(IAVectorFieldTest_3_1)

using IAVectorFieldTest_2_1 = IAVectorFieldTest<2, 1>;

VECTORFIELD_TESTS(IAVectorFieldTest_2_1)

using IAVectorFieldTest_1_1 = IAVectorFieldTest<1, 1>;

VECTORFIELD_TESTS(IAVectorFieldTest_1_1)

using IAVectorFieldTest_3_0 = IAVectorFieldTest_3_<0>;

VECTORFIELD_TESTS(IAVectorFieldTest_3_0)

using IAVectorFieldTest_2_0 = IAVectorFieldTest<2, 0>;

VECTORFIELD_TESTS(IAVectorFieldTest_2_0)

using IAVectorFieldTest_1_0 = IAVectorFieldTest<1, 0>;

VECTORFIELD_TESTS(IAVectorFieldTest_1_0)

int main(int argc, char **argv) {
  MppTest mppTest = MppTestBuilder(argc, argv).WithPPM();
  return mppTest.RUN_ALL_MPP_TESTS();
}
