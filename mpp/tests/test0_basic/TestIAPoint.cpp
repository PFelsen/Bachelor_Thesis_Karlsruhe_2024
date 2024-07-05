#include "Point.hpp"
#include "TestEnvironment.hpp"

#include <vector>

template<typename IAPOINT>
class IAPointTest : public ::testing::Test {
protected:
  int size = -1;
  int size_Space = -1;
  IAPOINT p, q;
  IAInterval value;
  double value_d;
  int value_i;

  bool timeDep = false;

  std::vector<IAInterval> values;
  std::vector<IAPOINT> points;

  IAPointTest() {
    mpp_geometry::SetTolerance(1e-15);
    size = p.SpaceDim() + p.TimeDim();
    size_Space = p.SpaceDim();
    for (int i = 0; i < size; ++i) {
      p[i] = RandomIAInterval();
      q[i] = RandomIAInterval();
    }
    value = RandomNonZeroIAInterval();
    value_d = RandomNonZeroDouble();
    value_i = RandomNonZeroInt();
    timeDep = p.TimeDim() == 1;

    IAPOINT tmp;
    values.resize(size);
    points.resize(size);
    for (int i = 0; i < size; ++i)
      values[i] = RandomIAInterval();
    points[0] = IAPOINT(values[0]);
    if constexpr (tmp.SpaceDim() + tmp.TimeDim() >= 2) points[1] = IAPOINT(values[0], values[1]);
    if constexpr (tmp.SpaceDim() + tmp.TimeDim() >= 3)
      points[2] = IAPOINT(values[0], values[1], values[2]);
    if constexpr (tmp.SpaceDim() + tmp.TimeDim() >= 4)
      points[3] = IAPOINT(values[0], values[1], values[2], values[3]);
  }

  void checkConstructor() {
    IAPOINT z0;
    IAPOINT z1(p);
    IAPOINT z2(mid(p));
    for (int i = 0; i < size; ++i) {
      EXPECT_EQ(z0[i], IAInterval());
      EXPECT_EQ(z1[i], p[i]);
      EXPECT_EQ(z2[i], IAInterval(mid(p[i])));
    }
    for (int i = 0; i < size; ++i) {
      for (int k = 0; k < i + 1; ++k)
        EXPECT_EQ(points[i][k], values[k]);
      for (int k = i + 1; k < size; ++k)
        EXPECT_EQ(points[i][k], IAInterval());
    }
  }

  virtual void checkAssignment() {
    IAPOINT z1 = p;
    IAPOINT z2 = mid(p);
    IAPOINT z3(p);
    z3 = value;
    IAPOINT z4(p);
    z4 = value_d;
    IAPOINT z5(p);
    z5 = value_i;
    for (int i = 0; i < size; ++i) {
      EXPECT_EQ(z1[i], p[i]);
      EXPECT_EQ(z2[i], IAInterval(mid(p[i])));
      EXPECT_EQ(z3[i], value);
      EXPECT_EQ(z4[i], IAInterval(value_d));
      EXPECT_EQ(z5[i], IAInterval(value_i));
    }
  }

  virtual void checkSummation() {
    auto r0 = p + q;
    auto r1 = mid(p) + q;
    auto r2 = p + mid(q);
    for (int i = 0; i < size_Space; ++i) {
      EXPECT_EQ(r0[i], p[i] + q[i]);
      EXPECT_EQ(r1[i], mid(p[i]) + q[i]);
      EXPECT_EQ(r2[i], p[i] + mid(q[i]));
    }
  }

  virtual void checkDifference() {
    auto r0 = p - q;
    auto r1 = mid(p) - q;
    auto r2 = p - mid(q);
    for (int i = 0; i < size_Space; ++i) {
      EXPECT_EQ(r0[i], p[i] - q[i]);
      EXPECT_EQ(r1[i], mid(p[i]) - q[i]);
      EXPECT_EQ(r2[i], p[i] - mid(q[i]));
    }
  }

  virtual void checkAdditiveInverse() {
    auto z = -p;
    for (int i = 0; i < size_Space; ++i)
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
    for (int i = 0; i < size_Space; ++i) {
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
    for (int i = 0; i < size_Space; ++i) {
      EXPECT_EQ(r1[i], p[i] / value);
      EXPECT_EQ(r2[i], p[i] / value_d);
      EXPECT_EQ(r3[i], p[i] / value_i);
      EXPECT_EQ(r4[i], mid(p[i]) / value);
    }
  }

  virtual void checkScalarProduct() {
    IAInterval p0, p1, p2;
    for (int i = 0; i < size_Space; ++i) {
      p0 += p[i] * q[i];
      p1 += mid(p[i]) * q[i];
      p2 += p[i] * mid(q[i]);
    }
    EXPECT_EQ(p0, p * q);
    EXPECT_EQ(p1, mid(p) * q);
    EXPECT_EQ(p2, p * mid(q));
  }

  virtual void checkNorm() {
    IAInterval sum;
    for (int i = 0; i < size_Space; ++i)
      sum += sqr(p[i]);
    EXPECT_EQ(norm(p), sqrt(sum));
  }

  virtual void checkDist() {
    IAInterval d;
    for (int i = 0; i < size_Space; ++i)
      d += sqr(p[i] - q[i]);
    EXPECT_EQ(dist(p, q), sqrt(d));
  }

  virtual void checkCrossProduct() {
    EXPECT_ANY_THROW(p ^ q);
    EXPECT_ANY_THROW(curl(p, q));
  };

  virtual void checkDet() {
    IAPOINT r;
    EXPECT_ANY_THROW(det(p, q, r));
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
    IAPOINT z;
    loader >> z;
    loader.close();
    EXPECT_EQ(p, z);
  }

  virtual void checkMid() {
    auto r = mid(p);
    for (int i = 0; i < size; ++i)
      EXPECT_DOUBLE_EQ(r[i], mid(p[i]));
  }

  virtual void checkInf() {
    auto r = inf(p);
    for (int i = 0; i < size; ++i)
      EXPECT_DOUBLE_EQ(r[i], inf(p[i]));
  }

  virtual void checkSup() {
    auto r = sup(p);
    for (int i = 0; i < size; ++i)
      EXPECT_DOUBLE_EQ(r[i], sup(p[i]));
  }
};

/// To avoid the same code in each test
#define POINT_TESTS(pointTestClass)                                                                \
                                                                                                   \
  TEST_F(pointTestClass, ConstructorTest) { checkConstructor(); }                                  \
                                                                                                   \
  TEST_F(pointTestClass, AssignmentTest) { checkAssignment(); }                                    \
                                                                                                   \
  TEST_F(pointTestClass, SummationTest) { checkSummation(); }                                      \
                                                                                                   \
  TEST_F(pointTestClass, DifferenceTest) { checkDifference(); }                                    \
                                                                                                   \
  TEST_F(pointTestClass, AdditiveInverseTest) { checkAdditiveInverse(); }                          \
                                                                                                   \
  TEST_F(pointTestClass, ScalarMultiplicationTest) { checkScalarMultiplication(); }                \
                                                                                                   \
  TEST_F(pointTestClass, ScalarDivisionTest) { checkScalarDivision(); }                            \
                                                                                                   \
  TEST_F(pointTestClass, ScalarProductTest) { checkScalarProduct(); }                              \
                                                                                                   \
  TEST_F(pointTestClass, NormTest) { checkNorm(); }                                                \
                                                                                                   \
  TEST_F(pointTestClass, DistTest) { checkDist(); }                                                \
                                                                                                   \
  TEST_F(pointTestClass, CrossProductTest) { checkCrossProduct(); }                                \
                                                                                                   \
  TEST_F(pointTestClass, DetTest) { checkDet(); }                                                  \
                                                                                                   \
  TEST_F(pointTestClass, TimeDependenceTest) { EXPECT_EQ(p.isTimeDep(), timeDep); }                \
                                                                                                   \
  TEST_F(pointTestClass, EqualTest) { checkEQ(); }                                                 \
                                                                                                   \
  TEST_F(pointTestClass, NonEqualTest) { checkNonEQ(); }                                           \
                                                                                                   \
  TEST_F(pointTestClass, SaveLoadTest) { checkSaveLoad(); }                                        \
                                                                                                   \
  TEST_F(pointTestClass, MidTest) { checkMid(); }                                                  \
                                                                                                   \
  TEST_F(pointTestClass, InfTest) { checkInf(); }                                                  \
                                                                                                   \
  TEST_F(pointTestClass, SupTest) { checkSup(); }

class IAPointTest_3_1 : public IAPointTest<PointT<IAInterval, 3, 1>> {
protected:
  void checkCrossProduct() override {
    PointT<IAInterval, 3, 1> z;
    z[0] = p[1] * q[2] - p[2] * q[1];
    z[1] = p[2] * q[0] - p[0] * q[2];
    z[2] = p[0] * q[1] - p[1] * q[0];
    EXPECT_EQ(z, p ^ q);
    EXPECT_EQ(z, curl(p, q));
  }

  void checkDet() override {
    PointT<IAInterval, 3, 1> z;
    for (int i = 0; i < 3; ++i)
      z[i] = RandomIAInterval();
    IAInterval d = p[0] * (q[1] * z[2] - q[2] * z[1]) + p[1] * (q[2] * z[0] - q[0] * z[2])
                   + p[2] * (q[0] * z[1] - q[1] * z[0]);
    EXPECT_EQ(det(p, q, z), d);
  }
};

POINT_TESTS(IAPointTest_3_1)

using IAPointTest_2_1 = IAPointTest<PointT<IAInterval, 2, 1>>;

POINT_TESTS(IAPointTest_2_1)

using IAPointTest_1_1 = IAPointTest<PointT<IAInterval, 1, 1>>;

POINT_TESTS(IAPointTest_1_1)

using IAPointTest_4_0 = IAPointTest<PointT<IAInterval, 4, 0>>;

POINT_TESTS(IAPointTest_4_0)

class IAPointTest_3_0 : public IAPointTest<PointT<IAInterval, 3, 0>> {
protected:
  void checkCrossProduct() override {
    PointT<IAInterval, 3, 0> z;
    z[0] = p[1] * q[2] - p[2] * q[1];
    z[1] = p[2] * q[0] - p[0] * q[2];
    z[2] = p[0] * q[1] - p[1] * q[0];
    EXPECT_EQ(z, p ^ q);
    EXPECT_EQ(z, curl(p, q));
  }

  void checkDet() override {
    PointT<IAInterval, 3, 0> z;
    for (int i = 0; i < 3; ++i)
      z[i] = RandomIAInterval();
    IAInterval d = p[0] * (q[1] * z[2] - q[2] * z[1]) + p[1] * (q[2] * z[0] - q[0] * z[2])
                   + p[2] * (q[0] * z[1] - q[1] * z[0]);
    EXPECT_EQ(det(p, q, z), d);
  }
};

POINT_TESTS(IAPointTest_3_0)

using IAPointTest_2_0 = IAPointTest<PointT<IAInterval, 2, 0>>;

POINT_TESTS(IAPointTest_2_0)

using IAPointTest_1_0 = IAPointTest<PointT<IAInterval, 1, 0>>;

POINT_TESTS(IAPointTest_1_0)

int main(int argc, char **argv) {
  MppTest mppTest = MppTestBuilder(argc, argv).WithPPM();
  return mppTest.RUN_ALL_MPP_TESTS();
}
