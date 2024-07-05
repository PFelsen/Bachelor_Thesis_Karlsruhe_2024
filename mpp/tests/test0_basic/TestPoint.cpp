#include "Point.hpp"
#include "TestEnvironment.hpp"

#include <vector>

template<typename POINT>
class PointTest : public ::testing::Test {
protected:
  int size = -1;
  int size_Space = -1;
  POINT p, q;
  double value;
  int value_i;

  bool timeDep = false;

  std::vector<double> values;
  std::vector<POINT> points;

  PointTest() {
    mpp_geometry::SetTolerance(1e-15);
    size = p.SpaceDim() + p.TimeDim();
    size_Space = p.SpaceDim();
    for (int i = 0; i < size; ++i) {
      p[i] = RandomDouble();
      q[i] = RandomDouble();
    }
    value = RandomNonZeroDouble();
    value_i = RandomNonZeroInt();
    timeDep = p.TimeDim() == 1;

    POINT tmp;
    values.resize(size);
    points.resize(size);
    for (int i = 0; i < size; ++i)
      values[i] = RandomDouble();
    points[0] = POINT(values[0]);
    if constexpr (tmp.SpaceDim() + tmp.TimeDim() >= 2) points[1] = POINT(values[0], values[1]);
    if constexpr (tmp.SpaceDim() + tmp.TimeDim() >= 3)
      points[2] = POINT(values[0], values[1], values[2]);
    if constexpr (tmp.SpaceDim() + tmp.TimeDim() >= 4)
      points[3] = POINT(values[0], values[1], values[2], values[3]);
  }

  void checkConstructor() {
    POINT z0;
    POINT z1(p);
    for (int i = 0; i < size; ++i) {
      EXPECT_DOUBLE_EQ(z0[i], 0.0);
      EXPECT_DOUBLE_EQ(z1[i], p[i]);
    }
    for (int i = 0; i < size; ++i) {
      for (int k = 0; k < i + 1; ++k)
        EXPECT_DOUBLE_EQ(points[i][k], values[k]);
      for (int k = i + 1; k < size; ++k)
        EXPECT_DOUBLE_EQ(points[i][k], 0.0);
    }
  }

  virtual void checkAssignment() {
    POINT z1 = p;
    POINT z2(p);
    z2 = value;
    POINT z3(p);
    z3 = value_i;
    for (int i = 0; i < size; ++i) {
      EXPECT_DOUBLE_EQ(z1[i], p[i]);
      EXPECT_DOUBLE_EQ(z2[i], value);
      EXPECT_DOUBLE_EQ(z3[i], double(value_i));
    }
  }

  virtual void checkSummation() {
    auto r0 = p + q;
    for (int i = 0; i < size_Space; ++i) {
      EXPECT_DOUBLE_EQ(r0[i], p[i] + q[i]);
    }
  }

  virtual void checkDifference() {
    auto r0 = p - q;
    for (int i = 0; i < size_Space; ++i) {
      EXPECT_DOUBLE_EQ(r0[i], p[i] - q[i]);
    }
  }

  virtual void checkAdditiveInverse() {
    auto z = -p;
    for (int i = 0; i < size_Space; ++i)
      EXPECT_DOUBLE_EQ(z[i], -p[i]);
  }

  virtual void checkScalarMultiplication() {
    auto r1 = p * value;
    auto r2 = value * p;
    auto r3 = p * value_i;
    auto r4 = value_i * p;
    for (int i = 0; i < size_Space; ++i) {
      EXPECT_DOUBLE_EQ(r1[i], p[i] * value);
      EXPECT_DOUBLE_EQ(r2[i], p[i] * value);
      EXPECT_DOUBLE_EQ(r3[i], p[i] * value_i);
      EXPECT_DOUBLE_EQ(r4[i], p[i] * value_i);
    }
  }

  virtual void checkScalarDivision() {
    auto r1 = p / value;
    auto r2 = p / value_i;
    for (int i = 0; i < size_Space; ++i) {
      EXPECT_DOUBLE_EQ(r1[i], p[i] / value);
      EXPECT_DOUBLE_EQ(r2[i], p[i] / value_i);
    }
  }

  virtual void checkScalarProduct() {
    double p0 = 0.0;
    for (int i = 0; i < size_Space; ++i) {
      p0 += p[i] * q[i];
    }
    EXPECT_DOUBLE_EQ(p0, p * q);
  }

  virtual void checkNorm() {
    double sum = 0.0;
    for (int i = 0; i < size_Space; ++i)
      sum += p[i] * p[i];
    EXPECT_DOUBLE_EQ(norm(p), sqrt(sum));
  }

  virtual void checkDist() {
    double d = 0.0;
    for (int i = 0; i < size_Space; ++i)
      d += (p[i] - q[i]) * (p[i] - q[i]);
    EXPECT_DOUBLE_EQ(dist(p, q), sqrt(d));
  }

  virtual void checkCrossProduct() {
    EXPECT_ANY_THROW(p ^ q);
    EXPECT_ANY_THROW(curl(p, q));
  };

  virtual void checkDet() {
    POINT r;
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
    POINT z;
    loader >> z;
    loader.close();
    EXPECT_EQ(p, z);
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
  TEST_F(pointTestClass, SaveLoadTest) { checkSaveLoad(); }

class PointTest_3_1 : public PointTest<PointT<double, 3, 1>> {
protected:
  void checkCrossProduct() override {
    PointT<double, 3, 1> z;
    z[0] = p[1] * q[2] - p[2] * q[1];
    z[1] = p[2] * q[0] - p[0] * q[2];
    z[2] = p[0] * q[1] - p[1] * q[0];
    EXPECT_EQ(z, p ^ q);
    EXPECT_EQ(z, curl(p, q));
  }

  void checkDet() override {
    PointT<double, 3, 1> z;
    for (int i = 0; i < 3; ++i)
      z[i] = RandomDouble();
    double d = p[0] * (q[1] * z[2] - q[2] * z[1]) + p[1] * (q[2] * z[0] - q[0] * z[2])
               + p[2] * (q[0] * z[1] - q[1] * z[0]);
    EXPECT_DOUBLE_EQ(det(p, q, z), d);
  }
};

POINT_TESTS(PointTest_3_1)

using PointTest_2_1 = PointTest<PointT<double, 2, 1>>;

POINT_TESTS(PointTest_2_1)

using PointTest_1_1 = PointTest<PointT<double, 1, 1>>;

POINT_TESTS(PointTest_1_1)

using PointTest_4_0 = PointTest<PointT<double, 4, 0>>;

POINT_TESTS(PointTest_4_0)

class PointTest_3_0 : public PointTest<PointT<double, 3, 0>> {
protected:
  void checkCrossProduct() override {
    PointT<double, 3, 0> z;
    z[0] = p[1] * q[2] - p[2] * q[1];
    z[1] = p[2] * q[0] - p[0] * q[2];
    z[2] = p[0] * q[1] - p[1] * q[0];
    EXPECT_EQ(z, p ^ q);
    EXPECT_EQ(z, curl(p, q));
  }

  void checkDet() override {
    PointT<double, 3, 0> z;
    for (int i = 0; i < 3; ++i)
      z[i] = RandomDouble();
    double d = p[0] * (q[1] * z[2] - q[2] * z[1]) + p[1] * (q[2] * z[0] - q[0] * z[2])
               + p[2] * (q[0] * z[1] - q[1] * z[0]);
    EXPECT_DOUBLE_EQ(det(p, q, z), d);
  }
};

POINT_TESTS(PointTest_3_0)

using PointTest_2_0 = PointTest<PointT<double, 2, 0>>;

POINT_TESTS(PointTest_2_0)

using PointTest_1_0 = PointTest<PointT<double, 1, 0>>;

POINT_TESTS(PointTest_1_0)

int main(int argc, char **argv) {
  MppTest mppTest = MppTestBuilder(argc, argv).WithPPM();
  return mppTest.RUN_ALL_MPP_TESTS();
}
