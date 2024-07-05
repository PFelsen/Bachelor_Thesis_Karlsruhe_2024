#include <vector>
#include "Quadrature.hpp"
#include "TestEnvironment.hpp"


using namespace ::testing;
using std::string;
using std::vector;

const double valueTol = 1e-14;

int factorial(int n) { return (n == 1 || n == 0) ? 1 : factorial(n - 1) * n; }

struct QuadratureTestParameter {
  CELLTYPE cellType;
  int exactUpTo;
  int numQuadraturePoints;
};

class QuadratureTest : public TestWithParam<QuadratureTestParameter> {
protected:
  const Quadrature &quad;
  int exactUpTo;

  int expectedNumQuadraturePoints;

  vector<double> quadratureValues;
  vector<double> expectedQuadratureValues;
public:
  QuadratureTest() : quad(GetQuadrature(GetParam().cellType, GetParam().exactUpTo)) {
    exactUpTo = GetParam().exactUpTo;
    expectedNumQuadraturePoints = GetParam().numQuadraturePoints;
  }

  int numValues() { return quadratureValues.size(); }

  double evalMonom(const Point &z, vector<int> exponent) {
    double value = pow(z[0], exponent[0]);
    for (int d = 1; d < exponent.size(); ++d)
      value *= pow(z[d], exponent[d]);
    return value;
  }

  double evalQuadrature(vector<int> exponent) {
    double value = 0.0;
    for (int q = 0; q < quad.size(); ++q)
      value += quad.Weight(q) * evalMonom(quad.QPoint(q), exponent);
    return value;
  }
};

/*
 * Tests for quadrature rules on intervals
 */
class IntQuadratureTest : public QuadratureTest {
public:
  IntQuadratureTest() {
    for (int d = 0; d <= exactUpTo; ++d) {
      vector<int> exponent(1);
      exponent[0] = d;
      quadratureValues.push_back(evalQuadrature(exponent));
      expectedQuadratureValues.push_back(1.0 / (d + 1));
    }
  }
};

TEST_P(IntQuadratureTest, NumQuadraturePointsTest) {
  EXPECT_EQ(quad.size(), expectedNumQuadraturePoints);
}

TEST_P(IntQuadratureTest, QuadratureValueTest) {
  for (int n = 0; n < numValues(); ++n)
    EXPECT_NEAR(quadratureValues[n], expectedQuadratureValues[n], valueTol);
}

INSTANTIATE_TEST_SUITE_P(
    QuadratureTest, IntQuadratureTest,
    Values(QuadratureTestParameter{INTERVAL, 0, 1}, QuadratureTestParameter{INTERVAL, 1, 1},
           QuadratureTestParameter{INTERVAL, 2, 2}, QuadratureTestParameter{INTERVAL, 3, 2},
           QuadratureTestParameter{INTERVAL, 4, 3}, QuadratureTestParameter{INTERVAL, 5, 3},
           QuadratureTestParameter{INTERVAL, 6, 4}, QuadratureTestParameter{INTERVAL, 7, 4},
           QuadratureTestParameter{INTERVAL, 8, 5}, QuadratureTestParameter{INTERVAL, 9, 5},
           QuadratureTestParameter{INTERVAL, 10, 6}, QuadratureTestParameter{INTERVAL, 11, 6},
           QuadratureTestParameter{INTERVAL, 12, 7}, QuadratureTestParameter{INTERVAL, 13, 7},
           QuadratureTestParameter{INTERVAL, 14, 8}, QuadratureTestParameter{INTERVAL, 15, 8},
           QuadratureTestParameter{INTERVAL, 16, 9}, QuadratureTestParameter{INTERVAL, 17, 9},
           QuadratureTestParameter{INTERVAL, 18, 10}, QuadratureTestParameter{INTERVAL, 19, 10}));

#if SpaceDimension >= 2

/*
 * Tests for quadrature rules on triangles
 */
class TriQuadratureTest : public QuadratureTest {
public:
  TriQuadratureTest() {
    for (int d = 0; d <= exactUpTo; ++d)
      for (int p = 0; p <= d; ++p) {
        vector<int> exponent(2);
        exponent[0] = p;
        exponent[1] = d - p;
        quadratureValues.push_back(evalQuadrature(exponent));
        expectedQuadratureValues.push_back(expectedValue(p, d - p));
      }
  }
protected:
  // int_T x^p y^q d(x,y)= p!*q!/(p+q+2)!
  double expectedValue(int p, int q) {
    double value = 1.0;
    if (p <= q)
      for (int i = 1; i <= q; ++i)
        value *= double(i) / (p + i);
    else
      for (int i = 1; i <= p; ++i)
        value *= double(i) / (q + i);
    return value * 1.0 / (p + q + 1) * 1.0 / (p + q + 2);
  }
};

TEST_P(TriQuadratureTest, NumQuadraturePointsTest) {
  EXPECT_EQ(quad.size(), expectedNumQuadraturePoints);
}

TEST_P(TriQuadratureTest, QuadratureValueTest) {
  for (int n = 0; n < numValues(); ++n)
    EXPECT_NEAR(quadratureValues[n], expectedQuadratureValues[n], valueTol);
}

INSTANTIATE_TEST_SUITE_P(
    QuadratureTest, TriQuadratureTest,
    Values(QuadratureTestParameter{TRIANGLE, 0, 1}, QuadratureTestParameter{TRIANGLE, 1, 1},
           QuadratureTestParameter{TRIANGLE, 2, 3}, QuadratureTestParameter{TRIANGLE, 3, 4},
           QuadratureTestParameter{TRIANGLE, 4, 6}, QuadratureTestParameter{TRIANGLE, 5, 7},
           QuadratureTestParameter{TRIANGLE, 6, 12}, QuadratureTestParameter{TRIANGLE, 7, 13},
           QuadratureTestParameter{TRIANGLE, 8, 16}, QuadratureTestParameter{TRIANGLE, 9, 19},
           QuadratureTestParameter{TRIANGLE, 10, 25}, QuadratureTestParameter{TRIANGLE, 11, 28},
           QuadratureTestParameter{TRIANGLE, 12, 33}, QuadratureTestParameter{TRIANGLE, 13, 37},
           QuadratureTestParameter{TRIANGLE, 14, 42}, QuadratureTestParameter{TRIANGLE, 15, 49},
           QuadratureTestParameter{TRIANGLE, 16, 55}, QuadratureTestParameter{TRIANGLE, 17, 60},
           QuadratureTestParameter{TRIANGLE, 18, 67}, QuadratureTestParameter{TRIANGLE, 19, 73},
           QuadratureTestParameter{TRIANGLE, 20, 79}));

/*
 * Tests for quadrature rules on quadrilaterals
 */
class QuadQuadratureTest : public QuadratureTest {
public:
  QuadQuadratureTest() {
    for (int d = 0; d <= exactUpTo; ++d)
      for (int p = 0; p <= d; ++p) {
        vector<int> exponent(2);
        exponent[0] = p;
        exponent[1] = d - p;
        quadratureValues.push_back(evalQuadrature(exponent));
        expectedQuadratureValues.push_back(1.0 / (p + 1) * 1.0 / (d - p + 1));
      }
  }
};

TEST_P(QuadQuadratureTest, NumQuadraturePointsTest) {
  EXPECT_EQ(quad.size(), expectedNumQuadraturePoints);
}

TEST_P(QuadQuadratureTest, QuadratureValueTest) {
  for (int n = 0; n < numValues(); ++n)
    EXPECT_NEAR(quadratureValues[n], expectedQuadratureValues[n], valueTol);
}

INSTANTIATE_TEST_SUITE_P(QuadratureTest, QuadQuadratureTest,
                         Values(QuadratureTestParameter{QUADRILATERAL, 0, 1},
                                QuadratureTestParameter{QUADRILATERAL, 1, 1},
                                QuadratureTestParameter{QUADRILATERAL, 2, 4},
                                QuadratureTestParameter{QUADRILATERAL, 3, 4},
                                QuadratureTestParameter{QUADRILATERAL, 4, 9},
                                QuadratureTestParameter{QUADRILATERAL, 5, 9},
                                QuadratureTestParameter{QUADRILATERAL, 6, 16},
                                QuadratureTestParameter{QUADRILATERAL, 7, 16},
                                QuadratureTestParameter{QUADRILATERAL, 8, 25},
                                QuadratureTestParameter{QUADRILATERAL, 9, 25},
                                QuadratureTestParameter{QUADRILATERAL, 10, 36},
                                QuadratureTestParameter{QUADRILATERAL, 11, 36},
                                QuadratureTestParameter{QUADRILATERAL, 12, 49},
                                QuadratureTestParameter{QUADRILATERAL, 13, 49},
                                QuadratureTestParameter{QUADRILATERAL, 14, 64},
                                QuadratureTestParameter{QUADRILATERAL, 15, 64},
                                QuadratureTestParameter{QUADRILATERAL, 16, 81},
                                QuadratureTestParameter{QUADRILATERAL, 17, 81},
                                QuadratureTestParameter{QUADRILATERAL, 18, 100},
                                QuadratureTestParameter{QUADRILATERAL, 19, 100}));

#endif
#if SpaceDimension >= 3

/*
 * Tests for quadrature rules on tetrahedrons
 */
class TetQuadratureTest : public QuadratureTest {
public:
  TetQuadratureTest() {
    for (int d = 0; d <= exactUpTo; ++d)
      for (int p = 0; p <= d; ++p)
        for (int r = 0; r <= d - p; ++r) {
          vector<int> exponent(3);
          exponent[0] = p;
          exponent[1] = r;
          exponent[2] = d - r - p;
          quadratureValues.push_back(evalQuadrature(exponent));
          expectedQuadratureValues.push_back(expectedValue(p, r, d - r - p));
        }
  }
protected:
  // int_T x^p y^q z^r d(x,y,z)= p!*q!*r!/(p+q+r+3)!
  double expectedValue(int p, int q, int r) {
    double value = 1.0;
    vector<int> factorials{p, q, r};
    std::sort(factorials.rbegin(), factorials.rend());

    for (int s = 1; s <= factorials[1]; ++s) {
      value *= double(s) / double(factorials[0] + s);
    }
    for (int t = 1; t <= factorials[2]; ++t) {
      value *= double(t) / double(factorials[0] + factorials[1] + t);
    }
    return value * 1.0 / (p + q + r + 1) * 1.0 / (p + q + r + 2) * 1.0 / (p + q + r + 3);
  }
};

TEST_P(TetQuadratureTest, NumQuadraturePointsTest) {
  EXPECT_EQ(quad.size(), expectedNumQuadraturePoints);
}

TEST_P(TetQuadratureTest, QuadratureValueTest) {
  for (int n = 0; n < numValues(); ++n)
    EXPECT_NEAR(quadratureValues[n], expectedQuadratureValues[n], valueTol);
}

INSTANTIATE_TEST_SUITE_P(
    QuadratureTest, TetQuadratureTest,
    Values(QuadratureTestParameter{TETRAHEDRON, 0, 1}, QuadratureTestParameter{TETRAHEDRON, 1, 1},
           QuadratureTestParameter{TETRAHEDRON, 2, 4}, QuadratureTestParameter{TETRAHEDRON, 3, 5},
           QuadratureTestParameter{TETRAHEDRON, 4, 11}, QuadratureTestParameter{TETRAHEDRON, 5, 14},
           QuadratureTestParameter{TETRAHEDRON, 6, 24}, QuadratureTestParameter{TETRAHEDRON, 7, 35},
           QuadratureTestParameter{TETRAHEDRON, 8, 46}, QuadratureTestParameter{TETRAHEDRON, 9, 59},
           QuadratureTestParameter{TETRAHEDRON, 10, 79},
           QuadratureTestParameter{TETRAHEDRON, 11, 96},
           QuadratureTestParameter{TETRAHEDRON, 12, 127},
           QuadratureTestParameter{TETRAHEDRON, 13, 149},
           QuadratureTestParameter{TETRAHEDRON, 14, 194},
           QuadratureTestParameter{TETRAHEDRON, 15, 246}));

/*
 * Tests for quadrature rules on prisms
 */
class PriQuadratureTest : public QuadratureTest {
public:
  PriQuadratureTest() {
    for (int d = 0; d <= exactUpTo; ++d)
      for (int p = 0; p <= d; ++p)
        for (int r = 0; r <= d - p; ++r) {
          vector<int> exponent(3);
          exponent[0] = p;
          exponent[1] = r;
          exponent[2] = d - r - p;
          quadratureValues.push_back(evalQuadrature(exponent));
          expectedQuadratureValues.push_back(expectedValue(p, r, d - r - p));
        }
  }
protected:
  // int_T x^p y^q d(x,y)= p!*q!/(p+q+2)!
  double expectedValue(int p, int q, int r) {
    double value = 1.0;
    if (p <= q)
      for (int i = 1; i <= q; ++i)
        value *= double(i) / (p + i);
    else
      for (int i = 1; i <= p; ++i)
        value *= double(i) / (q + i);
    return value * 1.0 / (p + q + 1) * 1.0 / (p + q + 2) * 1.0 / (r + 1);
  }
};

TEST_P(PriQuadratureTest, NumQuadraturePointsTest) {
  EXPECT_EQ(quad.size(), expectedNumQuadraturePoints);
}

TEST_P(PriQuadratureTest, QuadratureValueTest) {
  for (int n = 0; n < numValues(); ++n)
    EXPECT_NEAR(quadratureValues[n], expectedQuadratureValues[n], valueTol);
}

INSTANTIATE_TEST_SUITE_P(
    QuadratureTest, PriQuadratureTest,
    Values(QuadratureTestParameter{PRISM, 0, 1}, QuadratureTestParameter{PRISM, 1, 1},
           QuadratureTestParameter{PRISM, 2, 6}, QuadratureTestParameter{PRISM, 3, 8},
           QuadratureTestParameter{PRISM, 4, 18}, QuadratureTestParameter{PRISM, 5, 21},
           QuadratureTestParameter{PRISM, 6, 48}, QuadratureTestParameter{PRISM, 7, 52},
           QuadratureTestParameter{PRISM, 8, 80}, QuadratureTestParameter{PRISM, 9, 95},
           QuadratureTestParameter{PRISM, 10, 150}, QuadratureTestParameter{PRISM, 11, 168},
           QuadratureTestParameter{PRISM, 12, 231}, QuadratureTestParameter{PRISM, 13, 259},
           QuadratureTestParameter{PRISM, 14, 336}, QuadratureTestParameter{PRISM, 15, 392},
           QuadratureTestParameter{PRISM, 16, 495}, QuadratureTestParameter{PRISM, 17, 540},
           QuadratureTestParameter{PRISM, 18, 670}, QuadratureTestParameter{PRISM, 19, 730}));

/*
 * Tests for quadrature rules on hexahedrons
 */
class HexQuadratureTest : public QuadratureTest {
public:
  HexQuadratureTest() {
    for (int d = 0; d <= exactUpTo; ++d)
      for (int p = 0; p <= d; ++p)
        for (int r = 0; r <= d - p; ++r) {
          vector<int> exponent(3);
          exponent[0] = p;
          exponent[1] = r;
          exponent[2] = d - r - p;
          quadratureValues.push_back(evalQuadrature(exponent));
          expectedQuadratureValues.push_back(1.0 / (p + 1) * 1.0 / (r + 1) * 1.0 / (d - p - r + 1));
        }
  }
};

TEST_P(HexQuadratureTest, NumQuadraturePointsTest) {
  EXPECT_EQ(quad.size(), expectedNumQuadraturePoints);
}

TEST_P(HexQuadratureTest, QuadratureValueTest) {
  for (int n = 0; n < numValues(); ++n)
    EXPECT_NEAR(quadratureValues[n], expectedQuadratureValues[n], valueTol);
}

INSTANTIATE_TEST_SUITE_P(
    QuadratureTest, HexQuadratureTest,
    Values(
        QuadratureTestParameter{HEXAHEDRON, 0, 1}, QuadratureTestParameter{HEXAHEDRON, 1, 1},
        QuadratureTestParameter{HEXAHEDRON, 2, 8}, QuadratureTestParameter{HEXAHEDRON, 3, 8},
        QuadratureTestParameter{HEXAHEDRON, 4, 27}, QuadratureTestParameter{HEXAHEDRON, 5, 27},
        QuadratureTestParameter{HEXAHEDRON, 6, 64}, QuadratureTestParameter{HEXAHEDRON, 7, 64},
        QuadratureTestParameter{HEXAHEDRON, 8, 125}, QuadratureTestParameter{HEXAHEDRON, 9, 125},
        QuadratureTestParameter{HEXAHEDRON, 10, 216}, QuadratureTestParameter{HEXAHEDRON, 11, 216},
        QuadratureTestParameter{HEXAHEDRON, 12, 343}, QuadratureTestParameter{HEXAHEDRON, 13, 343},
        QuadratureTestParameter{HEXAHEDRON, 14, 512}, QuadratureTestParameter{HEXAHEDRON, 15, 512},
        QuadratureTestParameter{HEXAHEDRON, 16, 729}, QuadratureTestParameter{HEXAHEDRON, 17, 729},
        QuadratureTestParameter{HEXAHEDRON, 18, 1000},
        QuadratureTestParameter{HEXAHEDRON, 19, 1000}));

#endif

/*
 * Tests for old getter with string parameter
 */
struct StringGetterTestParameter {
  string name;
  int numQuadraturePoints;
};

class StringGetterTest : public TestWithParam<StringGetterTestParameter> {
protected:
  const Quadrature &quad;
  int expectedNumQuadraturePoints;
public:
  StringGetterTest() : quad(GetQuadrature(GetParam().name)) {
    expectedNumQuadraturePoints = GetParam().numQuadraturePoints;
  }
};

TEST_P(StringGetterTest, NumQuadraturePointsTest) {
  EXPECT_EQ(quad.size(), expectedNumQuadraturePoints);
}

INSTANTIATE_TEST_SUITE_P(
    StringGetterTest, StringGetterTest,
    Values(StringGetterTestParameter{"Qint1", 1}, StringGetterTestParameter{"Qint2", 2},
           StringGetterTestParameter{"Qint3", 3}, StringGetterTestParameter{"Qint4", 4},
           StringGetterTestParameter{"Qint5", 5}, StringGetterTestParameter{"Qint6", 6},
           StringGetterTestParameter{"Qint7", 7}, StringGetterTestParameter{"Qint8", 8},
           StringGetterTestParameter{"Qint9", 9}
#if SpaceDimension >= 2
           ,
           StringGetterTestParameter{"Qtri1", 1}, StringGetterTestParameter{"Qtri3", 3},
           StringGetterTestParameter{"Qtri4", 4}, StringGetterTestParameter{"Qtri6", 6},
           StringGetterTestParameter{"Qtri7", 7}, StringGetterTestParameter{"Qtri12", 12},
           StringGetterTestParameter{"Qtri13", 13}, StringGetterTestParameter{"Qtri16", 16},
           StringGetterTestParameter{"Qtri25", 25}, StringGetterTestParameter{"Qtri36", 33},
           StringGetterTestParameter{"Qquad1", 1}, StringGetterTestParameter{"Qquad4", 4},
           StringGetterTestParameter{"Qquad9", 9}, StringGetterTestParameter{"Qquad16", 16},
           StringGetterTestParameter{"Qquad25", 25}, StringGetterTestParameter{"Qquad36", 36},
           StringGetterTestParameter{"Qquad49", 49}, StringGetterTestParameter{"Qquad64", 64},
           StringGetterTestParameter{"Qquad81", 81}
#endif
#if SpaceDimension >= 3
           ,
           StringGetterTestParameter{"Qtet1", 1}, StringGetterTestParameter{"Qtet4", 4},
           StringGetterTestParameter{"Qtet11", 11}, StringGetterTestParameter{"Qtet15", 14},
           StringGetterTestParameter{"Qpri6", 6}, StringGetterTestParameter{"Qpri8", 8},
           StringGetterTestParameter{"Qhex1", 1}, StringGetterTestParameter{"Qhex8", 8},
           StringGetterTestParameter{"Qhex27", 27}, StringGetterTestParameter{"Qhex64", 64},
           StringGetterTestParameter{"Qhex125", 125}, StringGetterTestParameter{"Qhex216", 216},
           StringGetterTestParameter{"Qhex343", 343}, StringGetterTestParameter{"Qhex512", 512},
           StringGetterTestParameter{"Qhex729", 729}
#endif
           ));

int main(int argc, char **argv) {
  MppTest mppTest = MppTestBuilder(argc, argv).WithoutDefaultConfig();
  return mppTest.RUN_ALL_MPP_TESTS();
}
