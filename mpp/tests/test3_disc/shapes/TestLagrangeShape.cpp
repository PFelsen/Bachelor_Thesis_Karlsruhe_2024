#include "TestLagrangeShape.hpp"
#include "Cell.hpp"
#include "LagrangeShapes.hpp"
#include "ShapeTestPoints.hpp"

#include "TestEnvironment.hpp"

using namespace ::testing;

const double ansatzTol = 1e-14;
const double derivativeTol = 1e-5;
const double h = 1e-10;

class TestLagrangeShape : public TestWithParam<ShapeTestParameters> {
protected:
  int degree;
  const Cell *c = nullptr; // reference, do not delete!
  Shape *shape = nullptr;

  int shapeSize = -1;
  int expectedShapeSize = -1;

  int numNodalPoints = -1;
  int expectedNumNodalPoints = -1;

  vector<Point> nodalPoints{};
  vector<Point> nodalPointsSorted{};
  vector<Point> expectedNodalPoints{};

  TestLagrangeShape(CELLTYPE cellType) : degree(GetParam().degree) {
    c = ReferenceCell(cellType);
    shape = createLagrangeShape<double, SpaceDimension, TimeDimension>(cellType, degree);

    shapeSize = shape->size();
    expectedShapeSize = GetParam().expectedNodalPoints.size();

    numNodalPoints = shape->size();
    shape->NodalPoints(*c, nodalPoints);
    nodalPointsSorted = nodalPoints;
    expectedNodalPoints = GetParam().expectedNodalPoints;

    sort(nodalPointsSorted.begin(), nodalPointsSorted.end());
    sort(expectedNodalPoints.begin(), expectedNodalPoints.end());

    expectedNumNodalPoints = expectedNodalPoints.size();
  }

  ~TestLagrangeShape() override { delete shape; }

  void checkShapeSize() { EXPECT_EQ(shapeSize, expectedShapeSize); }

  void checkNumNodalPoints() { EXPECT_EQ(numNodalPoints, expectedNumNodalPoints); }

  void checkNodalPoints() {
    for (int i = 0; i < numNodalPoints; i++) {
      for (int d = 0; d < nodalPointsSorted[i].SpaceDim(); d++) {
        EXPECT_DOUBLE_EQ(expectedNodalPoints[i][d], nodalPointsSorted[i][d]);
      }
    }
  }

  void checkFunctions() {
    for (int i = 0; i < shapeSize; i++) {
      for (int j = 0; j < shapeSize; j++) {
        if (i == j) {
          EXPECT_NEAR((*shape)(nodalPoints[i], j), 1.0, ansatzTol);
        } else {
          EXPECT_NEAR((*shape)(nodalPoints[i], j), 0.0, ansatzTol);
        }
      }
    }
  }

  void checkDerivatives() {
    std::vector<Point> P = ShapeTestPoints(c->Type());
    for (int n = 0; n < P.size(); ++n) {
      for (int i = 0; i < shapeSize; ++i) {
        VectorField D = shape->LocalGradient(P[n], i);
        for (int k = 0; k < c->dim(); ++k) {
          Point ph = P[n];
          ph[k] += h;
          Point mh = P[n];
          mh[k] -= h;
          EXPECT_NEAR(((*shape)(ph, i) - (*shape)(mh, i)) / (2 * h), D[k], derivativeTol);
        }
      }
    }
  }
};

/// To avoid the same code in each test
#define LAGRANGE_TESTS(testclass, testCases)                                                       \
                                                                                                   \
  TEST_P(testclass, TestNumNodalPoints) { checkNumNodalPoints(); }                                 \
                                                                                                   \
  TEST_P(testclass, TestShapeSize) { checkShapeSize(); }                                           \
                                                                                                   \
  TEST_P(testclass, TestNodalPoints) { checkNodalPoints(); }                                       \
                                                                                                   \
  TEST_P(testclass, TestFunctions) { checkFunctions(); }                                           \
                                                                                                   \
  TEST_P(testclass, TestDerivatives) { checkDerivatives(); }                                       \
                                                                                                   \
  INSTANTIATE_TEST_CASE_P(TestLagrangeShape, testclass, ValuesIn(testCases));

class IntervalShapeTest : public TestLagrangeShape {
public:
  IntervalShapeTest() : TestLagrangeShape(INTERVAL) {}
};

LAGRANGE_TESTS(IntervalShapeTest, TestCasesInt)

#if SpaceDimension >= 2

class QuadraticShapeTest : public TestLagrangeShape {
public:
  QuadraticShapeTest() : TestLagrangeShape(QUADRILATERAL) {}
};

LAGRANGE_TESTS(QuadraticShapeTest, TestCasesQuad)

class TriangularShapeTest : public TestLagrangeShape {
public:
  TriangularShapeTest() : TestLagrangeShape(TRIANGLE) {}
};

LAGRANGE_TESTS(TriangularShapeTest, TestCasesTri)

#endif
#if SpaceDimension >= 3

class HexahedronShapeTest : public TestLagrangeShape {
public:
  HexahedronShapeTest() : TestLagrangeShape(HEXAHEDRON) {}
};

LAGRANGE_TESTS(HexahedronShapeTest, TestCasesHex)

class TetrahedronShapeTest : public TestLagrangeShape {
public:
  TetrahedronShapeTest() : TestLagrangeShape(TETRAHEDRON) {}
};

LAGRANGE_TESTS(TetrahedronShapeTest, TestCasesTet)

#endif

int main(int argc, char **argv) {
  MppTest mppTest = MppTestBuilder(argc, argv).WithPPM();
  return mppTest.RUN_ALL_MPP_TESTS();
}
