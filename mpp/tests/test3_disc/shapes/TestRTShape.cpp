#include "TestRTShape.hpp"
#include "Quadrature.hpp"
#include "ShapeTestPoints.hpp"
#include "cells/Cell.hpp"
#include "gtest/gtest.h"
#include "shapes/RTShapes.hpp"

#include "TestEnvironment.hpp"


using namespace ::testing;

const double ansatzTol = 1e-14;
const double derivativeTol = 1e-5;
const double h = 1e-10;

class TestRTShape : public TestWithParam<ShapeTestParameters> {
protected:
  int order;
  const Cell *c = nullptr; // reference, do not delete!
  Shape *shape = nullptr;
  std::vector<std::vector<std::pair<Point, double>>> transformedQuadrature;

  int shapeSize = -1;
  int expectedShapeSize = -1;
  int nodalPointsOnFaces = -1;

  int numNodalPoints = -1;
  int expectedNumNodalPoints = -1;

  vector<Point> nodalPoints = {};
  vector<Point> nodalPointsSorted = {};
  vector<Point> expectedNodalPoints = {};

  TestRTShape(CELLTYPE cellType) : order(GetParam().order) {
    c = ReferenceCell(cellType);
    shape = createRTShape<double, SpaceDimension, TimeDimension>(cellType, order);
    transformedQuadrature.resize(c->Faces());
    const Quadrature &quad = GetQuadrature(FaceCellType(cellType), order + 1);
    for (int face; face < transformedQuadrature.size(); ++face) {
      for (int q = 0; q < quad.size(); ++q) {
        transformedQuadrature[face].push_back(
            std::make_pair(c->FaceLocalToGlobal(face, quad.QPoint(q)),
                           c->LocalFaceArea(face) * quad.Weight(q)));
      }
    }

    shapeSize = shape->size();
    expectedShapeSize = GetParam().shapeSize;

    nodalPointsOnFaces = c->Faces() * (order + 1);
    shape->NodalPoints(*c, nodalPoints);
    numNodalPoints = nodalPoints.size();
    nodalPointsSorted = nodalPoints;
    expectedNodalPoints = GetParam().expectedNodalPoints;

    sort(nodalPointsSorted.begin(), nodalPointsSorted.end());
    sort(expectedNodalPoints.begin(), expectedNodalPoints.end());

    expectedNumNodalPoints = expectedNodalPoints.size();
  }

  ~TestRTShape() override { delete shape; }

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
    for (int i = 0; i < nodalPointsOnFaces; ++i) {
      for (int j = 0; j < nodalPointsOnFaces; ++j) {
        int face = j % c->Faces();
        if (i == j) {
          EXPECT_NEAR(shape->LocalVector(nodalPoints[j], i) * c->LocalFaceNormal(face)
                          * c->LocalFaceArea(face),
                      1.0, ansatzTol);
        } else {
          EXPECT_NEAR(shape->LocalVector(nodalPoints[j], i) * c->LocalFaceNormal(face)
                          * c->LocalFaceArea(face),
                      0.0, ansatzTol);
        }
      }
    }
    for (int i = nodalPointsOnFaces; i < shape->size(); ++i) {
      for (int j = 0; j < nodalPointsOnFaces; ++j) {
        int face = j % c->Faces();
        EXPECT_NEAR(shape->LocalVector(nodalPoints[j], i) * c->LocalFaceNormal(face)
                        * c->LocalFaceArea(face),
                    0.0, ansatzTol);
      }
    }
  }

  void checkDerivatives() {
    std::vector<Point> P = ShapeTestPoints(c->Type());
    for (int n = 0; n < P.size(); ++n) {
      for (int i = 0; i < shapeSize; ++i) {
        double computedDiv = shape->LocalDiv(P[n], i);
        double approxDiv = 0.0;
        for (int k = 0; k < c->dim(); ++k) {
          Point ph = P[n];
          ph[k] += h;
          Point mh = P[n];
          mh[k] -= h;
          approxDiv += (shape->LocalVector(ph, i)[k] - shape->LocalVector(mh, i)[k]) / (2 * h);
        }
        EXPECT_NEAR(approxDiv, computedDiv, derivativeTol);
      }
    }
  }
};

/// To avoid the same code in each test
#define RT_TESTS(testclass, testCases)                                                             \
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
  INSTANTIATE_TEST_CASE_P(TestRTShape, testclass, ValuesIn(testCases));

class IntervalShapeTest : public TestRTShape {
public:
  IntervalShapeTest() : TestRTShape(INTERVAL) {}
};

RT_TESTS(IntervalShapeTest, TestCasesInt)

#if SpaceDimension >= 2

class QuadraticShapeTest : public TestRTShape {
public:
  QuadraticShapeTest() : TestRTShape(QUADRILATERAL) {}
};

RT_TESTS(QuadraticShapeTest, TestCasesQuad)

class TriangularShapeTest : public TestRTShape {
public:
  TriangularShapeTest() : TestRTShape(TRIANGLE) {}
};

RT_TESTS(TriangularShapeTest, TestCasesTri)

#endif
#if SpaceDimension >= 3

class HexahedronShapeTest : public TestRTShape {
public:
  HexahedronShapeTest() : TestRTShape(HEXAHEDRON) {}
};

RT_TESTS(HexahedronShapeTest, TestCasesHex)

// class TetrahedronShapeTest : public TestRTShape {
// public:
//     TetrahedronShapeTest() : TestRTShape(TETRAHEDRON) {}
// };
//
// RT_TESTS(TetrahedronShapeTest, TestCasesTet)

#endif

int main(int argc, char **argv) {
  MppTest mppTest = MppTestBuilder(argc, argv).WithPPM();
  return mppTest.RUN_ALL_MPP_TESTS();
}
