#include "ArgyrisShapes.hpp"
#include "Cell.hpp"
#include "ShapeTestPoints.hpp"
#include "TestEnvironment.hpp"

const double ansatzTol = 1e-14;
const double derivativeTol = 1e-5;
const double h = 1e-10;

class TestArgyrisShape : public ::testing::Test {
protected:
  const Cell *c = nullptr; // reference, do not delete!
  Shape *shape = nullptr;

  int shapeSize = -1;
  int expectedShapeSize = 21;
  int numNodalPoints = -1;
  int expectedNumNodalPoints = -1;


  vector<Point> nodalPoints = {};
  vector<Point> expectedNodalPoints = {Point(0.0, 0.0), Point(1.0, 0.0), Point(0.0, 1.0),
                                       Point(0.5, 0.0), Point(0.5, 0.5), Point(0.0, 0.5)};

  TestArgyrisShape() {
    c = ReferenceCell(TRIANGLE);
    shape = new ArgyrisTri();

    shapeSize = shape->size();
    shape->NodalPoints(*c, nodalPoints);

    expectedNumNodalPoints = expectedNodalPoints.size();
    numNodalPoints = nodalPoints.size();
  }

  ~TestArgyrisShape() { delete shape; }

  void fillNodalPointData(vector<double> &values_i, int i, int nodalPoint) {
    values_i[nodalPoint * 6] = (*shape)(expectedNodalPoints[nodalPoint], i);
    VectorField g = shape->LocalGradient(expectedNodalPoints[nodalPoint], i);
    values_i[nodalPoint * 6 + 1] = g[0];
    values_i[nodalPoint * 6 + 2] = g[1];
    SymTensor h = shape->LocalHessian(expectedNodalPoints[nodalPoint], i);
    values_i[nodalPoint * 6 + 3] = h(0, 0);
    values_i[nodalPoint * 6 + 4] = h(0, 1);
    values_i[nodalPoint * 6 + 5] = h(1, 1);
  }
};
#if SpaceDimension >= 2

TEST_F(TestArgyrisShape, TestNumNodalPoints) { EXPECT_EQ(numNodalPoints, expectedNumNodalPoints); }

TEST_F(TestArgyrisShape, TestShapeSize) { EXPECT_EQ(shapeSize, expectedShapeSize); }

TEST_F(TestArgyrisShape, TestNodalPoints) {
  for (int i = 0; i < numNodalPoints; i++) {
    for (int d = 0; d < nodalPoints[i].SpaceDim(); d++) {
      EXPECT_DOUBLE_EQ(expectedNodalPoints[i][d], nodalPoints[i][d]);
    }
  }
}

TEST_F(TestArgyrisShape, TestFunctions) {
  vector<vector<double>> values(21);
  for (int i = 0; i < 21; ++i) {
    values[i].resize(21);
    for (int np = 0; np < 3; ++np)
      fillNodalPointData(values[i], i, np);
    VectorField g_0 = shape->LocalGradient(expectedNodalPoints[3], i);
    values[i][18] = -g_0[1];
    VectorField g_1 = shape->LocalGradient(expectedNodalPoints[4], i);
    values[i][19] = (g_1[0] + g_1[1]) / sqrt(2.0);
    VectorField g_2 = shape->LocalGradient(expectedNodalPoints[5], i);
    values[i][20] = -g_2[0];
  }

  for (int i = 0; i < 21; ++i) {
    for (int j = 0; j < 21; ++j)
      if (i == j) {
        EXPECT_NEAR(values[i][j], 1.0, ansatzTol);
      } else {
        EXPECT_NEAR(values[i][j], 0.0, ansatzTol);
      }
  }
}

TEST_F(TestArgyrisShape, TestDerivatives) {
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

TEST_F(TestArgyrisShape, TestHessian) {
  std::vector<Point> P = ShapeTestPoints(c->Type());
  for (int n = 0; n < P.size(); ++n) {
    for (int i = 0; i < shapeSize; ++i) {
      SymTensor H = shape->LocalHessian(P[n], i);
      for (int k = 0; k < c->dim(); ++k) {
        for (int l = 0; l < c->dim(); ++l) {
          Point ph = P[n];
          ph[l] += h;
          Point mh = P[n];
          mh[l] -= h;
          EXPECT_NEAR((shape->LocalGradient(ph, i)[k] - shape->LocalGradient(mh, i)[k]) / (2 * h),
                      H(k, l), derivativeTol);
        }
      }
    }
  }
}
#endif

int main(int argc, char **argv) {
  MppTest mppTest = MppTestBuilder(argc, argv).WithPPM();
  return mppTest.RUN_ALL_MPP_TESTS();
}