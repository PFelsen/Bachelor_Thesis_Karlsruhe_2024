#include "GaussLobattoShapes.hpp"
#include "TestEnvironment.hpp"

#include "Quadrature.hpp"

using namespace ::testing;

const double ansatzTol = 1e-14;

struct GLShapeParameters {
  CELLTYPE ctype;
  int degree;
  int expectedNodalPointsCount;
};

class TestGaussLobattoShape : public TestWithParam<GLShapeParameters> {
  ShapeT<double, SpaceDimension, TimeDimension> *shape;
  vector<Point> nodalPoints;
  Quadrature quad;
public:
  void SetUp() override {

    shape = createGaussLobattoShape<double, SpaceDimension, TimeDimension>(GetParam().ctype,
                                                                           GetParam().degree);
    const Cell *c = ReferenceCell(GetParam().ctype);
    shape->NodalPoints(*c, nodalPoints);
    quad = GetGLQuadrature(QUADRILATERAL, GetParam().degree + 1);
  }

  void TearDown() override { delete shape; }

  void checkNumNodalPoints() { EXPECT_EQ(nodalPoints.size(), GetParam().expectedNodalPointsCount); }

  void checkAnsatzFunctions() {
    for (int i = 0; i < shape->size(); i++) {
      for (int j = 0; j < shape->size(); j++) {
        if (i == j) {
          EXPECT_NEAR((*shape)(nodalPoints[i], j), 1.0, ansatzTol);
        } else {
          EXPECT_NEAR((*shape)(nodalPoints[i], j), 0.0, ansatzTol);
        }
        // mout << abs((*shape)(nodalPoints[i], j));
      }
      // mout << endl;
    }
  }

  void checkDerivativeFunctions() {
    const auto &NP = nodalPoints;
    for (int i = 0; i < shape->size(); i++) {
      for (int j = 0; j < shape->size(); j++) {
        VectorField grad = shape->LocalGradient(NP[i], j);
        double dx = ((*shape)(NP[i] + Point(1e-10, 0), j) - (*shape)(NP[i], j)) / 1e-10;
        double dy = ((*shape)(NP[i] + Point(0, 1e-10), j) - (*shape)(NP[i], j)) / 1e-10;
        VectorField gradNum = {dx, dy};
        /*std::cout << " ---" << i <<" "<<j << endl;
        std::cout << grad << endl;
        std::cout << gradNum << endl;*/
        double diff = norm(grad - gradNum);
        EXPECT_NEAR(diff, 0, 1e-5);
      }
    }
  }

  void checkShapeFitsToQuadrature() {
    for (int si = 0; si < shape->size(); si++) {
      for (int qi = 0; qi < quad.size(); qi++) {
        EXPECT_NEAR((*shape)(quad.QPoint(qi), si), si == qi ? 1.0 : 0.0, ansatzTol);
      }
    }
  }
};

TEST_P(TestGaussLobattoShape, TestNumNodalPoints) { checkNumNodalPoints(); }

TEST_P(TestGaussLobattoShape, TestAnsatzFunctions) { checkAnsatzFunctions(); }

TEST_P(TestGaussLobattoShape, TestDerivativeFunctions) { checkDerivativeFunctions(); }

TEST_P(TestGaussLobattoShape, TestShapeFitsToQuadrature) { checkShapeFitsToQuadrature(); }

INSTANTIATE_TEST_SUITE_P(
    TestGLQuadShapes, TestGaussLobattoShape,
    Values(
#if SpaceDimension >= 2
        GLShapeParameters{QUADRILATERAL, 1, 4}, GLShapeParameters{QUADRILATERAL, 2, 9},
        GLShapeParameters{QUADRILATERAL, 3, 16}, GLShapeParameters{QUADRILATERAL, 4, 25},
        GLShapeParameters{QUADRILATERAL, 5, 36}, GLShapeParameters{QUADRILATERAL, 6, 49}
#endif
        ));


GTEST_ALLOW_UNINSTANTIATED_PARAMETERIZED_TEST(TestGaussLobattoShape);

int main(int argc, char **argv) {
  MppTest mppTest = MppTestBuilder(argc, argv).WithPPM().WithScreenLogging();
  return mppTest.RUN_ALL_MPP_TESTS();
}