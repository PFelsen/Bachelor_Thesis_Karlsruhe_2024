#include <Quadrature.hpp>
#include "GaussLobattoShapes.hpp"
#include "LagrangeShapes.hpp"
#include "SpaceTimeQuadrature.hpp"
#include "SpaceTimeShape.hpp"
#include "TCell.hpp"
#include "TestEnvironment.hpp"
#include "gtest/gtest.h"

const double ansatzTol = 1e-14;

enum SHAPETYPE { LAGRANGE, GAUSSLOBATTO };

struct ShapeTestParameters {
  SHAPETYPE stype;
  int pdeg;
  CELLTYPE ctype;
  std::vector<Point> expectedNodalPoints;
  const Cell *cell;
};

class TestSpaceTimeShape : public TestWithParam<ShapeTestParameters> {
protected:
  SHAPETYPE stype;
  int pdeg;
  std::vector<Point> nodalPoints = {};
  std::vector<Point> nodalPointsSorted = {};
  std::vector<Point> expectedNodalPoints = {};

  Shape *sc;
  Shape *si;
  std::unique_ptr<SpaceTimeShape> shape;
  std::unique_ptr<SpaceTimeQuadrature> quad;
public:
  TestSpaceTimeShape() : stype(GetParam().stype), pdeg(GetParam().pdeg) {
    CELLTYPE celltype = GetParam().ctype;
    CELLTYPE spaceCelltype = SpaceCellType(celltype);
    Cell *cell = new TCell(ReferenceCell(spaceCelltype), TimeInterval(0, 1));

    quad = CreateSpaceTimeQuad(spaceCelltype);
    shape = CreateSpaceTimeShape(spaceCelltype);
    shape->NodalPoints(*cell, nodalPoints);
    nodalPointsSorted = nodalPoints;
    expectedNodalPoints = GetParam().expectedNodalPoints;
  }

  std::unique_ptr<SpaceTimeShape> CreateSpaceTimeShape(CELLTYPE spaceCelltype) {
    if (stype == LAGRANGE) {
      sc = createLagrangeShape<double, SpaceDimension, TimeDimension>(spaceCelltype, pdeg);
      si = createLagrangeShape<double, SpaceDimension, TimeDimension>(INTERVAL, pdeg);
    } else if (stype == GAUSSLOBATTO) {
      sc = createGaussLobattoShape<double, SpaceDimension, TimeDimension>(spaceCelltype, pdeg);
      si = createGaussLobattoShape<double, SpaceDimension, TimeDimension>(INTERVAL, pdeg);
    }
    return std::make_unique<SpaceTimeShape>(*sc, *si);
  }

  Quadrature CreateOneDQuad() {
    if (stype == LAGRANGE) {
      return Qint(2 * pdeg);
    } else if (stype == GAUSSLOBATTO) {
      return QintGaussLobatto(2 * pdeg);
    }
    Exit("no quad found")
  }

  std::unique_ptr<SpaceTimeQuadrature> CreateSpaceTimeQuad(CELLTYPE spaceCelltype) {
    Quadrature oned = CreateOneDQuad();
    if (spaceCelltype == INTERVAL) {
      return std::make_unique<SpaceTimeQuadrature>(oned, oned);
    } else if (spaceCelltype == QUADRILATERAL) {
      QTensor SQ = QTensor({oned, oned});
      return std::make_unique<SpaceTimeQuadrature>(SQ, oned);
    }
    Exit("no quad found for " + std::to_string(spaceCelltype));
  }

  ~TestSpaceTimeShape() override {}

  void checkNumNodalPoints() const {
    EXPECT_EQ(nodalPoints.size(), expectedNodalPoints.size());
    EXPECT_EQ(nodalPoints.size(), shape->size());
  }

  void checkNodalPoints() const {
    int tt = nodalPointsSorted[0].SpaceDim() + nodalPointsSorted[0].TimeDim();
    for (int i = 0; i < nodalPoints.size(); i++) {
      for (int d = 0; d < tt; d++) {
        if (stype == GAUSSLOBATTO) {
          EXPECT_DOUBLE_EQ(expectedNodalPoints[i][d], quad->QPoint(i)[d]);
        }
        EXPECT_DOUBLE_EQ(expectedNodalPoints[i][d], nodalPointsSorted[i][d]);
      }
    }
  }

  float delta(int i, int j) const { return (i == j) ? 1.0 : 0.0; }

  void checkAnsatzFunctions() const {
    for (int i = 0; i < nodalPoints.size(); i++) {
      for (int j = 0; j < nodalPoints.size(); j++) {
        if (stype == GAUSSLOBATTO) {
          EXPECT_NEAR((*shape)(quad->QPoint(i), j), delta(i, j), ansatzTol);
        }
        EXPECT_NEAR((*shape)(nodalPoints[i], j), delta(i, j), ansatzTol);
      }
    }
  }
};

TEST_P(TestSpaceTimeShape, TestNumNodalPoints) { checkNumNodalPoints(); }

TEST_P(TestSpaceTimeShape, TestNodalPoints) { checkNodalPoints(); }

TEST_P(TestSpaceTimeShape, TestAnsatzFunctions) { checkAnsatzFunctions(); }

INSTANTIATE_TEST_CASE_P(
    TestSpaceTimeShape, TestSpaceTimeShape,
    Values(
        ShapeTestParameters{GAUSSLOBATTO, 1, SPACETIME_INTERVAL,
                            std::vector<Point>{Point(0, 0, 0).WithT(0), Point{1, 0, 0}.WithT(0),
                                               Point{0, 0, 0}.WithT(1), Point{1, 0, 0}.WithT(1)}},
        ShapeTestParameters{LAGRANGE, 1, SPACETIME_INTERVAL,
                            std::vector<Point>{Point(0, 0, 0).WithT(0), Point{1, 0, 0}.WithT(0),
                                               Point{0, 0, 0}.WithT(1), Point{1, 0, 0}.WithT(1)}},
        ShapeTestParameters{GAUSSLOBATTO, 1, SPACETIME_QUADRILATERAL,
                            std::vector<Point>{Point(0, 0, 0).WithT(0), Point{1, 0, 0}.WithT(0),
                                               Point(0, 1, 0).WithT(0), Point{1, 1, 0}.WithT(0),
                                               Point{0, 0, 0}.WithT(1), Point{1, 0, 0}.WithT(1),
                                               Point{0, 1, 0}.WithT(1), Point{1, 1, 0}.WithT(1)}},
        ShapeTestParameters{LAGRANGE, 1, SPACETIME_QUADRILATERAL,
                            std::vector<Point>{Point(0, 0, 0).WithT(0), Point{1, 0, 0}.WithT(0),
                                               Point(0, 1, 0).WithT(0), Point{1, 1, 0}.WithT(0),
                                               Point{0, 0, 0}.WithT(1), Point{1, 0, 0}.WithT(1),
                                               Point{0, 1, 0}.WithT(1), Point{1, 1, 0}.WithT(1)}},
        ShapeTestParameters{LAGRANGE, 2, SPACETIME_QUADRILATERAL,
                            std::vector<Point>{Point(0, 0).WithT(0),       Point(0.5, 0).WithT(0),
                                               Point{1, 0}.WithT(0),       Point(0, 0.5).WithT(0),
                                               Point(0.5, 0.5).WithT(0),   Point{1, 0.5}.WithT(0),
                                               Point(0, 1).WithT(0),       Point(0.5, 1).WithT(0),
                                               Point{1, 1}.WithT(0),

                                               Point(0, 0).WithT(0.5),     Point(0.5, 0).WithT(0.5),
                                               Point{1, 0}.WithT(0.5),     Point(0, 0.5).WithT(0.5),
                                               Point(0.5, 0.5).WithT(0.5), Point{1, 0.5}.WithT(0.5),
                                               Point(0, 1).WithT(0.5),     Point(0.5, 1).WithT(0.5),
                                               Point{1, 1}.WithT(0.5),

                                               Point(0, 0).WithT(1),       Point(0.5, 0).WithT(1),
                                               Point{1, 0}.WithT(1),       Point(0, 0.5).WithT(1),
                                               Point(0.5, 0.5).WithT(1),   Point{1, 0.5}.WithT(1),
                                               Point(0, 1).WithT(1),       Point(0.5, 1).WithT(1),
                                               Point{1, 1}.WithT(1)}}));

int main(int argc, char **argv) {
  MppTest mppTest = MppTestBuilder(argc, argv).WithPPM().WithScreenLogging();
  return mppTest.RUN_ALL_MPP_TESTS();
}
