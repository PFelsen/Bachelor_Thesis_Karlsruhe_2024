#include "TestDiscretization.hpp"

#include "RTDiscretization.hpp"

class RTDiscretizationTest : public DiscretizationTest {
protected:
  RTDiscretizationTest(std::string meshesName) :
      DiscretizationTest(meshesName, "RT" + std::to_string(GetParam()),
                         "RaviartThomasDoF (order=" + std::to_string(GetParam()) + ")") {
    disc = new RTDiscretization(*meshes, GetParam());
  }
};

class RTDiscretizationIntervalTest : public RTDiscretizationTest {
protected:
  RTDiscretizationIntervalTest() : RTDiscretizationTest("Interval_DirichletBC") {
    shapeNames = {{INTERVAL, {{0, {"RT" + std::to_string(GetParam()) + "Int"}}}}};
    quadNames = {{INTERVAL, {{0, {"Qint: exact up to degree 2"}}}}};
    faceQuadNames = {{INTERVAL, {{0, "Qpoint"}}}};
  }
};

DISCRETIZATION_TESTS(RTDiscretizationIntervalTest)

INSTANTIATE_TEST_CASE_P(RTDiscretizationTest, RTDiscretizationIntervalTest, Values(0));

#if SpaceDimension >= 2

class RTDiscretizationTriangleTest : public RTDiscretizationTest {
protected:
  RTDiscretizationTriangleTest() : RTDiscretizationTest("Triangle_DirichletBC") {
    shapeNames = {{TRIANGLE, {{0, {"RT" + std::to_string(GetParam()) + "Tri"}}}}};
    quadNames = {
        {TRIANGLE, {{0, {"Qtri: exact up to degree " + std::to_string(2 * (GetParam() + 1))}}}}};
    faceQuadNames = {
        {TRIANGLE, {{0, "Qint: exact up to degree " + std::to_string(2 * (GetParam() + 1))}}}};
  }
};

DISCRETIZATION_TESTS(RTDiscretizationTriangleTest)

INSTANTIATE_TEST_CASE_P(RTDiscretizationTest, RTDiscretizationTriangleTest,
                        Values(0, 1, 2, 3, 4, 5, 6));

class RTDiscretizationQuadrilateralTest : public RTDiscretizationTest {
protected:
  RTDiscretizationQuadrilateralTest() : RTDiscretizationTest("Square_DirichletBC") {
    shapeNames = {{QUADRILATERAL, {{0, {"RT" + std::to_string(GetParam()) + "Quad"}}}}};
    quadNames = {{QUADRILATERAL, {{0, {"Qquad: exact up to degree 2 in x and 2 in y"}}}}};
    faceQuadNames = {{QUADRILATERAL, {{0, "Qint: exact up to degree 2"}}}};
  }
};

DISCRETIZATION_TESTS(RTDiscretizationQuadrilateralTest)

INSTANTIATE_TEST_CASE_P(RTDiscretizationTest, RTDiscretizationQuadrilateralTest, Values(0));

class RTDiscretizationQuadrilateralTriangleTest : public RTDiscretizationTest {
protected:
  RTDiscretizationQuadrilateralTriangleTest() : RTDiscretizationTest("SquareTriangle_DirichletBC") {
    shapeNames = {
        {TRIANGLE, {{0, {"RT" + std::to_string(GetParam()) + "Tri"}}}},
        {QUADRILATERAL, {{0, {"RT" + std::to_string(GetParam()) + "Quad"}}}},
    };
    quadNames = {{TRIANGLE,
                  {{0, {"Qtri: exact up to degree " + std::to_string(2 * (GetParam() + 1))}}}},
                 {QUADRILATERAL, {{0, {"Qquad: exact up to degree 2 in x and 2 in y"}}}}};
    faceQuadNames = {{TRIANGLE,
                      {{0, "Qint: exact up to degree " + std::to_string(2 * (GetParam() + 1))}}},
                     {QUADRILATERAL, {{0, "Qint: exact up to degree 2"}}}};
  }
};


DISCRETIZATION_TESTS(RTDiscretizationQuadrilateralTriangleTest)

INSTANTIATE_TEST_CASE_P(RTDiscretizationTest, RTDiscretizationQuadrilateralTriangleTest, Values(0));

#endif
#if SpaceDimension >= 3

class RTDiscretizationHexahedronTest : public RTDiscretizationTest {
protected:
  RTDiscretizationHexahedronTest() : RTDiscretizationTest("Hexahedron_DirichletBC") {
    shapeNames = {{HEXAHEDRON, {{0, {"RT" + std::to_string(GetParam()) + "Hex"}}}}};
    quadNames = {{HEXAHEDRON, {{0, {"Qhex: exact up to degree 2 in x and 2 in y and 2 in z"}}}}};
    faceQuadNames = {{HEXAHEDRON, {{0, "Qquad: exact up to degree 2 in x and 2 in y"}}}};
  }
};

DISCRETIZATION_TESTS(RTDiscretizationHexahedronTest)

INSTANTIATE_TEST_CASE_P(RTDiscretizationTest, RTDiscretizationHexahedronTest, Values(0));

#endif

int main(int argc, char **argv) {
  MppTest mppTest = MppTestBuilder(argc, argv).WithPPM();
  return mppTest.RUN_ALL_MPP_TESTS();
}