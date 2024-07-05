#include "TestDiscretization.hpp"

#include "LagrangeDiscretization.hpp"

class LagrangeDiscretizationTest : public DiscretizationTest {
protected:
  LagrangeDiscretizationTest(std::string meshesName) :
      DiscretizationTest(meshesName, "P" + std::to_string(GetParam()),
                         "LagrangeDoF (degree=" + std::to_string(GetParam()) + ")") {
    disc = new LagrangeDiscretization(*meshes, GetParam());
  }
};

class LagrangeDiscretizationIntervalTest : public LagrangeDiscretizationTest {
protected:
  LagrangeDiscretizationIntervalTest() : LagrangeDiscretizationTest("Interval_DirichletBC") {
    shapeNames = {{INTERVAL, {{0, {"P" + std::to_string(GetParam()) + "Int"}}}}};
    quadNames = {{INTERVAL, {{0, {"Qint: exact up to degree " + std::to_string(2 * GetParam())}}}}};
    faceQuadNames = {{INTERVAL, {{0, "Qpoint"}}}};
  }
};

DISCRETIZATION_TESTS(LagrangeDiscretizationIntervalTest)

INSTANTIATE_TEST_CASE_P(LagrangeDiscretizationTest, LagrangeDiscretizationIntervalTest,
                        Values(0, 1, 2, 3, 4, 5, 6, 7, 8));

#if SpaceDimension >= 2

class LagrangeDiscretizationTriangleTest : public LagrangeDiscretizationTest {
protected:
  LagrangeDiscretizationTriangleTest() : LagrangeDiscretizationTest("Triangle_DirichletBC") {
    shapeNames = {{TRIANGLE, {{0, {"P" + std::to_string(GetParam()) + "Tri"}}}}};
    quadNames = {{TRIANGLE, {{0, {"Qtri: exact up to degree " + std::to_string(2 * GetParam())}}}}};
    faceQuadNames = {
        {TRIANGLE, {{0, "Qint: exact up to degree " + std::to_string(2 * GetParam())}}}};
  }
};

DISCRETIZATION_TESTS(LagrangeDiscretizationTriangleTest)

INSTANTIATE_TEST_CASE_P(LagrangeDiscretizationTest, LagrangeDiscretizationTriangleTest,
                        Values(0, 1, 2, 3, 4, 5, 6));

class LagrangeDiscretizationQuadrilateralTest : public LagrangeDiscretizationTest {
protected:
  LagrangeDiscretizationQuadrilateralTest() : LagrangeDiscretizationTest("Square_DirichletBC") {
    shapeNames = {{QUADRILATERAL, {{0, {"P" + std::to_string(GetParam()) + "Quad"}}}}};
    quadNames = {{QUADRILATERAL,
                  {{0,
                    {"Qquad: exact up to degree " + std::to_string(2 * GetParam()) + " in x and "
                     + std::to_string(2 * GetParam()) + " in y"}}}}};
    faceQuadNames = {
        {QUADRILATERAL, {{0, "Qint: exact up to degree " + std::to_string(2 * GetParam())}}}};
  }
};

DISCRETIZATION_TESTS(LagrangeDiscretizationQuadrilateralTest)

INSTANTIATE_TEST_CASE_P(LagrangeDiscretizationTest, LagrangeDiscretizationQuadrilateralTest,
                        Values(0, 1, 2, 3, 4, 5, 6, 7, 8));

class LagrangeDiscretizationQuadrilateralTriangleTest : public LagrangeDiscretizationTest {
protected:
  LagrangeDiscretizationQuadrilateralTriangleTest() :
      LagrangeDiscretizationTest("SquareTriangle_DirichletBC") {
    shapeNames = {{QUADRILATERAL, {{0, {"P" + std::to_string(GetParam()) + "Quad"}}}},
                  {TRIANGLE, {{0, {"P" + std::to_string(GetParam()) + "Tri"}}}}};
    quadNames = {{QUADRILATERAL,
                  {{0,
                    {"Qquad: exact up to degree " + std::to_string(2 * GetParam()) + " in x and "
                     + std::to_string(2 * GetParam()) + " in y"}}}},
                 {TRIANGLE, {{0, {"Qtri: exact up to degree " + std::to_string(2 * GetParam())}}}}};
    faceQuadNames = {{QUADRILATERAL,
                      {{0, "Qint: exact up to degree " + std::to_string(2 * GetParam())}}},
                     {TRIANGLE,
                      {{0, "Qint: exact up to degree " + std::to_string(2 * GetParam())}}}};
  }
};

DISCRETIZATION_TESTS(LagrangeDiscretizationQuadrilateralTriangleTest)

INSTANTIATE_TEST_CASE_P(LagrangeDiscretizationTest, LagrangeDiscretizationQuadrilateralTriangleTest,
                        Values(0, 1, 2, 3, 4, 5, 6));

#endif
#if SpaceDimension >= 3

class LagrangeDiscretizationTetrahedronTest : public LagrangeDiscretizationTest {
protected:
  LagrangeDiscretizationTetrahedronTest() : LagrangeDiscretizationTest("Tetrahedron_DirichletBC") {
    shapeNames = {{TETRAHEDRON, {{0, {"P" + std::to_string(GetParam()) + "Tet"}}}}};
    quadNames = {
        {TETRAHEDRON, {{0, {"Qtet: exact up to degree " + std::to_string(2 * GetParam())}}}}};
    faceQuadNames = {
        {TETRAHEDRON, {{0, "Qtri: exact up to degree " + std::to_string(2 * GetParam())}}}};
  }
};

DISCRETIZATION_TESTS(LagrangeDiscretizationTetrahedronTest)

INSTANTIATE_TEST_CASE_P(LagrangeDiscretizationTest, LagrangeDiscretizationTetrahedronTest,
                        Values(0, 1, 2, 3));

class LagrangeDiscretizationHexahedronTest : public LagrangeDiscretizationTest {
protected:
  LagrangeDiscretizationHexahedronTest() : LagrangeDiscretizationTest("Hexahedron_DirichletBC") {
    shapeNames = {{HEXAHEDRON, {{0, {"P" + std::to_string(GetParam()) + "Hex"}}}}};
    quadNames = {{HEXAHEDRON,
                  {{0,
                    {"Qhex: exact up to degree " + std::to_string(2 * GetParam()) + " in x and "
                     + std::to_string(2 * GetParam()) + " in y and "
                     + std::to_string(2 * GetParam()) + " in z"}}}}};
    faceQuadNames = {{HEXAHEDRON,
                      {{0, "Qquad: exact up to degree " + std::to_string(2 * GetParam())
                               + " in x and " + std::to_string(2 * GetParam()) + " in y"}}}};
  }
};

DISCRETIZATION_TESTS(LagrangeDiscretizationHexahedronTest)

INSTANTIATE_TEST_CASE_P(LagrangeDiscretizationTest, LagrangeDiscretizationHexahedronTest,
                        Values(0, 1, 2, 3, 4, 5, 6, 7, 8));
#endif

int main(int argc, char **argv) {
  MppTest mppTest = MppTestBuilder(argc, argv).WithPPM();
  return mppTest.RUN_ALL_MPP_TESTS();
}