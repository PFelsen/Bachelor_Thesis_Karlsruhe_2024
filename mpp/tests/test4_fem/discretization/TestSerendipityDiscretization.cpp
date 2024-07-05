#include "TestDiscretization.hpp"

#include "SerendipityDiscretization.hpp"

class SerendipityDiscretizationTest : public DiscretizationTest {
  std::string getDiscName(int degree) {
    if (degree == 2) return "SerendipityP" + std::to_string(degree);
    return "P" + std::to_string(degree);
  }

  std::string getDoFName(int degree) {
    if (degree == 2) return "SerendipityDoF";
    return "LagrangeDoF (degree=" + std::to_string(degree) + ")";
  }
protected:
  SerendipityDiscretizationTest(std::string meshesName) :
      DiscretizationTest(meshesName, getDiscName(GetParam()), getDoFName(GetParam())) {
    disc = new SerendipityDiscretization(*meshes, GetParam());
  }
};

class SerendipityDiscretizationIntervalTest : public SerendipityDiscretizationTest {
protected:
  SerendipityDiscretizationIntervalTest() : SerendipityDiscretizationTest("Interval_DirichletBC") {
    shapeNames = {{INTERVAL, {{0, {"P" + std::to_string(GetParam()) + "Int"}}}}};
    quadNames = {{INTERVAL, {{0, {"Qint: exact up to degree " + std::to_string(2 * GetParam())}}}}};
    faceQuadNames = {{INTERVAL, {{0, "Qpoint"}}}};
  }
};

DISCRETIZATION_TESTS(SerendipityDiscretizationIntervalTest)

INSTANTIATE_TEST_CASE_P(SerendipityDiscretizationTest, SerendipityDiscretizationIntervalTest,
                        Values(0, 1, 2, 3, 4, 5, 6, 7, 8));

#if SpaceDimension >= 2

class SerendipityDiscretizationTriangleTest : public SerendipityDiscretizationTest {
protected:
  SerendipityDiscretizationTriangleTest() : SerendipityDiscretizationTest("Triangle_DirichletBC") {
    shapeNames = {{TRIANGLE, {{0, {"P" + std::to_string(GetParam()) + "Tri"}}}}};
    quadNames = {{TRIANGLE, {{0, {"Qtri: exact up to degree " + std::to_string(2 * GetParam())}}}}};
    faceQuadNames = {
        {TRIANGLE, {{0, "Qint: exact up to degree " + std::to_string(2 * GetParam())}}}};
  }
};

DISCRETIZATION_TESTS(SerendipityDiscretizationTriangleTest)

INSTANTIATE_TEST_CASE_P(SerendipityDiscretizationTest, SerendipityDiscretizationTriangleTest,
                        Values(0, 1, 2, 3, 4, 5, 6));

class SerendipityDiscretizationQuadrilateralTest : public SerendipityDiscretizationTest {
protected:
  SerendipityDiscretizationQuadrilateralTest() :
      SerendipityDiscretizationTest("Square_DirichletBC") {
    if (GetParam() == 2) {
      shapeNames = {{QUADRILATERAL, {{0, {"P" + std::to_string(GetParam()) + "QuadSerendipity"}}}}};
    } else {
      shapeNames = {{QUADRILATERAL, {{0, {"P" + std::to_string(GetParam()) + "Quad"}}}}};
    }
    quadNames = {{QUADRILATERAL,
                  {{0,
                    {"Qquad: exact up to degree " + std::to_string(2 * GetParam()) + " in x and "
                     + std::to_string(2 * GetParam()) + " in y"}}}}};
    faceQuadNames = {
        {QUADRILATERAL, {{0, "Qint: exact up to degree " + std::to_string(2 * GetParam())}}}};
  }
};

DISCRETIZATION_TESTS(SerendipityDiscretizationQuadrilateralTest)

INSTANTIATE_TEST_CASE_P(SerendipityDiscretizationTest, SerendipityDiscretizationQuadrilateralTest,
                        Values(0, 1, 2, 3, 4, 5, 6, 7, 8));

class SerendipityDiscretizationQuadrilateralTriangleTest : public SerendipityDiscretizationTest {
protected:
  SerendipityDiscretizationQuadrilateralTriangleTest() :
      SerendipityDiscretizationTest("SquareTriangle_DirichletBC") {
    if (GetParam() == 2) {
      shapeNames = {{QUADRILATERAL, {{0, {"P" + std::to_string(GetParam()) + "QuadSerendipity"}}}},
                    {TRIANGLE, {{0, {"P" + std::to_string(GetParam()) + "Tri"}}}}};
    } else {
      shapeNames = {{QUADRILATERAL, {{0, {"P" + std::to_string(GetParam()) + "Quad"}}}},
                    {TRIANGLE, {{0, {"P" + std::to_string(GetParam()) + "Tri"}}}}};
    }
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

DISCRETIZATION_TESTS(SerendipityDiscretizationQuadrilateralTriangleTest)

INSTANTIATE_TEST_CASE_P(SerendipityDiscretizationTest,
                        SerendipityDiscretizationQuadrilateralTriangleTest,
                        Values(0, 1, 2, 3, 4, 5, 6));

#endif
#if SpaceDimension >= 3

class SerendipityDiscretizationTetrahedronTest : public SerendipityDiscretizationTest {
protected:
  SerendipityDiscretizationTetrahedronTest() :
      SerendipityDiscretizationTest("Tetrahedron_DirichletBC") {
    shapeNames = {{TETRAHEDRON, {{0, {"P" + std::to_string(GetParam()) + "Tet"}}}}};
    quadNames = {
        {TETRAHEDRON, {{0, {"Qtet: exact up to degree " + std::to_string(2 * GetParam())}}}}};
    faceQuadNames = {
        {TETRAHEDRON, {{0, "Qtri: exact up to degree " + std::to_string(2 * GetParam())}}}};
  }
};

DISCRETIZATION_TESTS(SerendipityDiscretizationTetrahedronTest)

INSTANTIATE_TEST_CASE_P(SerendipityDiscretizationTest, SerendipityDiscretizationTetrahedronTest,
                        Values(0, 1, 2, 3));

class SerendipityDiscretizationHexahedronTest : public SerendipityDiscretizationTest {
protected:
  SerendipityDiscretizationHexahedronTest() :
      SerendipityDiscretizationTest("Hexahedron_DirichletBC") {
    if (GetParam() == 2) {
      shapeNames = {{HEXAHEDRON, {{0, {"P" + std::to_string(GetParam()) + "HexSerendipity"}}}}};
    } else {
      shapeNames = {{HEXAHEDRON, {{0, {"P" + std::to_string(GetParam()) + "Hex"}}}}};
    }
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

DISCRETIZATION_TESTS(SerendipityDiscretizationHexahedronTest)

INSTANTIATE_TEST_CASE_P(SerendipityDiscretizationTest, SerendipityDiscretizationHexahedronTest,
                        Values(0, 1, 2, 3, 4, 5, 6, 7, 8));

#endif

int main(int argc, char **argv) {
  MppTest mppTest = MppTestBuilder(argc, argv).WithPPM();
  return mppTest.RUN_ALL_MPP_TESTS();
}