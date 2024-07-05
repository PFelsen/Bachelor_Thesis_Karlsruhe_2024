#include "TestDiscretization.hpp"

#include "RTLagrangeDiscretization.hpp"

class RTLagrangeDiscretizationTest : public MixedDiscretizationTest {
protected:
  RTLagrangeDiscretizationTest(std::string meshesName) :
      MixedDiscretizationTest(meshesName,
                              "RT" + std::to_string(GetParam().first) + "_P"
                                  + std::to_string(GetParam().second),
                              "RaviartThomasDoF (order=" + std::to_string(GetParam().first) + "); "
                                  + "LagrangeDoF (degree=" + std::to_string(GetParam().second)
                                  + ")") {
    disc = new RTLagrangeDiscretization(*meshes, GetParam().first, GetParam().second);
  }
};

class RTLagrangeDiscretizationIntervalTest : public RTLagrangeDiscretizationTest {
protected:
  RTLagrangeDiscretizationIntervalTest() : RTLagrangeDiscretizationTest("Interval_DirichletBC") {
    shapeNames = {{INTERVAL,
                   {{0, {"RT" + std::to_string(GetParam().first) + "Int"}},
                    {1, {"P" + std::to_string(GetParam().second) + "Int"}}}}};
    int quadExactUpTo = std::max(2 * (GetParam().first + 1), 2 * GetParam().second);
    quadNames = {{INTERVAL,
                  {{0, {"Qint: exact up to degree " + std::to_string(quadExactUpTo)}},
                   {1, {"Qint: exact up to degree " + std::to_string(quadExactUpTo)}}}}};
    faceQuadNames = {{INTERVAL, {{0, "Qpoint"}, {1, "Qpoint"}}}};
  }
};

DISCRETIZATION_TESTS(RTLagrangeDiscretizationIntervalTest)

INSTANTIATE_TEST_CASE_P(RTLagrangeDiscretizationTest, RTLagrangeDiscretizationIntervalTest,
                        Values(std::make_pair(0, 0)));

#if SpaceDimension >= 2

class RTLagrangeDiscretizationTriangleTest : public RTLagrangeDiscretizationTest {
protected:
  RTLagrangeDiscretizationTriangleTest() : RTLagrangeDiscretizationTest("Triangle_DirichletBC") {
    shapeNames = {{TRIANGLE,
                   {{0, {"RT" + std::to_string(GetParam().first) + "Tri"}},
                    {1, {"P" + std::to_string(GetParam().second) + "Tri"}}}}};
    int quadExactUpTo = std::max(2 * (GetParam().first + 1), 2 * GetParam().second);
    quadNames = {{TRIANGLE,
                  {{0, {"Qtri: exact up to degree " + std::to_string(quadExactUpTo)}},
                   {1, {"Qtri: exact up to degree " + std::to_string(quadExactUpTo)}}}}};
    faceQuadNames = {{TRIANGLE,
                      {{0, "Qint: exact up to degree " + std::to_string(quadExactUpTo)},
                       {1, "Qint: exact up to degree " + std::to_string(quadExactUpTo)}}}};
  }
};

DISCRETIZATION_TESTS(RTLagrangeDiscretizationTriangleTest)

INSTANTIATE_TEST_CASE_P(RTLagrangeDiscretizationTest, RTLagrangeDiscretizationTriangleTest,
                        Values(std::make_pair(0, 0), std::make_pair(1, 0), std::make_pair(1, 1),
                               std::make_pair(2, 1)));

class RTLagrangeDiscretizationQuadrilateralTest : public RTLagrangeDiscretizationTest {
protected:
  RTLagrangeDiscretizationQuadrilateralTest() : RTLagrangeDiscretizationTest("Square_DirichletBC") {
    shapeNames = {{QUADRILATERAL,
                   {{0, {"RT" + std::to_string(GetParam().first) + "Quad"}},
                    {1, {"P" + std::to_string(GetParam().second) + "Quad"}}}}};
    int quadExactUpTo = std::max(2 * (GetParam().first + 1), 2 * GetParam().second);
    quadNames = {{QUADRILATERAL,
                  {{0,
                    {{"Qquad: exact up to degree " + std::to_string(quadExactUpTo) + " in x and "
                      + std::to_string(quadExactUpTo) + " in y"}}},
                   {1,
                    {"Qquad: exact up to degree " + std::to_string(quadExactUpTo) + " in x and "
                     + std::to_string(quadExactUpTo) + " in y"}}}}};
    faceQuadNames = {{QUADRILATERAL,
                      {{0, "Qint: exact up to degree " + std::to_string(quadExactUpTo)},
                       {1, "Qint: exact up to degree " + std::to_string(quadExactUpTo)}}}};
  }
};

DISCRETIZATION_TESTS(RTLagrangeDiscretizationQuadrilateralTest)

INSTANTIATE_TEST_CASE_P(RTLagrangeDiscretizationTest, RTLagrangeDiscretizationQuadrilateralTest,
                        Values(std::make_pair(0, 0)));

class RTLagrangeDiscretizationQuadrilateralTriangleTest : public RTLagrangeDiscretizationTest {
protected:
  RTLagrangeDiscretizationQuadrilateralTriangleTest() :
      RTLagrangeDiscretizationTest("SquareTriangle_DirichletBC") {
    shapeNames = {{TRIANGLE,
                   {{0, {"RT" + std::to_string(GetParam().first) + "Tri"}},
                    {1, {"P" + std::to_string(GetParam().second) + "Tri"}}}},
                  {QUADRILATERAL,
                   {{0, {"RT" + std::to_string(GetParam().first) + "Quad"}},
                    {1, {"P" + std::to_string(GetParam().second) + "Quad"}}}}};
    int quadExactUpTo = std::max(2 * (GetParam().first + 1), 2 * GetParam().second);
    quadNames = {{TRIANGLE,
                  {{0, {"Qtri: exact up to degree " + std::to_string(quadExactUpTo)}},
                   {1, {"Qtri: exact up to degree " + std::to_string(quadExactUpTo)}}}},
                 {QUADRILATERAL,
                  {{0,
                    {{"Qquad: exact up to degree " + std::to_string(quadExactUpTo) + " in x and "
                      + std::to_string(quadExactUpTo) + " in y"}}},
                   {1,
                    {"Qquad: exact up to degree " + std::to_string(quadExactUpTo) + " in x and "
                     + std::to_string(quadExactUpTo) + " in y"}}}}};
    faceQuadNames = {{TRIANGLE,
                      {{0, "Qint: exact up to degree " + std::to_string(quadExactUpTo)},
                       {1, "Qint: exact up to degree " + std::to_string(quadExactUpTo)}}},
                     {QUADRILATERAL,
                      {{0, "Qint: exact up to degree " + std::to_string(quadExactUpTo)},
                       {1, "Qint: exact up to degree " + std::to_string(quadExactUpTo)}}}};
  }
};

DISCRETIZATION_TESTS(RTLagrangeDiscretizationQuadrilateralTriangleTest)

INSTANTIATE_TEST_CASE_P(RTLagrangeDiscretizationTest,
                        RTLagrangeDiscretizationQuadrilateralTriangleTest,
                        Values(std::make_pair(0, 0)));

#endif
#if SpaceDimension >= 3

class RTLagrangeDiscretizationHexahedronTest : public RTLagrangeDiscretizationTest {
protected:
  RTLagrangeDiscretizationHexahedronTest() :
      RTLagrangeDiscretizationTest("Hexahedron_DirichletBC") {
    shapeNames = {{HEXAHEDRON,
                   {{0, {"RT" + std::to_string(GetParam().first) + "Hex"}},
                    {0, {"P" + std::to_string(GetParam().second) + "Hex"}}}}};
    int quadExactUpTo = std::max(2 * (GetParam().first + 1), 2 * GetParam().second);
    quadNames = {{HEXAHEDRON,
                  {{0,
                    {{"Qhex: exact up to degree " + std::to_string(quadExactUpTo) + " in x and "
                      + std::to_string(quadExactUpTo) + " in y and " + std::to_string(quadExactUpTo)
                      + " in z"}}},
                   {1,
                    {"Qhex: exact up to degree " + std::to_string(quadExactUpTo) + " in x and "
                     + std::to_string(quadExactUpTo) + " in y and " + std::to_string(quadExactUpTo)
                     + " in z"}}}}};
    faceQuadNames = {{HEXAHEDRON,
                      {{0, "Qquad: exact up to degree " + std::to_string(quadExactUpTo)
                               + " in x and " + std::to_string(quadExactUpTo) + " in y"},
                       {1, "Qquad: exact up to degree " + std::to_string(quadExactUpTo)
                               + " in x and " + std::to_string(quadExactUpTo) + " in y"}}}};
  }
};

DISCRETIZATION_TESTS(RTLagrangeDiscretizationHexahedronTest)

INSTANTIATE_TEST_CASE_P(RTLagrangeDiscretizationTest, RTLagrangeDiscretizationHexahedronTest,
                        Values(std::make_pair(0, 0)));

#endif

int main(int argc, char **argv) {
  MppTest mppTest = MppTestBuilder(argc, argv).WithPPM();
  return mppTest.RUN_ALL_MPP_TESTS();
}