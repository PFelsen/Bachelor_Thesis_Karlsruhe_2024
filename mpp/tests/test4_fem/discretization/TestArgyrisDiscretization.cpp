#include "TestDiscretization.hpp"

#include "ArgyrisDiscretization.hpp"

class ArgyrisDiscretizationTest : public DiscretizationTest {
protected:
  ArgyrisDiscretizationTest(std::string meshesName) :
      DiscretizationTest(meshesName, "Argyris", "ArgyrisDoF") {
    disc = new ArgyrisDiscretization(*meshes);
  }
};

class ArgyrisDiscretizationTriangleTest : public ArgyrisDiscretizationTest {
protected:
  ArgyrisDiscretizationTriangleTest() : ArgyrisDiscretizationTest("Triangle_DirichletBC") {
    shapeNames = {{TRIANGLE, {{0, {"ArgyrisTri"}}}}};
    quadNames = {{TRIANGLE, {{0, {"Qtri: exact up to degree 10"}}}}};
    faceQuadNames = {{TRIANGLE, {{0, "Qint: exact up to degree 10"}}}};
  }
};

DISCRETIZATION_TESTS(ArgyrisDiscretizationTriangleTest)

#if SpaceDimension == 1
GTEST_ALLOW_UNINSTANTIATED_PARAMETERIZED_TEST(ArgyrisDiscretizationTriangleTest);
#endif

#if SpaceDimension >= 2
INSTANTIATE_TEST_SUITE_P(ArgyrisDiscretizationTest, ArgyrisDiscretizationTriangleTest, Values(-1));
#endif

int main(int argc, char **argv) {
  MppTest mppTest = MppTestBuilder(argc, argv).WithPPM();
  return mppTest.RUN_ALL_MPP_TESTS();
}