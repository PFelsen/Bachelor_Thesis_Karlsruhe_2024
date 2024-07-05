#include "ArgyrisDoF.hpp"
#include "TestDoF.hpp"

class ArgyrisDoFTest : public DoFTest {
protected:
  ArgyrisDoFTest(CELLTYPE type) : DoFTest(type, "ArgyrisDoF") { doF = new ArgyrisDoF(); }
};

class ArgyrisDoFTriangleTest : public ArgyrisDoFTest {
protected:
  ArgyrisDoFTriangleTest() : ArgyrisDoFTest(TRIANGLE) {}
};

DOF_TESTS(ArgyrisDoFTriangleTest)

#if SpaceDimension == 1
GTEST_ALLOW_UNINSTANTIATED_PARAMETERIZED_TEST(ArgyrisDoFTriangleTest);
#endif

#if SpaceDimension >= 2
INSTANTIATE_TEST_SUITE_P(ArgyrisDoFTest, ArgyrisDoFTriangleTest,
                         Values(DoFData{-1,
                                        {Point(0.0, 0.0), Point(1.0, 0.0), Point(0.0, 1.0),
                                         Point(0.5, 0.0), Point(0.5, 0.5), Point(0.0, 0.5)},
                                        {6, 6, 6, 1, 1, 1},
                                        {{0, 1, 3}, {1, 2, 4}, {2, 0, 5}}}));
#endif

int main(int argc, char **argv) {
  MppTest mppTest = MppTestBuilder(argc, argv).WithPPM();
  return mppTest.RUN_ALL_MPP_TESTS();
}
