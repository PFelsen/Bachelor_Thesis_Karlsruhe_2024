#include "RTDoF.hpp"
#include "TestDoF.hpp"

class RTDoFTest : public DoFTest {
protected:
  RTDoFTest(CELLTYPE type) :
      DoFTest(type, "RaviartThomasDoF (order=" + std::to_string(GetParam().degree) + ")") {
    doF = new RTDoF(GetParam().degree);
  }
};

class RTDoFIntervalTest : public RTDoFTest {
protected:
  RTDoFIntervalTest() : RTDoFTest(INTERVAL) {}
};

DOF_TESTS(RTDoFIntervalTest)

INSTANTIATE_TEST_CASE_P(
    RTDoFTest, RTDoFIntervalTest,
    Values(DoFData{0, {Point(0.0), Point(1.0)}, vector<short>(2, 1), {{0}, {1}}}));

#if SpaceDimension >= 2

class RTDoFTriangleTest : public RTDoFTest {
protected:
  RTDoFTriangleTest() : RTDoFTest(TRIANGLE) {}
};

DOF_TESTS(RTDoFTriangleTest)

INSTANTIATE_TEST_CASE_P(
    RTDoFTest, RTDoFTriangleTest,
    Values(DoFData{0,
                   {Point(0.5, 0.0), Point(0.5, 0.5), Point(0.0, 0.5)},
                   vector<short>(3, 1),
                   {{0}, {1}, {2}}},
           DoFData{1,
                   {Point(1.0 / 3.0, 0.0), Point(2.0 / 3.0, 1.0 / 3.0), Point(0.0, 2.0 / 3.0),
                    Point(2.0 / 3.0, 0.0), Point(1.0 / 3.0, 2.0 / 3.0), Point(0.0, 1.0 / 3.0),
                    Point(1.0 / 3.0, 1.0 / 3.0)},
                   {1, 1, 1, 1, 1, 1, 2},
                   {{0, 3}, {1, 4}, {2, 5}}},
           DoFData{2,
                   {Point(0.25, 0.0), Point(0.75, 0.25), Point(0.0, 0.75), Point(0.5, 0.0),
                    Point(0.5, 0.5), Point(0.0, 0.5), Point(0.75, 0.0), Point(0.25, 0.75),
                    Point(0.0, 0.25), Point(0.25, 0.25), Point(0.5, 0.25), Point(0.25, 0.5)},
                   {1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2},
                   {{0, 3, 6}, {1, 4, 7}, {2, 5, 8}}},
           DoFData{3,
                   {Point(0.2, 0, 0), Point(0.8, 0.2, 0), Point(0, 0.8, 0), Point(0.4, 0, 0),
                    Point(0.6, 0.4, 0), Point(0, 0.6, 0), Point(0.6, 0, 0), Point(0.4, 0.6, 0),
                    Point(0, 0.4, 0), Point(0.8, 0, 0), Point(0.2, 0.8, 0), Point(0, 0.2, 0),
                    Point(0.2, 0.2, 0), Point(0.6, 0.2, 0), Point(0.2, 0.6, 0), Point(0.4, 0.2, 0),
                    Point(0.4, 0.4, 0), Point(0.2, 0.4, 0)},
                   {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2},
                   {{0, 3, 6, 9}, {1, 4, 7, 10}, {2, 5, 8, 11}}}));

class RTDoFQuadrilateralTest : public RTDoFTest {
protected:
  RTDoFQuadrilateralTest() : RTDoFTest(QUADRILATERAL) {}
};

DOF_TESTS(RTDoFQuadrilateralTest)

INSTANTIATE_TEST_CASE_P(RTDoFTest, RTDoFQuadrilateralTest,
                        Values(DoFData{0,
                                       {Point(0.5, 0.0), Point(1.0, 0.5), Point(0.5, 1.0),
                                        Point(0.0, 0.5)},
                                       vector<short>(4, 1),
                                       {{0}, {1}, {2}, {3}}}));

#endif
#if SpaceDimension >= 3

class RTDoFTetrahedronTest : public RTDoFTest {
protected:
  RTDoFTetrahedronTest() : RTDoFTest(TETRAHEDRON) {}
};

DOF_TESTS(RTDoFTetrahedronTest)

INSTANTIATE_TEST_CASE_P(
    RTDoFTest, RTDoFTetrahedronTest,
    Values(DoFData{0,
                   {Point(1.0 / 3.0, 1.0 / 3.0, 0.0), Point(1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0),
                    Point(0.0, 1.0 / 3.0, 1.0 / 3.0), Point(1.0 / 3.0, 0.0, 1.0 / 3.0)},
                   vector<short>(4, 1),
                   {{0}, {1}, {2}, {3}}}));

class RTDoFHexahedronTest : public RTDoFTest {
protected:
  RTDoFHexahedronTest() : RTDoFTest(HEXAHEDRON) {}
};

DOF_TESTS(RTDoFHexahedronTest)

INSTANTIATE_TEST_CASE_P(RTDoFTest, RTDoFHexahedronTest,
                        Values(DoFData{0,
                                       {Point(0.5, 0.5, 0.0), Point(0.5, 0.0, 0.5),
                                        Point(1.0, 0.5, 0.5), Point(0.5, 1.0, 0.5),
                                        Point(0.0, 0.5, 0.5), Point(0.5, 0.5, 1.0)},
                                       vector<short>(6, 1),
                                       {{0}, {1}, {2}, {3}, {4}, {5}}}));

#endif

int main(int argc, char **argv) {
  MppTest mppTest = MppTestBuilder(argc, argv).WithPPM();
  return mppTest.RUN_ALL_MPP_TESTS();
}
