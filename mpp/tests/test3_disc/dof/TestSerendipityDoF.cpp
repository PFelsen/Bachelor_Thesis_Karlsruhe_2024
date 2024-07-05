#include "SerendipityDoF.hpp"
#include "TestDoF.hpp"

class SerendipityDoFTest : public DoFTest {
protected:
  SerendipityDoFTest(CELLTYPE type) : DoFTest(type, "SerendipityDoF") {
    doF = new SerendipityDoF();
  }
};

class SerendipityDoFIntervalTest : public SerendipityDoFTest {
protected:
  SerendipityDoFIntervalTest() : SerendipityDoFTest(INTERVAL) {}
};

DOF_TESTS(SerendipityDoFIntervalTest)

INSTANTIATE_TEST_SUITE_P(SerendipityDoFTest, SerendipityDoFIntervalTest,
                         Values(DoFData{2,
                                        {Point(0.0 / 2.0), Point(1.0 / 2.0), Point(2.0 / 2.0)},
                                        vector<short>(3, 1),
                                        {{0}, {2}}}));

#if SpaceDimension >= 2

class SerendipityDoFTriangleTest : public SerendipityDoFTest {
protected:
  SerendipityDoFTriangleTest() : SerendipityDoFTest(TRIANGLE) {}
};

DOF_TESTS(SerendipityDoFTriangleTest)

INSTANTIATE_TEST_CASE_P(SerendipityDoFTest, SerendipityDoFTriangleTest,
                        Values(DoFData{2,
                                       {Point(0.0, 0.0), Point(1.0, 0.0), Point(0.0, 1.0),
                                        Point(0.5, 0.0), Point(0.5, 0.5), Point(0.0, 0.5)},
                                       vector<short>(6, 1),
                                       {{0, 3, 1}, {1, 4, 2}, {2, 5, 0}}}));

class SerendipityDoFQuadrilateralTest : public SerendipityDoFTest {
protected:
  SerendipityDoFQuadrilateralTest() : SerendipityDoFTest(QUADRILATERAL) {}
};

DOF_TESTS(SerendipityDoFQuadrilateralTest)

INSTANTIATE_TEST_CASE_P(SerendipityDoFTest, SerendipityDoFQuadrilateralTest,
                        Values(DoFData{2,
                                       {Point(0.0, 0.0), Point(1.0, 0.0), Point(1.0, 1.0),
                                        Point(0.0, 1.0), Point(0.5, 0.0), Point(1.0, 0.5),
                                        Point(0.5, 1.0), Point(0.0, 0.5)},
                                       vector<short>(8, 1),
                                       {{0, 1, 4}, {1, 2, 5}, {2, 3, 6}, {3, 0, 7}}}));

#endif
#if SpaceDimension >= 3

class SerendipityDoFTetrahedronTest : public SerendipityDoFTest {
protected:
  SerendipityDoFTetrahedronTest() : SerendipityDoFTest(TETRAHEDRON) {}
};

DOF_TESTS(SerendipityDoFTetrahedronTest)

INSTANTIATE_TEST_CASE_P(
    SerendipityDoFTest, SerendipityDoFTetrahedronTest,
    Values(
        DoFData{2,
                {Point(0.0, 0.0, 0.0), Point(0.5, 0.0, 0.0), Point(1.0, 0.0, 0.0),
                 Point(0.0, 0.5, 0.0), Point(0.5, 0.5, 0.0), Point(0.0, 1.0, 0.0),
                 Point(0.0, 0.0, 0.5), Point(0.5, 0.0, 0.5), Point(0.0, 0.5, 0.5),
                 Point(0.0, 0.0, 1.0)},
                vector<short>(10, 1),
                {{0, 3, 5, 1, 4, 2}, {2, 4, 5, 7, 8, 9}, {0, 6, 9, 3, 8, 5}, {0, 1, 2, 6, 7, 9}}}));

class SerendipityDoFHexahedronTest : public SerendipityDoFTest {
protected:
  SerendipityDoFHexahedronTest() : SerendipityDoFTest(HEXAHEDRON) {}
};

DOF_TESTS(SerendipityDoFHexahedronTest)

INSTANTIATE_TEST_CASE_P(
    SerendipityDoFTest, SerendipityDoFHexahedronTest,
    Values(DoFData{2,
                   {Point(0.0, 0.0, 0.0), Point(1.0, 0.0, 0.0), Point(1.0, 1.0, 0.0),
                    Point(0.0, 1.0, 0.0), Point(0.0, 0.0, 1.0), Point(1.0, 0.0, 1.0),
                    Point(1.0, 1.0, 1.0), Point(0.0, 1.0, 1.0), Point(0.5, 0.0, 0.0),
                    Point(1.0, 0.5, 0.0), Point(0.5, 1.0, 0.0), Point(0.0, 0.5, 0.0),
                    Point(0.0, 0.0, 0.5), Point(1.0, 0.0, 0.5), Point(1.0, 1.0, 0.5),
                    Point(0.0, 1.0, 0.5), Point(0.5, 0.0, 1.0), Point(1.0, 0.5, 1.0),
                    Point(0.5, 1.0, 1.0), Point(0.0, 0.5, 1.0)},
                   vector<short>(20, 1),
                   {{0, 3, 2, 1, 8, 9, 10, 11},
                    {0, 1, 5, 4, 8, 12, 13, 16},
                    {1, 2, 6, 5, 9, 13, 14, 17},
                    {2, 3, 7, 6, 10, 14, 15, 18},
                    {3, 0, 4, 7, 11, 15, 12, 19},
                    {4, 5, 6, 7, 16, 17, 18, 19}}}));

#endif

int main(int argc, char **argv) {
  MppTest mppTest = MppTestBuilder(argc, argv).WithPPM();
  return mppTest.RUN_ALL_MPP_TESTS();
}