#ifndef TESTCHECKORIENTATION_HPP
#define TESTCHECKORIENTATION_HPP

#include "Point.hpp"

// correctly oriented corners are required for the test
struct CheckOrientationTestParameters {
  std::vector<Point> corners;
  bool expectSuccess;
};

void PrintTo(const CheckOrientationTestParameters &param, std::ostream *os) {
  *os << "(";
  for (const Point &p : param.corners) {
    *os << p << ", ";
  }
  *os << param.expectSuccess << ")";
}

static std::vector<CheckOrientationTestParameters>
    TestCasesInt{CheckOrientationTestParameters{{Point{0.0}, Point{1.0}}, true},
                 CheckOrientationTestParameters{{Point{0.5}, Point{0.75}}, true},
                 CheckOrientationTestParameters{{Point{0.0}, Point{GeometricTolerance / 10}},
                                                false},
#if SpaceDimension >= 2
                 CheckOrientationTestParameters{{Point{0.5, 1.0}, Point{0.75, 1.0}}, false},
#endif
#if SpaceDimension >= 3
                 CheckOrientationTestParameters{{Point{0.5, 1.0, 1.0}, Point{0.75, 1.0, 1.0}},
                                                false},
#endif
                 CheckOrientationTestParameters{{Point{0.0}, Point{0.5}, Point{1.0}}, false}};

static std::vector<CheckOrientationTestParameters>
    TestCasesTri{CheckOrientationTestParameters{{Point{0.0, 0.0}, Point{1.0, 0.0}, Point{0.0, 1.0}},
                                                true},
                 CheckOrientationTestParameters{{Point{0.0, 0.0}, Point{1.0, 0.0}, Point{0.5, 0.0}},
                                                false},
                 CheckOrientationTestParameters{{Point{0.0, 0.0}, Point{1.0, 0.0},
                                                 Point{GeometricTolerance / 10, 0.0}},
                                                false},
                 CheckOrientationTestParameters{{Point{0.0, 0.0}, Point{1.0, 0.0},
                                                 Point{0.5, GeometricTolerance / 10}},
                                                false},
                 CheckOrientationTestParameters{{Point{-1.0, 3.0}, Point{2.5, -1.0},
                                                 Point{0.5, 1.5}},
                                                true},
#if SpaceDimension >= 3
                 CheckOrientationTestParameters{{Point{0.0, 0.0, 1.0}, Point{1.0, 0.0, 1.0},
                                                 Point{0.0, 1.0, 1.0}},
                                                false},
#endif
                 CheckOrientationTestParameters{{Point{0.0, 0.0}, Point{1.0, 0.0}, Point{0.5, 0.5},
                                                 Point{0.0, 1.0}},
                                                false}};

static std::vector<CheckOrientationTestParameters>
    TestCasesQuad{CheckOrientationTestParameters{{Point{0.0, 0.0}, Point{1.0, 0.0}, Point{1.0, 1.0},
                                                  Point{0.0, 1.0}},
                                                 true},
                  CheckOrientationTestParameters{{Point{0.0, 0.0}, Point{1.0, 0.0}, Point{0.5, 0.5},
                                                  Point{0.0, 1.0}},
                                                 false},
                  CheckOrientationTestParameters{{Point{0.0, 0.0}, Point{1.0, 0.0},
                                                  Point{0.25, 0.25}, Point{0.0, 1.0}},
                                                 false},
#if SpaceDimension >= 3
                  CheckOrientationTestParameters{{Point{0.0, 0.0, 1.0}, Point{1.0, 0.0, 1.0},
                                                  Point{1.0, 1.0, 1.0}, Point{0.0, 1.0, 1.0}},
                                                 false},
#endif
                  CheckOrientationTestParameters{{Point{0.0, 0.0}, Point{1.0, 0.0}, Point{1.0, 1.0},
                                                  Point{0.0, 1.0}, Point{0.0, 0.5}},
                                                 false},
                  CheckOrientationTestParameters{{Point{0.0, 0.0}, Point{1.0, 0.0}, Point{1.0, 1.0},
                                                  Point{0.0, GeometricTolerance / 10}},
                                                 false},
                  CheckOrientationTestParameters{{Point{1.0, -1.0}, Point{3.5, 1.0},
                                                  Point{1.0, 2.0}, Point{-0.5, 1.0}},
                                                 true},
                  CheckOrientationTestParameters{{Point{1.0, -1.0}, Point{0.5, 1.0},
                                                  Point{1.0, 2.0}, Point{-0.5, 1.0}},
                                                 false}};

static std::vector<CheckOrientationTestParameters>
    TestCasesTet{CheckOrientationTestParameters{{Point{0.0, 0.0, 0.0}, Point{1.0, 0.0, 0.0},
                                                 Point{0.0, 1.0, 0.0}, Point{0.0, 0.0, 1.0}},
                                                true},
                 CheckOrientationTestParameters{{Point{0.0, 0.0, 0.0}, Point{1.0, 0.0, 0.0},
                                                 Point{0.0, 1.0, 0.0},
                                                 Point{0.0, 0.0, GeometricTolerance / 10}},
                                                false},
                 CheckOrientationTestParameters{{Point{0.0, 0.0, 0.0},
                                                 Point{GeometricTolerance / 10, 0.0, 0.0},
                                                 Point{0.0, 1.0, 0.0}, Point{0.0, 0.0, 1.0}},
                                                false},
                 CheckOrientationTestParameters{{Point{0.0, 0.0, 0.0}, Point{1.0, 0.0, 0.0},
                                                 Point{0.0, 1.0, 0.0}, Point{0.5, 0.0, 0.0}},
                                                false},
                 CheckOrientationTestParameters{{Point{0.0, 0.0, 0.0}, Point{1.0, 0.0, 0.0},
                                                 Point{0.0, 1.0, 0.0}, Point{0.0, 0.0, 1.0},
                                                 Point{0.25, 0.25, 0.25}},
                                                false},
                 CheckOrientationTestParameters{{Point{1.0, 0.0, 0.0}, Point{1.0, 1.0, 0.0},
                                                 Point{1.0, 1.0, 1.0}, Point{2.0, 1.0, 1.0}},
                                                true},
                 CheckOrientationTestParameters{{Point{0.0, 0.0, 1.0}, Point{1.0, 1.0, 1.0},
                                                 Point{0.0, 1.0, 1.0}, Point{1.0, 1.0, 2.0}},
                                                true}};

static std::vector<CheckOrientationTestParameters>
    TestCasesHex{CheckOrientationTestParameters{{Point{0.0, 0.0, 0.0}, Point{1.0, 0.0, 0.0},
                                                 Point{1.0, 1.0, 0.0}, Point{0.0, 1.0, 0.0},
                                                 Point{0.0, 0.0, 1.0}, Point{1.0, 0.0, 1.0},
                                                 Point{1.0, 1.0, 1.0}, Point{0.0, 1.0, 1.0}},
                                                true},
                 CheckOrientationTestParameters{{Point{0.0, 0.0, 0.0}, Point{1.0, 0.0, 0.0},
                                                 Point{1.0, 1.0, 0.0}, Point{0.0, 1.0, 0.0},
                                                 Point{0.0, 0.0, 1.0}, Point{1.0, 0.0, 1.0},
                                                 Point{1.0, 1.0, 1.0}, Point{0.0, 1.0, 1.0},
                                                 Point{0.5, 0.5, 0.5}},
                                                false},
                 CheckOrientationTestParameters{{Point{0.0, 0.0, 0.0}, Point{1.0, 0.0, 0.0},
                                                 Point{1.0, 1.0, 0.0},
                                                 Point{0.0, GeometricTolerance / 10, 0.0},
                                                 Point{0.0, 0.0, 1.0}, Point{1.0, 0.0, 1.0},
                                                 Point{1.0, 1.0, 1.0},
                                                 Point{0.0, GeometricTolerance / 10, 1.0}},
                                                false},
                 CheckOrientationTestParameters{{Point{0.0, 0.0, 0.0}, Point{1.0, 0.0, 0.0},
                                                 Point{0.5, 0.5, 0.0}, Point{0.0, 1.0, 0.0},
                                                 Point{0.0, 0.0, 1.0}, Point{1.0, 0.0, 1.0},
                                                 Point{0.5, 0.5, 1.0}, Point{0.0, 1.0, 1.0}},
                                                false},
                 CheckOrientationTestParameters{{Point{0.0, 0.0, 0.0}, Point{1.0, 0.0, 0.0},
                                                 Point{0.25, 0.25, 0.0}, Point{0.0, 1.0, 0.0},
                                                 Point{0.0, 0.0, 1.0}, Point{1.0, 0.0, 1.0},
                                                 Point{0.25, 0.25, 1.0}, Point{0.0, 1.0, 1.0}},
                                                false},
                 CheckOrientationTestParameters{{Point{180.0, -18.5, 0.0}, Point{240.0, -23.0, 0.0},
                                                 Point{240.0, 23.0, 0.0}, Point{180.0, 18.5, 0.0},
                                                 Point{180.0, -18.5, 100.0},
                                                 Point{240.0, -23.0, 100.0},
                                                 Point{240.0, 23.0, 100.0},
                                                 Point{180.0, 18.5, 100.0}},
                                                true},
                 CheckOrientationTestParameters{{Point{-8.0 / 63.0, 11.0 / 21.0, 8.0 / 9.0},
                                                 Point{73.0 / 37.0, -23.0 / 37, -66.0 / 37.0},
                                                 Point{-5.0 / 7.0, -2.0 / 7.0, -36.0 / 7.0},
                                                 Point{-24.0 / 13.0, 11.0 / 13.0, -8.0 / 13.0},
                                                 Point{38.0 / 17.0, 33.0 / 17.0, -13.0 / 17.0},
                                                 Point{118.0 / 37.0, 13.0 / 37.0, -93.0 / 37.0},
                                                 Point{4.0 / 7.0, 19.0 / 28.0, -81.0 / 14.0},
                                                 Point{4.0 / 7.0, 33.0 / 14.0, -17.0 / 7.0}},
                                                true},
                 CheckOrientationTestParameters{{Point{-8.0 / 63.0, 11.0 / 21.0, 8.0 / 9.0},
                                                 Point{73.0 / 37.0, -23.0 / 37, -66.0 / 37.0},
                                                 Point{-5.0 / 7.0, -2.0 / 7.0, -36.0 / 7.0},
                                                 Point{-24.0 / 13.0, 11.0 / 13.0, -8.0 / 13.0},
                                                 Point{38.0 / 17.0, 33.0 / 17.0, -13.0 / 17.0},
                                                 Point{118.0 / 37.0, 13.0 / 37.0, -93.0 / 37.0},
                                                 Point{4.0 / 7.0, 19.0 / 28.0, -81.0 / 14.0},
                                                 Point{4.0 / 7.0, 33.0 / 14.0, -18.0 / 7.0}},
                                                false}};

#endif // TESTCHECKORIENTATION_HPP
