#ifndef TESTRTSHAPE_HPP
#define TESTRTSHAPE_HPP

#include "Point.hpp"

struct ShapeTestParameters {
  int order;
  int shapeSize;
  std::vector<Point> expectedNodalPoints;
};

static std::vector<ShapeTestParameters> TestCasesInt{
    ShapeTestParameters{0, 2, {Point(0.0), Point(1.0)}}};

static std::vector<ShapeTestParameters> TestCasesQuad{
    ShapeTestParameters{0,
                        4,
                        {Point(1.0 / 2.0, 0.0), Point(1.0, 1.0 / 2.0), Point(1.0 / 2.0, 1.0),
                         Point(0.0, 1.0 / 2.0)}}};

static std::vector<ShapeTestParameters>
    TestCasesTri{ShapeTestParameters{0,
                                     3,
                                     {Point(1.0 / 2.0, 0.0), Point(1.0 / 2.0, 1.0 / 2.0),
                                      Point(0.0, 1.0 / 2.0)}},
                 ShapeTestParameters{1,
                                     8,
                                     {Point(1.0 / 3.0, 0.0), Point(2.0 / 3.0, 0.0),
                                      Point(2.0 / 3.0, 1.0 / 3.0), Point(1.0 / 3.0, 2.0 / 3.0),
                                      Point(0.0, 2.0 / 3.0), Point(0.0, 1.0 / 3.0),
                                      Point(1.0 / 3.0, 1.0 / 3.0), Point()}},
                 ShapeTestParameters{2,
                                     15,
                                     {Point(1.0 / 4.0, 0.0), Point(1.0 / 2.0, 0.0),
                                      Point(3.0 / 4.0, 0.0), Point(0.0, 1.0 / 4.0),
                                      Point(1.0 / 4.0, 1.0 / 4.0), Point(1.0 / 2.0, 1.0 / 4.0),
                                      Point(3.0 / 4.0, 1.0 / 4.0), Point(0.0, 1.0 / 2.0),
                                      Point(1.0 / 4.0, 1.0 / 2.0), Point(1.0 / 2.0, 1.0 / 2.0),
                                      Point(0.0, 3.0 / 4.0), Point(1.0 / 4.0, 3.0 / 4.0), Point(),
                                      Point(), Point()}},
                 ShapeTestParameters{3,
                                     24,
                                     {Point(0.0 / 5.0, 1.0 / 5.0),
                                      Point(0.0 / 5.0, 2.0 / 5.0),
                                      Point(0.0 / 5.0, 3.0 / 5.0),
                                      Point(0.0 / 5.0, 4.0 / 5.0),
                                      Point(1.0 / 5.0, 0.0 / 5.0),
                                      Point(1.0 / 5.0, 1.0 / 5.0),
                                      Point(1.0 / 5.0, 2.0 / 5.0),
                                      Point(1.0 / 5.0, 3.0 / 5.0),
                                      Point(1.0 / 5.0, 4.0 / 5.0),
                                      Point(2.0 / 5.0, 0.0 / 5.0),
                                      Point(2.0 / 5.0, 1.0 / 5.0),
                                      Point(2.0 / 5.0, 2.0 / 5.0),
                                      Point(2.0 / 5.0, 3.0 / 5.0),
                                      Point(3.0 / 5.0, 0.0 / 5.0),
                                      Point(3.0 / 5.0, 1.0 / 5.0),
                                      Point(3.0 / 5.0, 2.0 / 5.0),
                                      Point(4.0 / 5.0, 0.0 / 5.0),
                                      Point(4.0 / 5.0, 1.0 / 5.0),
                                      Point(),
                                      Point(),
                                      Point(),
                                      Point(),
                                      Point(),
                                      Point()}}};

static std::vector<ShapeTestParameters> TestCasesHex{
    ShapeTestParameters{0,
                        6,
                        {Point(1.0 / 2.0, 1.0 / 2.0, 0.0), Point(1.0 / 2.0, 0.0, 1.0 / 2.0),
                         Point(1.0, 1.0 / 2.0, 1.0 / 2.0), Point(1.0 / 2.0, 1.0, 1.0 / 2.0),
                         Point(0.0, 1.0 / 2.0, 1.0 / 2.0), Point(1.0 / 2.0, 1.0 / 2.0, 1.0)}}};

static std::vector<ShapeTestParameters> TestCasesTet{
    ShapeTestParameters{0,
                        4,
                        {Point(1.0 / 3.0, 1.0 / 3.0, 0.0), Point(1.0 / 3.0, 0.0, 1.0 / 3.0),
                         Point(0.0, 1.0 / 3.0, 1.0 / 3.0),
                         Point(1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0)}}};

#endif // TESTRTSHAPE_HPP
