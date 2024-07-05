#ifndef TESTCELL_H
#define TESTCELL_H
#include <fstream>
#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include "TestCellData.hpp"
#include "cells/Cell.hpp"


using namespace ::testing;

class CellTest : public TestWithParam<CellTestParameters> {
protected:
  CELLTYPE cType;
  short sd;
  Cell *c = nullptr;
  const vector<Point> &refCorners = {};
  const vector<Point> &refEdges = {};
  const vector<Point> &refFaces = {};
  const Point &refCenter;

  explicit CellTest(const CellTestParameters &refParam) :
      cType(refParam.celltype), refCorners(refParam.corners), refEdges(refParam.edges),
      refFaces(refParam.faces), refCenter(refParam.center) {}

  void SetUp() override {
    if (!GetParam().isDegenerated) constructCell();
  }

  ~CellTest() override { delete c; }

  void AssertContained(const vector<Point> &referencePoints, vector<Point> &testPoints) {
    int n = referencePoints.size();
    ASSERT_EQ(n, testPoints.size());
    vector<bool> contained(n, false);
    for (const auto &p : testPoints) {
      // Ensures that the corners of testCell are in corners and that each element of corners
      // is actually a corner of testCell
      for (int j = 0; j < n; j++) {
        contained[j] = (contained[j] || p == referencePoints[j]);
      }
    }
    ASSERT_THAT(contained, Each(true));
  }

  Point midPoint(vector<Point> points) {
    Point mp(0.0, 0.0, 0.0);
    for (auto p : points)
      mp += p;
    return mp / points.size();
  }

  void constructCell() {
    sd = rand() % static_cast<int>(1000);
    c = CreateCell(cType, sd, GetParam().corners);
  }
};

inline CELLTYPE getCellType(const string &typeName) {
  if (typeName == "Interval") return INTERVAL;
  if (typeName == "Triangle") return TRIANGLE;
  if (typeName == "Quadrilateral") return QUADRILATERAL;
  if (typeName == "Tetrahedron") return TETRAHEDRON;
  if (typeName == "Hexahedron") return HEXAHEDRON;
  else Exit("Celltype Name unknown.");
}

static vector<CellTestParameters> intervalTestCases;

static vector<CellTestParameters> triangleTestCases;

static vector<CellTestParameters> quadTestCases;

static vector<CellTestParameters> tetrahedronTestCases;

static vector<CellTestParameters> hexahedronTestCases;

void constructTestData(const std::vector<CellTestParameters> &data) {
  for (auto &cellData : data) {
    switch (cellData.celltype) {
    case INTERVAL:
      intervalTestCases.push_back(cellData);
      break;
#if SpaceDimension >= 2
    case TRIANGLE:
      triangleTestCases.push_back(cellData);
      break;
    case QUADRILATERAL:
      quadTestCases.push_back(cellData);
      break;
#endif
#if SpaceDimension >= 3
    case TETRAHEDRON:
      tetrahedronTestCases.push_back(cellData);
      break;
    case HEXAHEDRON:
      hexahedronTestCases.push_back(cellData);
      break;
#endif
    default:
      continue;
    }
  }
}

void readTestData(const string &fileName) {
  std::ifstream dataFile("../../tests/cells/" + fileName);
  if (dataFile.is_open()) {
    string line;
    string value;
    int c = 0, e = 0, f = 0;
    while (getline(dataFile, line)) {
      std::istringstream linestream(line);
      // First entry is cell type.
      getline(linestream, value, ' ');
      CELLTYPE tp = getCellType(value);
      // Second entry is degenration of cell.
      getline(linestream, value, ' ');
      bool isDegen = (std::stod(value) > 0);

      switch (tp) {
      case INTERVAL:
        c = 2;
        e = 2;
        f = 2;
        break;
      case TRIANGLE:
        c = 3;
        e = 3;
        f = 3;
        break;
      case QUADRILATERAL:
        c = 4;
        e = 4;
        f = 4;
        break;
      case TETRAHEDRON:
        c = 4;
        e = 6;
        f = 4;
        break;
      case HEXAHEDRON:
        c = 8;
        e = 12;
        f = 6;
        break;
      default:
        continue;
      }
      int m = max(c, e);
      m = max(m, f);
      double pointsInLine[3 * (c + e + f + 1)];
      for (int j = 0; getline(linestream, value, ' '); j++) {
        pointsInLine[j] = std::stod(value);
      }
      vector<Point> corners(c);
      vector<Point> edges(e);
      vector<Point> faces(f);
      for (int i = 0; i < m; i++) {
        if (i < c)
          corners[i] = Point(pointsInLine[3 * i], pointsInLine[3 * i + 1], pointsInLine[3 * i + 2]);
        if (i < e)
          edges[i] = Point(pointsInLine[3 * (c + i)], pointsInLine[3 * (c + i) + 1],
                           pointsInLine[3 * (c + i) + 2]);
        if (i < f)
          faces[i] = Point(pointsInLine[3 * (c + e + i)], pointsInLine[3 * (c + e + i) + 1],
                           pointsInLine[3 * (c + e + i) + 2]);
      }
      Point center(pointsInLine[3 * (c + e + f)], pointsInLine[3 * (c + e + f) + 1],
                   pointsInLine[3 * (c + e + f) + 2]);

      switch (tp) {
      case INTERVAL:
        intervalTestCases.push_back({tp, isDegen, corners, edges, faces, center});
        break;
      case TRIANGLE:
        triangleTestCases.push_back({tp, isDegen, corners, edges, faces, center});
        break;
      case QUADRILATERAL:
        quadTestCases.push_back({tp, isDegen, corners, edges, faces, center});
        break;
      case TETRAHEDRON:
        tetrahedronTestCases.push_back({tp, isDegen, corners, edges, faces, center});
        break;
      case HEXAHEDRON:
        hexahedronTestCases.push_back({tp, isDegen, corners, edges, faces, center});
        break;
      default:
        continue;
      }
    }
    dataFile.close();
  } else {
    Exit("No test data found.")
  }
}

static const vector<CellTestParameters>
    referenceCells{{INTERVAL,
                    false,
                    {Point(0.0, 0.0, 0.0), Point(1.0, 0.0, 0.0)},
                    {Point(0.0, 0.0, 0.0), Point(1.0, 0.0, 0.0)},
                    {Point(0.0, 0.0, 0.0), Point(1.0, 0.0, 0.0)},
                    Point(0.5, 0.0, 0.0)},
                   {TRIANGLE,
                    false,
                    {Point(0.0, 0.0, 0.0), Point(1.0, 0.0, 0.0), Point(0.0, 1.0, 0.0)},
                    {Point(0.5, 0.0, 0.0), Point(0.5, 0.5, 0.0), Point(0.0, 0.5, 0.0)},
                    {Point(0.5, 0.0, 0.0), Point(0.5, 0.5, 0.0), Point(0.0, 0.5, 0.0)},
                    Point(1. / 3., 1. / 3., 0.0)},
                   {QUADRILATERAL,
                    false,
                    {Point(0.0, 0.0, 0.0), Point(1.0, 0.0, 0.0), Point(1.0, 1.0, 0.0),
                     Point(0.0, 1.0, 0.0)},
                    {Point(0.5, 0.0, 0.0), Point(1, 0.5, 0.0), Point(0.5, 1.0, 0.0),
                     Point(0.0, 0.5, 0.0)},
                    {Point(0.5, 0.0, 0.0), Point(1, 0.5, 0.0), Point(0.5, 1.0, 0.0),
                     Point(0.0, 0.5, 0.0)},
                    Point(0.5, 0.5, 0.0)}
#if SpaceDimension >= 3
                   ,
                   {TETRAHEDRON,
                    false,
                    {Point(0.0, 0.0, 0.0), Point(1.0, 0.0, 0.0), Point(0.0, 1.0, 0.0),
                     Point(0.0, 0.0, 1.0)},
                    {Point(0.5, 0.0, 0.0), Point(0.5, 0.5, 0.0), Point(0.0, 0.5, 0.0),
                     Point(0.0, 0.0, 0.5), Point(0.5, 0.0, 0.5), Point(0.0, 0.5, 0.5)},
                    {Point(1. / 3., 1. / 3., 0.0), Point(1. / 3., 0.0, 1. / 3.),
                     Point(0.0, 1. / 3., 1. / 3.), Point(1. / 3., 1. / 3., 1. / 3.)},
                    Point(0.25, 0.25, 0.25)},
                   {HEXAHEDRON,
                    false,
                    // Corners
                    {Point(0.0, 0.0, 0.0), Point(1.0, 0.0, 0.0), Point(1.0, 1.0, 0.0),
                     Point(0.0, 1.0, 0.0), Point(0.0, 0.0, 1.0), Point(1.0, 0.0, 1.0),
                     Point(1.0, 1.0, 1.0), Point(0.0, 1.0, 1.0)},
                    // Edges
                    {Point(0.5, 0.0, 0.0), Point(1, 0.5, 0.0), Point(0.5, 1.0, 0.0),
                     Point(0.0, 0.5, 0.0), Point(1.0, 0.0, 0.5), Point(1, 0.0, 0.5),
                     Point(1.0, 1.0, 0.5), Point(0.0, 1.0, 0.5), Point(0.5, 0.0, 1.0),
                     Point(1, 0.5, 1.0), Point(0.5, 1.0, 1.0), Point(0.0, 0.5, 1.0)},
                    // Faces
                    {Point(0.5, 0.5, 0.0), Point(0.5, 0.0, 0.5), Point(1.0, 0.5, 0.5),
                     Point(0.5, 1.0, 0.5), Point(0.0, 0.5, 0.5), Point(0.5, 0.5, 1.0)},
                    // Center
                    Point(0.5, 0.5, 0.5)}
#endif
    };

inline const CellTestParameters &getReferenceParameters(CELLTYPE cType) {
  switch (cType) {
  case INTERVAL:
    return referenceCells[0];
  case TRIANGLE:
    return referenceCells[1];
  case QUADRILATERAL:
    return referenceCells[2];
  case TETRAHEDRON:
    return referenceCells[3];
  case HEXAHEDRON:
    return referenceCells[4];
  default:
    Exit("Reference Cell Type unknown.")
  }
}

#endif // TESTCELL_H
