#include "TestCell.hpp"
#include "TestEnvironment.hpp"

#include "CheckOrientation.hpp"

#include <utility>

int getCellDimension(CELLTYPE cType) {
  if (cType == INTERVAL) return 1;
  if (cType == TRIANGLE || cType == QUADRILATERAL) return 2;
  if (cType == TETRAHEDRON || cType == HEXAHEDRON) return 3;
  Exit("Cell Type unknown.")
}

class CellConstructionTest : public CellTest {
public:
  CellConstructionTest() : CellTest(getReferenceParameters(GetParam().celltype)) {}

  double FaceArea(vector<Point> facePoints) {
    switch (facePoints.size()) {
    case 1:
      return 1.0;
    case 2:
      return dist(facePoints[0], facePoints[1]);
    case 3:
      return norm((facePoints[1] - facePoints[0]) ^ (facePoints[2] - facePoints[0])) / 2.0;
    case 4:
      return norm((facePoints[1] - facePoints[0]) ^ (facePoints[3] - facePoints[0]));
    default:
      Exit("No default areay for polygons with n>4.")
    }
  }
};

/**
 * Checks wether a cell contains exactly the specified corners.
 * This implicitly tests that all corners are fetched with the
 *      Cell::Corner(int i)
 * method.
 */
TEST_P(CellConstructionTest, TestCorners) {
  int n = c->Corners();
  vector<Point> cellCorners(n);

  for (int i = 0; i < n; i++) {
    cellCorners[i] = c->Corner(i);
  }
  AssertContained(GetParam().corners, cellCorners);
}

/**
 * Checks wether a cell contains exactly the specified edges.
 * This implicitly tests that all edges are fetched with the
 *      Cell::Edge(int i)
 * method.
 */
TEST_P(CellConstructionTest, TestEdges) {
  int n = c->Edges();
  vector<Point> cellEdges(n);

  for (int i = 0; i < n; i++) {
    cellEdges[i] = c->Edge(i);
  }
  AssertContained(GetParam().edges, cellEdges);
}

/**
 * Checks wether a cell contains exactly the specified faces.
 * This implicitly tests that all faces are fetched with the
 *      Cell::Face(int i)
 * method.
 */
TEST_P(CellConstructionTest, TestFaces) {
  int n = c->Faces();
  vector<Point> cellFaces(n);

  for (int i = 0; i < n; i++) {
    cellFaces[i] = c->Face(i);
  }
  AssertContained(GetParam().faces, cellFaces);
}

TEST_P(CellConstructionTest, TestEdgeCorners) {
  int n = c->Edges();
  vector<Point> eCorners(2);
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < 2; j++) {
      EXPECT_EQ(c->EdgeCorner(i, j), c->Corner(c->edgecorner(i, j)));
      eCorners[j] = c->EdgeCorner(i, j);
    }
    EXPECT_EQ(c->Edge(i), midPoint(eCorners));
  }
}

TEST_P(CellConstructionTest, TestFaceCorners) {
  int n = c->Faces();
  vector<Point> fCorners(c->FaceCorners(0));
  for (int i = 0; i < n; i++) {
    fCorners.resize(c->FaceCorners(i));
    for (int j = 0; j < c->FaceCorners(i); j++) {
      EXPECT_EQ(c->FaceCorner(i, j), c->Corner(c->facecorner(i, j)));
      fCorners[j] = c->FaceCorner(i, j);
    }
    EXPECT_EQ(c->Face(i), midPoint(fCorners));
  }
}

TEST_P(CellConstructionTest, TestFaceArea) {
  vector<Point> fCorners(2);
  for (int f = 0; f < c->Faces(); f++) {
    fCorners.resize(c->FaceCorners(f));
    for (int fc = 0; fc < c->FaceCorners(f); fc++)
      fCorners[fc] = c->FaceCorner(f, fc);
    EXPECT_EQ(FaceArea(fCorners), c->FaceArea(f));
  }
}

/**
 * Checks that the corners, edges, faces and center of the reference cell are correctly
 * initialized and accessed.
 */
TEST_P(CellConstructionTest, TestReferenceCell) {
  /// Test Local Corners
  int n = c->Corners();
  vector<Point> testPoints(n);
  for (int i = 0; i < n; i++) {
    testPoints[i] = c->LocalCorner(i);
  }
  AssertContained(refCorners, testPoints);

  /// Test Local Edges
  n = c->Edges();
  testPoints.resize(n);
  for (int i = 0; i < n; i++) {
    testPoints[i] = c->LocalEdge(i);
  }
  AssertContained(refEdges, testPoints);

  /// Test Local Faces
  n = c->Faces();
  testPoints.resize(n);
  vector<Point> fCorners(2);
  for (int i = 0; i < n; i++) {
    testPoints[i] = c->LocalFace(i);

    /// Test Local Face Area
    fCorners.resize(c->FaceCorners(i));
    for (int fc = 0; fc < c->FaceCorners(i); fc++)
      fCorners[fc] = c->LocalCorner(c->facecorner(i, fc));
    EXPECT_EQ(FaceArea(fCorners), c->LocalFaceArea(i));
  }
  AssertContained(refFaces, testPoints);

  /// Test Local Center
  EXPECT_EQ(refCenter, c->LocalCenter());
}

TEST_P(CellConstructionTest, TestLocalToGlobalViceVersa) {
  // TODO: Fix GlobalToLocal in Hex
  if (c->ReferenceType() == CELLTYPE::HEXAHEDRON) { return; }
  for (int i = 0; i < c->Corners(); ++i) {
    std::cout << DOUT(c->Corner(i)) << endl;
    std::cout << DOUT(c->GlobalToLocal(c->Corner(i))) << endl;
    std::cout << DOUT(c->LocalToGlobal(c->GlobalToLocal(c->Corner(i)))) << endl;
    EXPECT_EQ(c->LocalToGlobal(c->GlobalToLocal(c->Corner(i))), c->Corner(i)) << c;
    EXPECT_EQ(c->GlobalToLocal(c->LocalToGlobal(c->LocalCorner(i))), c->LocalCorner(i)) << c;
  }
  for (int i = 0; i < c->Edges(); ++i) {
    EXPECT_EQ(c->LocalToGlobal(c->GlobalToLocal(c->Edge(i))), c->Edge(i));
    EXPECT_EQ(c->GlobalToLocal(c->LocalToGlobal(c->LocalEdge(i))), c->LocalEdge(i));
  }
  for (int i = 0; i < c->Faces(); ++i) {
    EXPECT_EQ(c->LocalToGlobal(c->GlobalToLocal(c->Face(i))), c->Face(i));
    EXPECT_EQ(c->GlobalToLocal(c->LocalToGlobal(c->LocalFace(i))), c->LocalFace(i));
  }
  EXPECT_EQ(c->LocalToGlobal(c->GlobalToLocal(c->Center())), c->Center());
  EXPECT_EQ(c->GlobalToLocal(c->LocalToGlobal(c->LocalCenter())), c->LocalCenter());
}

TEST_P(CellConstructionTest, TestLocalToGlobal) {
  for (int i = 0; i < c->Corners(); ++i) {
    EXPECT_EQ(c->LocalToGlobal(c->LocalCorner(i)), c->Corner(i));
  }
  for (int i = 0; i < c->Edges(); ++i) {
    EXPECT_EQ(c->LocalToGlobal(c->LocalEdge(i)), c->Edge(i));
  }
  for (int i = 0; i < c->Faces(); ++i) {
    EXPECT_EQ(c->LocalToGlobal(c->LocalFace(i)), c->Face(i));
  }
  EXPECT_EQ(c->LocalToGlobal(c->LocalCenter()), c->Center());
}

TEST_P(CellConstructionTest, TestGlobalToLocal) {
  // TODO: Fix GlobalToLocal in Quad and Hex
  if (c->ReferenceType() == CELLTYPE::HEXAHEDRON) { return; }
  for (int i = 0; i < c->Corners(); ++i) {
    // mout << c->Corner(i) << endl;
    EXPECT_EQ(c->GlobalToLocal(c->Corner(i)), c->LocalCorner(i));
  }
  for (int i = 0; i < c->Edges(); ++i) {
    EXPECT_EQ(c->GlobalToLocal(c->Edge(i)), c->LocalEdge(i));
  }
  for (int i = 0; i < c->Faces(); ++i) {
    EXPECT_EQ(c->GlobalToLocal(c->Face(i)), c->LocalFace(i));
  }
  EXPECT_EQ(c->GlobalToLocal(c->Center()), c->LocalCenter());
}

TEST_P(CellConstructionTest, TestFaceLocalToGlobal) {
  for (int face = 0; face < c->Faces(); ++face) {
    if (c->ReferenceType() == INTERVAL) {
      EXPECT_EQ(c->FaceLocalToGlobal(face, Point(1.0)), c->FaceCorner(face, 0));
    } else if (c->ReferenceType() == TRIANGLE || c->ReferenceType() == QUADRILATERAL) {
      EXPECT_EQ(c->FaceLocalToGlobal(face, Point(0.0)), c->FaceCorner(face, 0));
      EXPECT_EQ(c->FaceLocalToGlobal(face, Point(1.0)), c->FaceCorner(face, 1));
      EXPECT_EQ(c->FaceLocalToGlobal(face, Point(0.5)), c->Face(face));
    } else if (c->ReferenceType() == TETRAHEDRON) {
      EXPECT_EQ(c->FaceLocalToGlobal(face, Point(0.0, 0.0)), c->FaceCorner(face, 0));
      EXPECT_EQ(c->FaceLocalToGlobal(face, Point(1.0, 0.0)), c->FaceCorner(face, 1));
      EXPECT_EQ(c->FaceLocalToGlobal(face, Point(0.0, 1.0)), c->FaceCorner(face, 2));
      EXPECT_EQ(c->FaceLocalToGlobal(face, Point(0.5, 0.0)), c->FaceEdge(face, 0));
      EXPECT_EQ(c->FaceLocalToGlobal(face, Point(0.5, 0.5)), c->FaceEdge(face, 1));
      EXPECT_EQ(c->FaceLocalToGlobal(face, Point(0.0, 0.5)), c->FaceEdge(face, 2));
      EXPECT_EQ(c->FaceLocalToGlobal(face, Point(1.0 / 3.0, 1.0 / 3.0)), c->Face(face));
    } else if (c->ReferenceType() == HEXAHEDRON) {
      //      EXPECT_EQ(c->FaceLocalToGlobal(face, Point(0.0, 0.0)), c->FaceCorner(face, 0));
      //      EXPECT_EQ(c->FaceLocalToGlobal(face, Point(1.0, 0.0)), c->FaceCorner(face, 1));
      //      EXPECT_EQ(c->FaceLocalToGlobal(face, Point(1.0, 1.0)), c->FaceCorner(face, 2));
      //      EXPECT_EQ(c->FaceLocalToGlobal(face, Point(0.0, 1.0)), c->FaceCorner(face, 3));
      //      EXPECT_EQ(c->FaceLocalToGlobal(face, Point(0.5, 0.0)), c->FaceEdge(face, 0));
      //      EXPECT_EQ(c->FaceLocalToGlobal(face, Point(1.0, 0.5)), c->FaceEdge(face, 1));
      //      EXPECT_EQ(c->FaceLocalToGlobal(face, Point(0.5, 1.0)), c->FaceEdge(face, 2));
      //      EXPECT_EQ(c->FaceLocalToGlobal(face, Point(0.0, 0.5)), c->FaceEdge(face, 3));
      //      EXPECT_EQ(c->FaceLocalToGlobal(face, Point(0.25, 0.25)), c->Face(face));
    } else {
      THROW("TestFaceLocalToGlobal not implemented for cell type "
            + std::to_string(c->ReferenceType()))
    }
  }
}

TEST_P(CellConstructionTest, TestFaceLocalToLocal) {
  for (int face = 0; face < c->Faces(); ++face) {
    if (c->ReferenceType() == INTERVAL) {
      EXPECT_EQ(c->FaceLocalToLocal(face, Point(1.0)), c->LocalFaceCorner(face, 0));
    } else if (c->ReferenceType() == TRIANGLE || c->ReferenceType() == QUADRILATERAL) {
      EXPECT_EQ(c->FaceLocalToLocal(face, Point(0.0)), c->LocalFaceCorner(face, 0));
      EXPECT_EQ(c->FaceLocalToLocal(face, Point(1.0)), c->LocalFaceCorner(face, 1));
      EXPECT_EQ(c->FaceLocalToLocal(face, Point(0.5)), c->LocalFace(face));
    } else if (c->ReferenceType() == TETRAHEDRON) {
      EXPECT_EQ(c->FaceLocalToLocal(face, Point(0.0, 0.0)), c->LocalFaceCorner(face, 0));
      EXPECT_EQ(c->FaceLocalToLocal(face, Point(1.0, 0.0)), c->LocalFaceCorner(face, 1));
      EXPECT_EQ(c->FaceLocalToLocal(face, Point(0.0, 1.0)), c->LocalFaceCorner(face, 2));
      EXPECT_EQ(c->FaceLocalToLocal(face, Point(0.5, 0.0)), c->LocalFaceEdge(face, 0));
      EXPECT_EQ(c->FaceLocalToLocal(face, Point(0.5, 0.5)), c->LocalFaceEdge(face, 1));
      EXPECT_EQ(c->FaceLocalToLocal(face, Point(0.0, 0.5)), c->LocalFaceEdge(face, 2));
      EXPECT_EQ(c->FaceLocalToLocal(face, Point(1.0 / 3.0, 1.0 / 3.0)), c->LocalFace(face));
    } else if (c->ReferenceType() == HEXAHEDRON) {
      //      EXPECT_EQ(c->FaceLocalToLocal(face, Point(0.0, 0.0)), c->LocalFaceCorner(face, 0));
      //      EXPECT_EQ(c->FaceLocalToLocal(face, Point(1.0, 0.0)), c->LocalFaceCorner(face, 1));
      //      EXPECT_EQ(c->FaceLocalToLocal(face, Point(1.0, 1.0)), c->LocalFaceCorner(face, 2));
      //      EXPECT_EQ(c->FaceLocalToLocal(face, Point(0.0, 1.0)), c->LocalFaceCorner(face, 3));
      //      EXPECT_EQ(c->FaceLocalToLocal(face, Point(0.5, 0.0)), c->LocalFaceEdge(face, 0));
      //      EXPECT_EQ(c->FaceLocalToLocal(face, Point(1.0, 0.5)), c->LocalFaceEdge(face, 1));
      //      EXPECT_EQ(c->FaceLocalToLocal(face, Point(0.5, 1.0)), c->LocalFaceEdge(face, 2));
      //      EXPECT_EQ(c->FaceLocalToLocal(face, Point(0.0, 0.5)), c->LocalFaceEdge(face, 3));
      //      EXPECT_EQ(c->FaceLocalToLocal(face, Point(0.25, 0.25)), c->LocalFace(face));
    } else {
      THROW("TestFaceLocalToLocal not implemented for cell type "
            + std::to_string(c->ReferenceType()))
    }
  }
}

TEST_P(CellConstructionTest, TestFaceLocalToLocalANDGlobal) {
  for (int face = 0; face < c->Faces(); ++face) {
    if (c->ReferenceType() == INTERVAL) {
      EXPECT_EQ(c->LocalToGlobal(c->FaceLocalToLocal(face, Point(1.0))),
                c->FaceLocalToGlobal(face, Point(1.0)));
    } else if (c->ReferenceType() == TRIANGLE || c->ReferenceType() == QUADRILATERAL) {
      EXPECT_EQ(c->LocalToGlobal(c->FaceLocalToLocal(face, Point(0.0))),
                c->FaceLocalToGlobal(face, Point(0.0)));
      EXPECT_EQ(c->LocalToGlobal(c->FaceLocalToLocal(face, Point(0.5))),
                c->FaceLocalToGlobal(face, Point(0.5)));
      EXPECT_EQ(c->LocalToGlobal(c->FaceLocalToLocal(face, Point(1.0))),
                c->FaceLocalToGlobal(face, Point(1.0)));
    } else if (c->ReferenceType() == TETRAHEDRON) {
      EXPECT_EQ(c->LocalToGlobal(c->FaceLocalToLocal(face, Point(0.0, 0.0))),
                c->FaceLocalToGlobal(face, Point(0.0, 0.0)));
      EXPECT_EQ(c->LocalToGlobal(c->FaceLocalToLocal(face, Point(1.0, 0.0))),
                c->FaceLocalToGlobal(face, Point(1.0, 0.0)));
      EXPECT_EQ(c->LocalToGlobal(c->FaceLocalToLocal(face, Point(0.0, 1.0))),
                c->FaceLocalToGlobal(face, Point(0.0, 1.0)));
      EXPECT_EQ(c->LocalToGlobal(c->FaceLocalToLocal(face, Point(0.5, 0.0))),
                c->FaceLocalToGlobal(face, Point(0.5, 0.0)));
      EXPECT_EQ(c->LocalToGlobal(c->FaceLocalToLocal(face, Point(0.5, 0.5))),
                c->FaceLocalToGlobal(face, Point(0.5, 0.5)));
      EXPECT_EQ(c->LocalToGlobal(c->FaceLocalToLocal(face, Point(0.0, 0.5))),
                c->FaceLocalToGlobal(face, Point(0.0, 0.5)));
      EXPECT_EQ(c->LocalToGlobal(c->FaceLocalToLocal(face, Point(1.0 / 3.0, 1.0 / 3.0))),
                c->FaceLocalToGlobal(face, Point(1.0 / 3.0, 1.0 / 3.0)));
    } else if (c->ReferenceType() == HEXAHEDRON) {
      EXPECT_EQ(c->LocalToGlobal(c->FaceLocalToLocal(face, Point(0.0, 0.0))),
                c->FaceLocalToGlobal(face, Point(0.0, 0.0)));
      EXPECT_EQ(c->LocalToGlobal(c->FaceLocalToLocal(face, Point(1.0, 0.0))),
                c->FaceLocalToGlobal(face, Point(1.0, 0.0)));
      EXPECT_EQ(c->LocalToGlobal(c->FaceLocalToLocal(face, Point(1.0, 1.0))),
                c->FaceLocalToGlobal(face, Point(1.0, 1.0)));
      EXPECT_EQ(c->LocalToGlobal(c->FaceLocalToLocal(face, Point(0.0, 1.0))),
                c->FaceLocalToGlobal(face, Point(0.0, 1.0)));
      EXPECT_EQ(c->LocalToGlobal(c->FaceLocalToLocal(face, Point(0.5, 0.0))),
                c->FaceLocalToGlobal(face, Point(0.5, 0.0)));
      EXPECT_EQ(c->LocalToGlobal(c->FaceLocalToLocal(face, Point(1.0, 0.5))),
                c->FaceLocalToGlobal(face, Point(1.0, 0.5)));
      EXPECT_EQ(c->LocalToGlobal(c->FaceLocalToLocal(face, Point(0.5, 1.0))),
                c->FaceLocalToGlobal(face, Point(0.5, 1.0)));
      EXPECT_EQ(c->LocalToGlobal(c->FaceLocalToLocal(face, Point(0.0, 0.5))),
                c->FaceLocalToGlobal(face, Point(0.0, 0.5)));
      EXPECT_EQ(c->LocalToGlobal(c->FaceLocalToLocal(face, Point(0.25, 0.25))),
                c->FaceLocalToGlobal(face, Point(0.25, 0.25)));
    } else {
      THROW("TestFaceLocalToLocal not implemented for cell type "
            + std::to_string(c->ReferenceType()))
    }
  }
}

TEST_P(CellConstructionTest, TestDimension) {
  int expectedDim = getCellDimension(GetParam().celltype);
  EXPECT_EQ(expectedDim, c->dim());
  EXPECT_EQ(expectedDim < 3, c->plane());
}

TEST_P(CellConstructionTest, TestSubdomain) { EXPECT_EQ(sd, c->Subdomain()); }

TEST_P(CellConstructionTest, TestOrientation) {
  std::vector<Point> corners = c->copyCorners();
  std::vector<int> cornerIndices(corners.size());
  std::iota(cornerIndices.begin(), cornerIndices.end(), 0);
  CheckOrientation::Check(c->Type(), corners, cornerIndices);
  for (int i = 0; i < corners.size(); ++i) {
    EXPECT_EQ(cornerIndices[i], i);
    EXPECT_EQ(corners[i], c->Corner(i));
  }
}

TEST_P(CellConstructionTest, TestRefinementRules) {
  vector<Point> z;
  const vector<Rule> &rules = c->Refine(z);
  for (const Rule &rule : rules) {
    Cell *childCell = CreateCell(rule.type(), c->Subdomain(), rule.getCorners(z));
    std::vector<Point> corners = childCell->copyCorners();
    std::vector<int> cornerIndices(corners.size());
    std::iota(cornerIndices.begin(), cornerIndices.end(), 0);
    CheckOrientation::Check(childCell->Type(), corners, cornerIndices);
    for (int i = 0; i < corners.size(); ++i) {
      EXPECT_EQ(cornerIndices[i], i);
      EXPECT_EQ(corners[i], childCell->Corner(i));
    }
    delete childCell;
  }
}

INSTANTIATE_TEST_SUITE_P(IntervalConstructionTests, CellConstructionTest,
                         ValuesIn(intervalTestCases));

#if SpaceDimension >= 2
INSTANTIATE_TEST_SUITE_P(TriangleConstructionTests, CellConstructionTest,
                         ValuesIn(triangleTestCases));

INSTANTIATE_TEST_SUITE_P(QuadConstructionTests, CellConstructionTest, ValuesIn(quadTestCases));
#endif
#if SpaceDimension >= 3
INSTANTIATE_TEST_SUITE_P(TetrahedronConstructionTests, CellConstructionTest,
                         ValuesIn(tetrahedronTestCases));

INSTANTIATE_TEST_SUITE_P(HexahedronConstructionTests, CellConstructionTest,
                         ValuesIn(hexahedronTestCases));
#endif

int main(int argc, char **argv) {
  constructTestData(cellTestData);
  MppTest mppTest = MppTestBuilder(argc, argv).WithoutDefaultConfig().WithPPM().WithScreenLogging();
  return mppTest.RUN_ALL_MPP_TESTS();
}