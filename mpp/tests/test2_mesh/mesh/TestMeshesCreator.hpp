#ifndef TESTMESHESCREATOR_HPP
#define TESTMESHESCREATOR_HPP

#include "MeshesCreator.hpp"

#include "TestEnvironment.hpp"

class TestMeshesCreator : public TestWithParam<int> {
protected:
  int verbose = 0;

  int commSplit;

  int pLevel = 0;

  int level;

  Meshes *meshes;

  const Mesh &mesh;

  double maxMeshWidth;

  double minMeshWidth;

  TestMeshesCreator(const std::string &meshName, int commSplit = 0) :
      commSplit(commSplit), level(GetParam()),
      meshes(
          MeshesCreator(meshName).WithoutDistribute().WithPLevel(pLevel).WithLevel(level).Create()),
      mesh((*meshes)[meshes->FineLevel().WithCommSplit(commSplit)]),
      maxMeshWidth(mesh.MaxMeshWidth()), minMeshWidth(mesh.MinMeshWidth()) {
    Config::Get("TestMeshesCreatorVerbose", verbose);
  }

  int EulerCharacteristic() {
    if (mesh.dim() == 1) {
      return 1 - mesh.CellCountGeometry() + mesh.VertexCountGeometry();
    } else
      return 1 + mesh.CellCountGeometry() - mesh.EdgeCountGeometry() + mesh.VertexCountGeometry();
  }

  void TearDown() override {
    PPM->Barrier(0);
    delete meshes;
  }
};

/*
 * Interval
 */

class TestMeshesCreatorInterval : public TestMeshesCreator {
protected:
  explicit TestMeshesCreatorInterval(int commSplit = 0) :
      TestMeshesCreator("Interval", commSplit) {}
};

class TestMeshesCreatorIntervalWithSplit : public TestMeshesCreatorInterval {
protected:
  TestMeshesCreatorIntervalWithSplit() : TestMeshesCreatorInterval(1) {}
};

class TestMeshesCreatorIntervalWithDoubleSplit : public TestMeshesCreatorInterval {
protected:
  TestMeshesCreatorIntervalWithDoubleSplit() : TestMeshesCreatorInterval(2) {}
};

#define MESHESCREATOR_INTERVAL_TESTS(MeshesCreatorIntervalTestClass)                               \
                                                                                                   \
  INSTANTIATE_TEST_SUITE_P(TestMeshesCreator, MeshesCreatorIntervalTestClass, Values(0, 1, 2, 3)); \
                                                                                                   \
  TEST_P(MeshesCreatorIntervalTestClass, TestCellCountGeometry) {                                  \
    EXPECT_EQ(mesh.CellCountGeometry(), pow(2, level));                                            \
  }                                                                                                \
                                                                                                   \
  TEST_P(MeshesCreatorIntervalTestClass, TestIntervalDim) { EXPECT_EQ(mesh.dim(), 1); }            \
                                                                                                   \
  TEST_P(MeshesCreatorIntervalTestClass, TestVertexCountGeometry) {                                \
    EXPECT_EQ(mesh.VertexCountGeometry(), pow(2, level) + 1);                                      \
  }                                                                                                \
                                                                                                   \
  TEST_P(MeshesCreatorIntervalTestClass, TestEdgeCountGeometry) {                                  \
    EXPECT_EQ(mesh.EdgeCountGeometry(), pow(2, level) + 1);                                        \
  }                                                                                                \
                                                                                                   \
  TEST_P(MeshesCreatorIntervalTestClass, TestFaceCountGeometry) {                                  \
    EXPECT_EQ(mesh.FaceCountGeometry(), pow(2, level) + 1);                                        \
  }                                                                                                \
                                                                                                   \
  TEST_P(MeshesCreatorIntervalTestClass, TestMeshWidthMin) {                                       \
    EXPECT_DOUBLE_EQ(minMeshWidth, pow(2, -level));                                                \
  }                                                                                                \
                                                                                                   \
  TEST_P(MeshesCreatorIntervalTestClass, TestMeshWidthMax) {                                       \
    EXPECT_DOUBLE_EQ(maxMeshWidth, pow(2, -level));                                                \
  }                                                                                                \
                                                                                                   \
  TEST_P(MeshesCreatorIntervalTestClass, TestEulerCharacteristic) {                                \
    EXPECT_EQ(EulerCharacteristic(), 2);                                                           \
  }

/*
 * Triangle
 */

class TestMeshesCreatorTriangle : public TestMeshesCreator {
protected:
  explicit TestMeshesCreatorTriangle(int commSplit = 0) :
      TestMeshesCreator("Triangle", commSplit) {}
};

class TestMeshesCreatorTriangleWithSplit : public TestMeshesCreatorTriangle {
protected:
  TestMeshesCreatorTriangleWithSplit() : TestMeshesCreatorTriangle(1) {}
};

class TestMeshesCreatorTriangleWithDoubleSplit : public TestMeshesCreatorTriangle {
protected:
  TestMeshesCreatorTriangleWithDoubleSplit() : TestMeshesCreatorTriangle(2) {}
};

#define MESHESCREATOR_TRIANGLE_TESTS(MeshesCreatorTriangleTestClass)                               \
  ;                                                                                                \
                                                                                                   \
  INSTANTIATE_TEST_SUITE_P(TestMeshesCreator, MeshesCreatorTriangleTestClass, Values(0, 1, 2, 3)); \
                                                                                                   \
  TEST_P(MeshesCreatorTriangleTestClass, TestDim) { EXPECT_EQ(mesh.dim(), 2); }                    \
                                                                                                   \
  TEST_P(MeshesCreatorTriangleTestClass, TestCellCountGeometry) {                                  \
    EXPECT_EQ(mesh.CellCountGeometry(), pow(4, level));                                            \
  }                                                                                                \
                                                                                                   \
  TEST_P(MeshesCreatorTriangleTestClass,                                                           \
         TestVertexCountGeometry){VerboseWarning("Not implemented", 1)}                            \
                                                                                                   \
  TEST_P(MeshesCreatorTriangleTestClass,                                                           \
         TestEdgeCountGeometry){VerboseWarning("Not implemented", 1)}                              \
                                                                                                   \
  TEST_P(MeshesCreatorTriangleTestClass,                                                           \
         TestFaceCountGeometry){VerboseWarning("Not implemented", 1)}                              \
                                                                                                   \
  TEST_P(MeshesCreatorTriangleTestClass, TestMeshWidthMin) {                                       \
    EXPECT_EQ(minMeshWidth, pow(2, -level));                                                       \
  }                                                                                                \
                                                                                                   \
  TEST_P(MeshesCreatorTriangleTestClass, TestMeshWidthMax) {                                       \
    EXPECT_EQ(maxMeshWidth, sqrt(2) * pow(2, -level));                                             \
  }                                                                                                \
                                                                                                   \
  TEST_P(MeshesCreatorTriangleTestClass, TestEulerCharacteristic) {                                \
    EXPECT_EQ(EulerCharacteristic(), 2);                                                           \
  }

/*
 * Square
 */

class TestMeshesCreatorSquare : public TestMeshesCreator {
protected:
  explicit TestMeshesCreatorSquare(int commSplit = 0) : TestMeshesCreator("Square", commSplit) {}
};

class TestMeshesCreatorSquareWithSplit : public TestMeshesCreatorSquare {
protected:
  TestMeshesCreatorSquareWithSplit() : TestMeshesCreatorSquare(1) {}
};

class TestMeshesCreatorSquareWithDoubleSplit : public TestMeshesCreatorSquare {
protected:
  TestMeshesCreatorSquareWithDoubleSplit() : TestMeshesCreatorSquare(2) {}
};

#define MESHESCREATOR_SQUARE_TESTS(MeshesCreatorSquareTestClass)                                   \
  ;                                                                                                \
                                                                                                   \
  INSTANTIATE_TEST_SUITE_P(TestMeshesCreator, MeshesCreatorSquareTestClass, Values(0, 1, 2, 3));   \
                                                                                                   \
  TEST_P(MeshesCreatorSquareTestClass, TestDim) { EXPECT_EQ(mesh.dim(), 2); }                      \
                                                                                                   \
  TEST_P(MeshesCreatorSquareTestClass, TestCellCountGeometry) {                                    \
    EXPECT_EQ(mesh.CellCountGeometry(), pow(4, level));                                            \
  }                                                                                                \
                                                                                                   \
  TEST_P(MeshesCreatorSquareTestClass, TestVertexCountGeometry) {                                  \
    EXPECT_EQ(mesh.VertexCountGeometry(), pow((pow(2, level) + 1), 2));                            \
  }                                                                                                \
                                                                                                   \
  TEST_P(MeshesCreatorSquareTestClass, TestEdgeCountGeometry) {                                    \
    EXPECT_EQ(mesh.EdgeCountGeometry(),                                                            \
              4 + 3 * 2 * (pow(2, level) - 1) + 2 * pow((pow(2, level) - 1), 2));                  \
  }                                                                                                \
                                                                                                   \
  TEST_P(MeshesCreatorSquareTestClass, TestFaceCountGeometry) {                                    \
    EXPECT_EQ(mesh.FaceCountGeometry(),                                                            \
              4 + 3 * 2 * (pow(2, level) - 1) + 2 * pow((pow(2, level) - 1), 2));                  \
  }                                                                                                \
                                                                                                   \
  TEST_P(MeshesCreatorSquareTestClass, TestMeshWidthMin) {                                         \
    EXPECT_EQ(minMeshWidth, pow(2, -level));                                                       \
  }                                                                                                \
                                                                                                   \
  TEST_P(MeshesCreatorSquareTestClass, TestMeshWidthMax) {                                         \
    EXPECT_EQ(maxMeshWidth, pow(2, -level));                                                       \
  }                                                                                                \
                                                                                                   \
  TEST_P(MeshesCreatorSquareTestClass, TestEulerCharacteristic) {                                  \
    EXPECT_EQ(EulerCharacteristic(), 2);                                                           \
  }

/*
 * SquareTriangle
 */

class TestMeshesCreatorSquareTriangle : public TestMeshesCreator {
protected:
  explicit TestMeshesCreatorSquareTriangle(int commSplit = 0) :
      TestMeshesCreator("SquareTriangle", commSplit) {}
};

class TestMeshesCreatorSquareTriangleWithSplit : public TestMeshesCreatorSquareTriangle {
protected:
  TestMeshesCreatorSquareTriangleWithSplit() : TestMeshesCreatorSquareTriangle(1) {}
};

class TestMeshesCreatorSquareTriangleWithDoubleSplit : public TestMeshesCreatorSquareTriangle {
protected:
  TestMeshesCreatorSquareTriangleWithDoubleSplit() : TestMeshesCreatorSquareTriangle(2) {}
};

#define MESHESCREATOR_SQUARETRIANGLE_TESTS(MeshesCreatorSquareTriangleTestClass) ;

INSTANTIATE_TEST_SUITE_P(TestMeshesCreator, TestMeshesCreatorSquareTriangle, Values(0, 1, 2, 3));

TEST_P(TestMeshesCreatorSquareTriangle, TestDim) { EXPECT_EQ(mesh.dim(), 2); }

TEST_P(TestMeshesCreatorSquareTriangle, TestCellCountGeometry) {
  EXPECT_EQ(mesh.CellCountGeometry(), 2 * pow(4, level));
}

TEST_P(TestMeshesCreatorSquareTriangle,
       TestVertexCountGeometry){VerboseWarning("Not implemented", 1)}

TEST_P(TestMeshesCreatorSquareTriangle, TestEdgeCountGeometry){VerboseWarning("Not implemented", 1)}

TEST_P(TestMeshesCreatorSquareTriangle, TestFaceCountGeometry){VerboseWarning("Not implemented", 1)}

TEST_P(TestMeshesCreatorSquareTriangle, TestMeshWidthMin) {
  EXPECT_EQ(minMeshWidth, pow(2, -level));
}

TEST_P(TestMeshesCreatorSquareTriangle, TestMeshWidthMax) {
  EXPECT_EQ(maxMeshWidth, sqrt(2) * pow(2, -level));
}

/*
 * Tetrahedron
 */

class TestMeshesCreatorTetrahedron : public TestMeshesCreator {
protected:
  explicit TestMeshesCreatorTetrahedron(int commSplit = 0) :
      TestMeshesCreator("Tetrahedron", commSplit) {}
};

class TestMeshesCreatorTetrahedronWithSplit : public TestMeshesCreatorTetrahedron {
protected:
  TestMeshesCreatorTetrahedronWithSplit() : TestMeshesCreatorTetrahedron(1) {}
};

class TestMeshesCreatorTetrahedronWithDoubleSplit : public TestMeshesCreatorTetrahedron {
protected:
  TestMeshesCreatorTetrahedronWithDoubleSplit() : TestMeshesCreatorTetrahedron(2) {}
};

#define MESHESCREATOR_TETRAHEDRON_TESTS(MeshesCreatorTetrahedronTestClass)                         \
  ;                                                                                                \
                                                                                                   \
  INSTANTIATE_TEST_SUITE_P(TestMeshesCreator, MeshesCreatorTetrahedronTestClass,                   \
                           Values(0, 1, 2, 3));                                                    \
                                                                                                   \
  TEST_P(MeshesCreatorTetrahedronTestClass, TestDim) { EXPECT_EQ(mesh.dim(), 3); }                 \
                                                                                                   \
  TEST_P(MeshesCreatorTetrahedronTestClass, TestCellCountGeometry) {                               \
    EXPECT_EQ(mesh.CellCountGeometry(), pow(8, level));                                            \
  }                                                                                                \
                                                                                                   \
  TEST_P(MeshesCreatorTetrahedronTestClass,                                                        \
         TestVertexCountGeometry){VerboseWarning("Not implemented", 1)}                            \
                                                                                                   \
  TEST_P(MeshesCreatorTetrahedronTestClass,                                                        \
         TestEdgeCountGeometry){VerboseWarning("Not implemented", 1)}                              \
                                                                                                   \
  TEST_P(MeshesCreatorTetrahedronTestClass,                                                        \
         TestFaceCountGeometry){VerboseWarning("Not implemented", 1)}                              \
                                                                                                   \
  TEST_P(MeshesCreatorTetrahedronTestClass, TestMeshWidthMin) {                                    \
    EXPECT_EQ(minMeshWidth, pow(2, -level));                                                       \
  }                                                                                                \
                                                                                                   \
  TEST_P(MeshesCreatorTetrahedronTestClass, TestMeshWidthMax) {                                    \
    VerboseWarning("Not implemented", 1)                                                           \
  }

/*
 * Hexahedron
 */

class TestMeshesCreatorHexahedron : public TestMeshesCreator {
protected:
  explicit TestMeshesCreatorHexahedron(int commSplit = 0) :
      TestMeshesCreator("Hexahedron", commSplit) {}
};

class TestMeshesCreatorHexahedronWithSplit : public TestMeshesCreatorHexahedron {
protected:
  TestMeshesCreatorHexahedronWithSplit() : TestMeshesCreatorHexahedron(1) {}
};

class TestMeshesCreatorHexahedronWithDoubleSplit : public TestMeshesCreatorHexahedron {
protected:
  TestMeshesCreatorHexahedronWithDoubleSplit() : TestMeshesCreatorHexahedron(2) {}
};

#define MESHESCREATOR_HEXAHEDRON_TESTS(MeshesCreatorHexahedronTestClass)                           \
  ;                                                                                                \
                                                                                                   \
  INSTANTIATE_TEST_SUITE_P(TestMeshesCreator, MeshesCreatorHexahedronTestClass,                    \
                           Values(0, 1, 2, 3));                                                    \
                                                                                                   \
  TEST_P(MeshesCreatorHexahedronTestClass, TestDim) { EXPECT_EQ(mesh.dim(), 3); }                  \
                                                                                                   \
  TEST_P(MeshesCreatorHexahedronTestClass, TestCellCountGeometry) {                                \
    EXPECT_EQ(mesh.CellCountGeometry(), pow(8, level));                                            \
  }                                                                                                \
                                                                                                   \
  TEST_P(MeshesCreatorHexahedronTestClass,                                                         \
         TestVertexCountGeometry){VerboseWarning("Not implemented", 1)}                            \
                                                                                                   \
  TEST_P(MeshesCreatorHexahedronTestClass,                                                         \
         TestEdgeCountGeometry){VerboseWarning("Not implemented", 1)}                              \
                                                                                                   \
  TEST_P(MeshesCreatorHexahedronTestClass,                                                         \
         TestFaceCountGeometry){VerboseWarning("Not implemented", 1)}                              \
                                                                                                   \
  TEST_P(MeshesCreatorHexahedronTestClass, TestMeshWidthMin) {                                     \
    EXPECT_EQ(minMeshWidth, pow(2, -level));                                                       \
  }                                                                                                \
                                                                                                   \
  TEST_P(MeshesCreatorHexahedronTestClass, TestMeshWidthMax) {                                     \
    EXPECT_EQ(maxMeshWidth, pow(2, -level));                                                       \
  }

/*
 * SpaceTimeSquare
 */

class TestMeshesCreatorSpaceTimeSquare : public TestMeshesCreator {
protected:
  explicit TestMeshesCreatorSpaceTimeSquare(int commSplit = 0) :
      TestMeshesCreator("SpaceTimeSquare", commSplit) {}
};

class TestMeshesCreatorSpaceTimeSquareWithSplit : public TestMeshesCreatorSpaceTimeSquare {
protected:
  TestMeshesCreatorSpaceTimeSquareWithSplit() : TestMeshesCreatorSpaceTimeSquare(1) {}
};

class TestMeshesCreatorSpaceTimeSquareWithDoubleSplit : public TestMeshesCreatorSpaceTimeSquare {
protected:
  TestMeshesCreatorSpaceTimeSquareWithDoubleSplit() : TestMeshesCreatorSpaceTimeSquare(2) {}
};

#define MESHESCREATOR_SPACETIMESQUARE_TESTS(MeshesCreatorSpaceTimeSquareTestClass)                 \
  ;                                                                                                \
                                                                                                   \
  INSTANTIATE_TEST_SUITE_P(TestMeshesCreator, MeshesCreatorSpaceTimeSquareTestClass,               \
                           Values(0, 1, 2, 3));                                                    \
                                                                                                   \
  TEST_P(MeshesCreatorSpaceTimeSquareTestClass, TestDim) { EXPECT_EQ(mesh.dim(), 2); }             \
                                                                                                   \
  TEST_P(MeshesCreatorSpaceTimeSquareTestClass, TestCellCountGeometry) {                           \
    EXPECT_EQ(mesh.CellCountGeometry(), pow(8, level));                                            \
  }                                                                                                \
                                                                                                   \
  TEST_P(MeshesCreatorSpaceTimeSquareTestClass,                                                    \
         TestVertexCountGeometry){VerboseWarning("Not implemented", 1)}                            \
                                                                                                   \
  TEST_P(MeshesCreatorSpaceTimeSquareTestClass,                                                    \
         TestEdgeCountGeometry){VerboseWarning("Not implemented", 1)}                              \
                                                                                                   \
  TEST_P(MeshesCreatorSpaceTimeSquareTestClass,                                                    \
         TestFaceCountGeometry){VerboseWarning("Not implemented", 1)}                              \
                                                                                                   \
  TEST_P(MeshesCreatorSpaceTimeSquareTestClass, TestMeshWidthMin) {                                \
    EXPECT_EQ(minMeshWidth, pow(2, -level));                                                       \
  }                                                                                                \
                                                                                                   \
  TEST_P(MeshesCreatorSpaceTimeSquareTestClass, TestMeshWidthMax) {                                \
    EXPECT_EQ(maxMeshWidth, pow(2, -level));                                                       \
  }


#endif // TESTMESHESCREATOR_HPP
