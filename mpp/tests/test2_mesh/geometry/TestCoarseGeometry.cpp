
#include "TestCoarseGeometry.hpp"

#include <Celltype.hpp>
#include <CoarseGeometry.hpp>
#include <MeshesCreator.hpp>
#include <gtest/gtest.h>

#if SpaceDimension >= 3

class CoarseGeometryTest : public testing::Test {
protected:
  CoarseGeometry geoTest;

  CoarseGeometryTest() :
      geoTest(Coordinates{{0.0, 0.0, 0.0},
                          {1.0, 0.0, 0.0},
                          {1.0, 1.0, 0.0},
                          {0.0, 1.0, 0.0},
                          {0.0, 0.0, 1.0},
                          {1.0, 0.0, 1.0},
                          {1.0, 1.0, 1.0},
                          {0.0, 1.0, 1.0}},
              CellIds{{HEXAHEDRON, 10, {0, 1, 2, 3, 4, 5, 6, 7}}},
              FaceIds{{1, {0, 1, 2, 3}},
                      {1, {0, 1, 5, 4}},
                      {1, {1, 2, 6, 5}},
                      {1, {2, 3, 7, 6}},
                      {1, {3, 0, 4, 7}},
                      {1, {4, 5, 6, 7}}},
              VertexDataList{DataContainer({1.0, 2.0, 3.0}), DataContainer({1.0, 2.0, 3.0}),
                             DataContainer({1.0, 2.0, 3.0}), DataContainer({1.0, 2.0, 3.0}),
                             DataContainer({1.0, 2.0, 3.0}), DataContainer({1.0, 2.0, 3.0}),
                             DataContainer({1.0, 2.0, 3.0}), DataContainer({1.0, 2.0, 3.0})},
              CellDataList{DataContainer({0.0, 1.0, 2.0})}, "TestCube") {}
};

TEST_F(CoarseGeometryTest, TestLegacyReading) {
  CoarseGeometry geoLegacy("TestCubeLegacy");
  EXPECT_GEO_EQ(geoLegacy, geoTest);
}

TEST_F(CoarseGeometryTest, TestVtuReading) {
  CoarseGeometry geoLegacy("TestCubeVtu");
  EXPECT_GEO_EQ(geoLegacy, geoTest);
}

TEST_F(CoarseGeometryTest, TestLegacyWithData) {
  CoarseGeometry geoLegacy("TestCubeLegacyWithData");
  EXPECT_DATAGEO_EQ(geoLegacy, geoTest);
}

TEST_F(CoarseGeometryTest, TestVtuWithData) {
  CoarseGeometry geoLegacy("TestCubeVtuWithData");
  EXPECT_DATAGEO_EQ(geoLegacy, geoTest);
}

TEST_F(CoarseGeometryTest, TestVtuWithCelltype) {
  CoarseGeometry geoLegacy("TestCubeCelltype");
  EXPECT_GEO_EQ(geoLegacy, geoTest);
}

TEST_F(CoarseGeometryTest, TestVtuWithIntegerCelltype) {
  CoarseGeometry geoLegacy("TestCubeCelltypeAsInt");
  EXPECT_GEO_EQ(geoLegacy, geoTest);
}

TEST_F(CoarseGeometryTest, TestVtuWithSubdomain) {
  CoarseGeometry geoLegacy("TestCubeSubdomain");
  EXPECT_GEO_EQ(geoLegacy, geoTest);
}

TEST_F(CoarseGeometryTest, TestVtuWithBoundary) {
  CoarseGeometry geoLegacy("TestCubeBoundary");
  EXPECT_GEO_EQ(geoLegacy, geoTest);
}

TEST_F(CoarseGeometryTest, TestVtuWithAll) {
  CoarseGeometry geoLegacy("TestCubeAll");
  EXPECT_GEO_EQ(geoLegacy, geoTest);
}

TEST_F(CoarseGeometryTest, TestVtuWithDefaultData) {
  CoarseGeometry geoLegacy("TestCubeDefaultData");
  EXPECT_DATAGEO_EQ(geoLegacy, geoTest);
}

TEST_F(CoarseGeometryTest, TestWriting) {
  CoarseGeometry geoLegacy("TestCubeLegacy");
  CoarseGeometry geoVtu("TestCubeLegacy");
  EXPECT_GEO_EQ(geoLegacy, geoVtu);
}

TEST_F(CoarseGeometryTest, TestWritingWithData) {
  CoarseGeometry geoLegacy("TestCubeLegacyWithData");
  CoarseGeometry geoVtu("TestCubeLegacyWithData");
  EXPECT_GEO_EQ(geoLegacy, geoVtu);
}
#endif

int main(int argc, char **argv) {
  Config::SetGeoPath(string(ProjectBuildDir) + "/data/geo");
  MppTest mppTest = MppTestBuilder(argc, argv)
                        .WithPPM()
                        .WithConfigEntry("Mesh", "TestCubeLegacy")
                        .WithGeoPath(string(ProjectMppDir) + "/tests/test2_mesh/geometry/")
                        .WithScreenLogging();
  return mppTest.RUN_ALL_MPP_TESTS();
}