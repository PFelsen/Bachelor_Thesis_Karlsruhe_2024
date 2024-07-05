
#include "TestMeshCollections.hpp"

TEST_P(TestMesh, LoopThroughCells) { LoopThroughCells(); }

TEST_P(TestMesh, LoopTroughVertices) { LoopTroughVertices(); }

TEST_P(TestMesh, LoopThroughEdges) { LoopThroughEdges(); }

TEST_P(TestMesh, LoopThroughFaces) { LoopThroughFaces(); }

TEST_P(TestMesh, LoopThroughBNDFaces) { LoopThroughBNDFaces(); }

#if SpaceDimension >= 2

INSTANTIATE_TEST_SUITE_P(MeshCollectionTest, TestMesh, testing::Values(testMeshParams2D));
#endif

INSTANTIATE_TEST_SUITE_P(UnitMeshCollectionTest, TestMesh, testing::ValuesIn(unitGeometries));

int main(int argc, char **argv) {
  return MppTest(MppTestBuilder(argc, argv).WithPPM().WithScreenLogging()).RUN_ALL_MPP_TESTS();
}