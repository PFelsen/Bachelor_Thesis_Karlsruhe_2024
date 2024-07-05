#include "TestMeshGeneration.hpp"

#include "TestMesh.hpp"

#include <Mesh.cpp>

TEST_P(TestMesh, GenerateFromGeometry) { ASSERT_MATCHES(mesh, geometry); }

TEST_P(TestMesh, GenerateFromMesh) {
  Mesh newMesh(mesh);
  ASSERT_MESH_EQ(mesh, newMesh);

  // Verify its a deep copy
  ASSERT_NE(newMesh.cells(), newMesh.cells_end());
  const auto &oldCell = mesh.cells();
  const auto &newCell = newMesh.find_cell(oldCell());
  ASSERT_NE(oldCell->second, newCell->second);

  ASSERT_NE(newMesh.edges(), newMesh.edges_end());
  const auto &oldEdge = mesh.edges();
  const auto &newEdge = newMesh.find_edge(oldEdge());
  ASSERT_NE(oldEdge->second, newEdge->second);

  ASSERT_NE(newMesh.vertices(), newMesh.vertices_end());
  const auto &oldVertex = mesh.vertices();
  const auto &newVertex = newMesh.find_vertex(oldVertex());
  ASSERT_NE(oldVertex->second, newVertex->second);
}

INSTANTIATE_TEST_SUITE_P(UnitGeometryGenerationTest, TestMesh, ValuesIn(unitGeometries));

int main(int argc, char **argv) {
  MppTest mppTest = MppTestBuilder(argc, argv).WithPPM();
  return mppTest.RUN_ALL_MPP_TESTS();
}