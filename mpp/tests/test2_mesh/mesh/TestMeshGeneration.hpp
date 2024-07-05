#ifndef TESTMESHGENERATION_H
#define TESTMESHGENERATION_H

#include "TestMesh.hpp"

extern void match(Mesh &mesh, CoarseGeometry &geo) {
  // Check vertex count and positions
  ASSERT_EQ(mesh.VertexCount(), geo.GetCoordinates().size());
  for (auto &coord : geo.GetCoordinates()) {
    ASSERT_TRUE(mesh.find_vertex(coord) != mesh.vertices_end());
  }

  // Check cell count and corners
  ASSERT_EQ(mesh.CellCount(), geo.GetCellIds().size());
  for (auto &cellID : geo.GetCellIds()) {
    std::vector<Point> corners = geo.Corners(cellID);

    Point midPoint = CalculateMidPoint(corners);

    auto meshCell = mesh.find_cell(midPoint);
    ASSERT_TRUE(meshCell != mesh.cells_end());

    //        ASSERT_EQ(meshCell.Type(), GetCellType(cellID.type(), meshCell.AsVector()));
    ASSERT_EQ(meshCell.Subdomain(), cellID.Subdomain());

    for (int i = 0; i < cellID.CornerIndices().size(); ++i) {
      ASSERT_EQ(meshCell.Corner(i), corners[i]);
    }
  }

  // TODO: Check face count and corners
  ASSERT_EQ(mesh.FaceCount(), geo.GetFaceIds().size());
  for (auto &faceID : geo.GetFaceIds()) {
    Point midPoint = geo.FaceCenter(faceID);
    auto meshFace = mesh.find_face(midPoint);
    ASSERT_TRUE(meshFace != mesh.faces_end());
  }
}

#define ASSERT_MATCHES(mesh, geo) match(mesh, *geo)


static vector<MeshTestParameters>
    unitGeometries{// Unit Interval
                   MeshTestParameters{std::vector<vector<double>>{{0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}},
                                      std::vector<vector<int>>{{2, 10, 0, 1}},
                                      std::vector<vector<int>>{{1, 5, 0}, {1, 15, 1}}}
#if SpaceDimension >= 2
                   ,
                   // Unit Triangle
                   MeshTestParameters{std::vector<vector<double>>{{0.0, 0.0, 0.0},
                                                                  {1.0, 0.0, 0.0},
                                                                  {0.0, 1.0, 0.0}},
                                      std::vector<vector<int>>{{3, 10, 0, 1, 2}},
                                      std::vector<vector<int>>{{2, 5, 0, 1},
                                                               {2, 15, 1, 2},
                                                               {2, 25, 2, 0}}},
                   // Unit Square
                   MeshTestParameters{std::vector<vector<double>>{{0.0, 0.0, 0.0},
                                                                  {1.0, 0.0, 0.0},
                                                                  {1.0, 1.0, 0.0},
                                                                  {0.0, 1.0, 0.0}},
                                      std::vector<vector<int>>{{4, 10, 0, 1, 2, 3}},
                                      std::vector<vector<int>>{{2, 5, 0, 1},
                                                               {2, 15, 1, 2},
                                                               {2, 25, 2, 3},
                                                               {2, 35, 3, 0}}}
#endif
#if SpaceDimension >= 3
                   ,
                   // Unit Tetrahedron
                   MeshTestParameters{std::vector<vector<double>>{{0.0, 0.0, 0.0},
                                                                  {1.0, 0.0, 0.0},
                                                                  {0.0, 1.0, 0.0},
                                                                  {0.0, 0.0, 1.0}},
                                      std::vector<vector<int>>{{4, 10, 0, 1, 2, 3}},
                                      std::vector<vector<int>>{{3, 5, 0, 1, 2},
                                                               {3, 15, 0, 1, 3},
                                                               {3, 25, 1, 2, 3},
                                                               {3, 35, 2, 0, 3}}},
                   // Unit Hexahedron
                   MeshTestParameters{std::vector<vector<double>>{{0.0, 0.0, 0.0},
                                                                  {1.0, 0.0, 0.0},
                                                                  {1.0, 1.0, 0.0},
                                                                  {0.0, 1.0, 0.0},
                                                                  {0.0, 0.0, 1.0},
                                                                  {1.0, 0.0, 1.0},
                                                                  {1.0, 1.0, 1.0},
                                                                  {0.0, 1.0, 1.0}},
                                      std::vector<vector<int>>{{8, 10, 0, 1, 2, 3, 4, 5, 6, 7}},
                                      std::vector<vector<int>>{{4, 5, 0, 1, 2, 3},
                                                               {4, 15, 0, 1, 5, 4},
                                                               {4, 25, 1, 2, 6, 5},
                                                               {4, 35, 2, 3, 7, 6},
                                                               {4, 45, 3, 0, 4, 7},
                                                               {4, 55, 4, 5, 6, 7}}}
#endif
    };

#endif // TESTMESHGENERATION_H
