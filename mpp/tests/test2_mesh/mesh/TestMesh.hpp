#ifndef TESTMESH_H
#define TESTMESH_H

#include <algorithm>
#include <vector>

#include <Celltype.hpp>
#include <CoarseGeometry.hpp>
#include <Mesh.hpp>
#include <Point.hpp>
#include <gtest/gtest.h>

#include "Edge.hpp"
#include "Face.hpp"
#include "Identify.hpp"
#include "ProcSet.hpp"
#include "TestEnvironment.hpp"

// TODO: Change Input of tests to 'new' structures (Coordinates,CellIds,FaceIds) instead of using
// the old versions(vec<vec<double>>,vec<vec<int>>,vec<vec<int>>) and converting it.

CELLTYPE getCellType(int calculatedSpaceDim, int tp, int numberOfCorners) {
  if (calculatedSpaceDim == 1) return INTERVAL;
  else if (calculatedSpaceDim == 2) {
    if (tp == 3 && numberOfCorners == 3) return TRIANGLE;
    if (tp == 4 && numberOfCorners == 4) return QUADRILATERAL;
    if (tp == 4 && numberOfCorners == 8) return QUADRILATERAL2;
  } else {
    if (tp == 4 && numberOfCorners == 4) return TETRAHEDRON;
    if (tp == 8 && numberOfCorners == 8) return HEXAHEDRON;
    if (tp == 8 && numberOfCorners == 20) return HEXAHEDRON20;
    if (tp == 8 && numberOfCorners == 27) return HEXAHEDRON27;
  }
  THROW("No corresponding CellType")
}

struct MeshTestParameters {
  std::vector<Point> coordinates;
  vector<FaceId> faceIds;
  vector<CellId> cellIds;
  int calculatedSpaceDim;
  CELLTYPE type;

  std::vector<Point> GetCoords() const { return coordinates; }

  vector<FaceId> GetFaceIds() const { return faceIds; }

  vector<CellId> GetCellIds() const { return cellIds; }

  void InitCoordinates(const vector<vector<double>> &corners) {
    for (const auto &point : corners) {
      double z[]{0.0, 0.0, 0.0, 0.0};
      for (int j = 0; j < min(3, int(point.size())); ++j)
        z[j] = point[j];
      coordinates.emplace_back(Point(z));
    }
    calculatedSpaceDim = 1;
    if constexpr (SpaceDimension > 1) {
      for (int i = 0; i < coordinates.size(); ++i)
        if (coordinates[i][1] != 0) calculatedSpaceDim = 2;
    }
    if constexpr (SpaceDimension > 2) {
      for (int i = 0; i < coordinates.size(); ++i)
        if (coordinates[i][2] != 0) calculatedSpaceDim = 3;
    }
  }

  void InitFaceIds(const vector<vector<int>> &values) {
    for (const auto &value : values) {
      if (value.size() < 3) continue;
      int n = value[0];
      int f = value[1];
      std::vector<int> id_vec(value.size() - 2);

      Assert(n == id_vec.size());

      for (int j = 2; j < value.size(); ++j)
        id_vec[j - 2] = value[j];
      faceIds.emplace_back(f, std::move(id_vec));
    }
  }

  void InitCellIds(const vector<vector<int>> &values) {
    for (const auto &value : values) {
      if (value.size() < 3) continue;
      int n = value[0];
      int f = value[1];
      std::vector<int> id_vec(value.size() - 2);

      for (int j = 2; j < value.size(); ++j)
        id_vec[j - 2] = value[j];

      CELLTYPE type = getCellType(calculatedSpaceDim, n, (int)id_vec.size());
      cellIds.emplace_back(type, f, std::move(id_vec));
    }
  }

  MeshTestParameters(std::vector<std::vector<double>> corners, std::vector<std::vector<int>> cells,
                     std::vector<std::vector<int>> faces) {
    InitCoordinates(corners);
    InitCellIds(cells);
    InitFaceIds(faces);
  }
};

class TestMesh : public TestWithParam<MeshTestParameters> {
  MeshSettings createSettingsByCoarseGeometry(std::shared_ptr<CoarseGeometry> geometry) {
    MeshSettings settings;
    settings.coarseGeometry = std::move(geometry);
    return settings;
  }
protected:
  std::shared_ptr<CoarseGeometry> geometry;
  MeshSettings settings;
  Mesh mesh;

  explicit TestMesh() :
      geometry(std::make_shared<CoarseGeometry>(GetParam().GetCoords(), GetParam().GetCellIds(),
                                                GetParam().GetFaceIds())),
      settings(createSettingsByCoarseGeometry(geometry)), mesh(settings, 0) {
    mesh.fill();
  }

  void LoopThroughBNDFaces() {
    int iteratorCount = 0;
    for (auto it = mesh.bnd_faces(); it != mesh.bnd_faces_end(); ++it) {
      iteratorCount++;
    }
    // Todo is this right?
    EXPECT_EQ(iteratorCount, geometry->GetFaceIds().size());
  }

  void LoopThroughFaces() {
    int iteratorCount = 0;
    for (auto it = mesh.faces(); it != mesh.faces_end(); ++it) {
      iteratorCount++;
    }
    EXPECT_GE(iteratorCount, geometry->GetFaceIds().size());
  }

  void LoopThroughEdges() {
    int iteratorCount = 0;
    for (auto it = mesh.edges(); it != mesh.edges_end(); ++it) {
      iteratorCount++;
    }
    // Todo is this right?
    EXPECT_GE(iteratorCount, geometry->GetFaceIds().size());
  }

  void LoopTroughVertices() {
    int iteratorCount = 0;
    for (auto it = mesh.vertices(); it != mesh.vertices_end(); ++it) {
      iteratorCount++;
    }
    EXPECT_EQ(iteratorCount, geometry->GetCoordinates().size());
  }

  void LoopThroughCells() {
    int iteratorCount = 0;
    for (auto it = mesh.cells(); it != mesh.cells_end(); ++it) {
      iteratorCount++;
    }
    EXPECT_EQ(iteratorCount, geometry->GetCellIds().size());
  }
};

static Point CalculateMidPoint(const std::vector<Point> &points) {
  Point midPoint;
  for (int i = 0; i < points.size(); ++i)
    midPoint += points[i];
  return midPoint / points.size();
}

void match(const Mesh &mesh, const Mesh &other) {
  ASSERT_EQ(mesh.Level(), other.Level());

  // Check vertex count and positions
  ASSERT_EQ(mesh.VertexCount(), other.VertexCount());
  for (vertex v = mesh.vertices(); v != mesh.vertices_end(); ++v) {
    ASSERT_TRUE(other.find_vertex(v()) != other.vertices_end());
  }

  // Check cell count and corners
  ASSERT_EQ(mesh.CellCount(), other.CellCount());
  for (cell c = mesh.cells(); c != mesh.cells_end(); ++c) {
    cell otherCell = other.find_cell(c());
    ASSERT_TRUE(otherCell != other.cells_end());

    ASSERT_EQ(c.Type(), otherCell.Type());
    ASSERT_EQ(c.Subdomain(), otherCell.Subdomain());

    for (int i = 0; i < c.Corners(); ++i) {
      ASSERT_EQ(c.Corner(i), otherCell.Corner(i));
    }
  }

  // TODO: Check face count and corners
  ASSERT_EQ(mesh.FaceCount(), other.FaceCount());
  for (face f = mesh.faces(); f != mesh.faces_end(); ++f) {
    auto otherFace = mesh.find_face(f());
    ASSERT_TRUE(otherFace != other.faces_end());
    ASSERT_EQ(f.Left(), otherFace.Left());
    ASSERT_EQ(f.Right(), otherFace.Right());
  }

  ASSERT_EQ(mesh.EdgeCount(), other.EdgeCount());
  for (edge e = mesh.edges(); e != mesh.edges_end(); ++e) {
    auto otherEdge = mesh.find_edge(e());
    ASSERT_NE(otherEdge, other.edges_end());
    ASSERT_EQ(e.Left(), otherEdge.Left());
    ASSERT_EQ(e.Right(), otherEdge.Right());
  }

  ASSERT_EQ(mesh.ProcSetsCount(), other.ProcSetsCount());
  for (procset ps = mesh.procsets(); ps != mesh.procsets_end(); ++ps) {
    auto otherPs = mesh.find_procset(ps());
    ASSERT_NE(otherPs, other.procsets_end());
    ASSERT_EQ(ps.size(), otherPs.size());
  }

  ASSERT_EQ(mesh.BoundaryFaceCount(), other.BoundaryFaceCount());
  for (bnd_face bf = mesh.bnd_faces(); bf != mesh.bnd_faces_end(); ++bf) {
    auto otherBf = other.find_bnd_face(bf());
    ASSERT_NE(otherBf, other.bnd_faces_end());
    ASSERT_EQ(bf.Part(), otherBf.Part());
  }

  ASSERT_EQ(mesh.identifySets.size(), other.identifySets.size());
  for (identifyset is = mesh.identifysets(); is != mesh.identifysets_end(); ++is) {
    auto otherIs = other.find_identifyset(is());
    ASSERT_NE(otherIs, other.identifysets_end());
    ASSERT_EQ(is.size(), otherIs.size());
  }

  ASSERT_EQ(mesh.steps(), other.steps());

  ASSERT_EQ(mesh.GetTimeMidpointValues().size(), other.GetTimeMidpointValues().size());
}

#define ASSERT_MESH_EQ(mesh, other) match(mesh, other)

#endif // TESTMESH_H
