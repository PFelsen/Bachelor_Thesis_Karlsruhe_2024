#include "LagrangeDoF.hpp"
#include "MeshesCreator.hpp"
#include "RTDoF.hpp"
#include "TestEnvironment.hpp"
#include "VectorMatrixBase.hpp"

#include <algorithm>

struct TestParam {
  std::string distributeName;
  std::string dofName;
  int degree;
  int numberOfDofsAtNodalPoint;

  friend void PrintTo(const TestParam &param, std::ostream *os) {
    *os << "DistName: " << param.distributeName << " dofName: " << param.dofName
        << " degree: " << param.degree << " numberOfDofsAtNodalPoint"
        << param.numberOfDofsAtNodalPoint;
  }
};

class MatrixGraphTest : public TestWithParam<TestParam> {
protected:
  int pLevel = 2;
  int level = 2;
  Meshes *meshes = nullptr;
  const Mesh &mesh;
  MatrixGraph *graph = nullptr;

  int numberOfRows = -1;
  int numberOfDofsAtNodalPoint = -1;
  int vectorSize = -1;
  int matrixSize = -1;

  MatrixGraphTest(std::string meshName) :
      meshes(MeshesCreator(meshName)
                 .WithDistribute(GetParam().distributeName)
                 .WithPLevel(pLevel)
                 .WithLevel(level)
                 .Create()),
      mesh(meshes->fine()) {
    computeNumberOfRows(meshName);
    numberOfDofsAtNodalPoint = GetParam().numberOfDofsAtNodalPoint;
    vectorSize = numberOfRows * numberOfDofsAtNodalPoint;
    std::unique_ptr<IDoF> _dof = nullptr;
    if (GetParam().dofName == "LagrangeDoF") {
      _dof = std::make_unique<LagrangeDoF>(GetParam().degree, GetParam().numberOfDofsAtNodalPoint);
    } else if (GetParam().dofName == "RTDoF") {
      _dof = std::make_unique<RTDoF>(GetParam().degree, GetParam().numberOfDofsAtNodalPoint);
    } else {
      Exit("DoF not implemented in MatrixGraphTest")
    }
    graph = new MatrixGraph(mesh, std::move(_dof));
    computeMatrixSize();
  }

  ~MatrixGraphTest() {
    delete graph;
    delete meshes;
  }
private:
  bool entryExists(const std::vector<Point> &v, const Point &p) const {
    return std::find(v.begin(), v.end(), p) != v.end();
  }

  void computeNumberOfRows(const std::string &meshName) {
    if (GetParam().dofName == "LagrangeDoF") {
      switch (GetParam().degree) {
      case 0:
        numberOfRows = mesh.CellCount();
        break;
      case 1:
        numberOfRows = mesh.VertexCount();
        break;
      case 2:
        numberOfRows = 0;
        switch (mesh.dim()) {
        case 1:
          numberOfRows = mesh.EdgeCount() + mesh.CellCount();
          break;
        case 3:
        case 2:
          if (meshName == "Hexahedron") numberOfRows += mesh.CellCount() + mesh.FaceCount();
          if (meshName == "Square") numberOfRows += mesh.CellCount();
          numberOfRows += mesh.VertexCount() + mesh.EdgeCount();
          break;
        default:
          Exit("Celltype not implemented")
        }
        break;
      default:
        Exit("Degree not implemented")
      }
    } else if (GetParam().dofName == "RTDoF") {
      if (GetParam().degree != 0) Exit("Degree not implemented") numberOfRows = mesh.FaceCount();
    } else {
      Exit("DoF not implemented in MatrixGraphTest")
    }
  }

  void computeMatrixSize() {
    matrixSize = 0;
    addCenters();
    addVertices();
    if (mesh.dim() > 1) addFaces();
    if (mesh.dim() > 2) addEdges();
    matrixSize *= numberOfDofsAtNodalPoint * numberOfDofsAtNodalPoint;
  }

  void addCenters() {
    for (cell cll = mesh.cells(); cll != mesh.cells_end(); ++cll) {
      std::vector<Point> nodalPointsOnCell = graph->GetDoF().GetNodalPoints(*cll);
      if (entryExists(nodalPointsOnCell, cll->first)) matrixSize += nodalPointsOnCell.size();
    }
  }

  bool addNodalPoints(std::vector<Point> &nodalPoints, const cell &cll,
                      const Point &consideredPoint) {
    std::vector<Point> nP = graph->GetDoF().GetNodalPoints(*cll);
    if (!entryExists(nP, consideredPoint)) // consideredPoint itself is no nodal point
      return false;
    for (int k = 0; k < nP.size(); ++k) {
      if (!entryExists(nodalPoints, nP[k])) nodalPoints.push_back(nP[k]);
    }
    return true;
  }

  void addVertices() {
    for (auto vtx = graph->vertices(); vtx != graph->vertices_end(); ++vtx) {
      // find all nodal points of cells containing vtx
      std::vector<Point> nodalPoints{};
      for (cell cll = mesh.cells(); cll != mesh.cells_end(); ++cll) {
        for (int i = 0; i < cll->second->Corners(); ++i) {
          if (vtx->first == cll->second->Corner(i)) {
            if (!addNodalPoints(nodalPoints, cll, vtx->first))
              return; // vtx itself is no nodal point
            break;
          }
        }
      }
      matrixSize += nodalPoints.size();
    }
  }

  void addFaces() {
    for (auto fce = graph->faces(); fce != graph->faces_end(); ++fce) {
      // find all nodal points of cells containing fce
      std::vector<Point> nodalPoints{};
      if (fce->second.Left() != Infty) {
        cell c = graph->find_cell(fce->second.Left());
        if (c != mesh.cells_end())
          if (!addNodalPoints(nodalPoints, c, fce->first)) return;
      }
      if (fce->second.Right() != Infty) {
        cell c = graph->find_cell(fce->second.Right());
        if (c != mesh.cells_end())
          if (!addNodalPoints(nodalPoints, c, fce->first)) return;
      }
      matrixSize += nodalPoints.size();
    }
  }

  void addEdges() {
    for (auto edg = graph->edges(); edg != graph->edges_end(); ++edg) {
      // find all nodal points of cells containing edg
      std::vector<Point> nodalPoints{};
      for (cell cll = mesh.cells(); cll != mesh.cells_end(); ++cll) {
        for (int i = 0; i < cll->second->Edges(); ++i) {
          if (edg->first == cll->second->Edge(i)) {
            if (!addNodalPoints(nodalPoints, cll, edg->first))
              return; // edg itself is no nodal point
            break;
          }
        }
      }
      matrixSize += nodalPoints.size();
    }
  }
};

#define MATRIXGRAPH_TESTS_LAGRANGE(MeshName, DistributeName)                                       \
  class MatrixGraph##DistributeName##MeshName##TestLagrange : public MatrixGraphTest {             \
  public:                                                                                          \
    MatrixGraph##DistributeName##MeshName##TestLagrange() : MatrixGraphTest(#MeshName) {}          \
  };                                                                                               \
                                                                                                   \
  TEST_P(MatrixGraph##DistributeName##MeshName##TestLagrange, NumberOfRowsTest) {                  \
    EXPECT_EQ(graph->rowsize(), numberOfRows);                                                     \
  }                                                                                                \
  TEST_P(MatrixGraph##DistributeName##MeshName##TestLagrange, NumberOfNodalPointsTest) {           \
    for (int i = 0; i < numberOfRows; ++i)                                                         \
      EXPECT_EQ(graph->Dof(i), numberOfDofsAtNodalPoint);                                          \
  }                                                                                                \
                                                                                                   \
  TEST_P(MatrixGraph##DistributeName##MeshName##TestLagrange, VectorSizeTest) {                    \
    EXPECT_EQ(graph->size(), vectorSize);                                                          \
  }                                                                                                \
                                                                                                   \
  TEST_P(MatrixGraph##DistributeName##MeshName##TestLagrange, VectorIndexTest) {                   \
    for (int i = 0; i < numberOfRows; ++i)                                                         \
      EXPECT_EQ(graph->Index(i), i *numberOfDofsAtNodalPoint);                                     \
  }                                                                                                \
                                                                                                   \
  TEST_P(MatrixGraph##DistributeName##MeshName##TestLagrange, MatrixSizeTest) {                    \
    EXPECT_EQ(graph->Size(), matrixSize);                                                          \
  }                                                                                                \
                                                                                                   \
  INSTANTIATE_TEST_SUITE_P(MatrixGraphTest, MatrixGraph##DistributeName##MeshName##TestLagrange,   \
                           Values(TestParam{#DistributeName, "LagrangeDoF", 0, 1},                 \
                                  TestParam{#DistributeName, "LagrangeDoF", 0, 2},                 \
                                  TestParam{#DistributeName, "LagrangeDoF", 1, 1},                 \
                                  TestParam{#DistributeName, "LagrangeDoF", 1, 2},                 \
                                  TestParam{#DistributeName, "LagrangeDoF", 2, 1},                 \
                                  TestParam{#DistributeName, "LagrangeDoF", 2, 2}));

#define MATRIXGRAPH_TESTS_RT(MeshName, DistributeName)                                             \
  class MatrixGraph##DistributeName##MeshName##TestRT : public MatrixGraphTest {                   \
  public:                                                                                          \
    MatrixGraph##DistributeName##MeshName##TestRT() : MatrixGraphTest(#MeshName) {}                \
  };                                                                                               \
                                                                                                   \
  TEST_P(MatrixGraph##DistributeName##MeshName##TestRT, NumberOfRowsTest) {                        \
    EXPECT_EQ(graph->rowsize(), numberOfRows);                                                     \
  }                                                                                                \
  TEST_P(MatrixGraph##DistributeName##MeshName##TestRT, NumberOfNodalPointsTest) {                 \
    for (int i = 0; i < numberOfRows; ++i)                                                         \
      EXPECT_EQ(graph->Dof(i), numberOfDofsAtNodalPoint);                                          \
  }                                                                                                \
                                                                                                   \
  TEST_P(MatrixGraph##DistributeName##MeshName##TestRT, VectorSizeTest) {                          \
    EXPECT_EQ(graph->size(), vectorSize);                                                          \
  }                                                                                                \
                                                                                                   \
  TEST_P(MatrixGraph##DistributeName##MeshName##TestRT, VectorIndexTest) {                         \
    for (int i = 0; i < numberOfRows; ++i)                                                         \
      EXPECT_EQ(graph->Index(i), i *numberOfDofsAtNodalPoint);                                     \
  }                                                                                                \
                                                                                                   \
  TEST_P(MatrixGraph##DistributeName##MeshName##TestRT, MatrixSizeTest) {                          \
    EXPECT_EQ(graph->Size(), matrixSize);                                                          \
  }                                                                                                \
                                                                                                   \
  INSTANTIATE_TEST_SUITE_P(MatrixGraphTest, MatrixGraph##DistributeName##MeshName##TestRT,         \
                           Values(TestParam{#DistributeName, "RTDoF", 0, 1},                       \
                                  TestParam{#DistributeName, "RTDoF", 0, 2}));

MATRIXGRAPH_TESTS_LAGRANGE(Interval, Stripes)
MATRIXGRAPH_TESTS_LAGRANGE(Interval, RCB)

MATRIXGRAPH_TESTS_LAGRANGE(Triangle, Stripes)
MATRIXGRAPH_TESTS_LAGRANGE(Triangle, RCB)

MATRIXGRAPH_TESTS_RT(Triangle, Stripes)
MATRIXGRAPH_TESTS_RT(Triangle, RCB)

MATRIXGRAPH_TESTS_LAGRANGE(Square, Stripes)
MATRIXGRAPH_TESTS_LAGRANGE(Square, RCB)

MATRIXGRAPH_TESTS_LAGRANGE(Tetrahedron, Stripes)
MATRIXGRAPH_TESTS_LAGRANGE(Tetrahedron, RCB)

// TODO: Fix Tetrahedron RT-Tests.

// MATRIXGRAPH_TESTS_RT(Tetrahedron, Stripes)
// MATRIXGRAPH_TESTS_RT(Tetrahedron, RCB)

MATRIXGRAPH_TESTS_LAGRANGE(Hexahedron, Stripes)
MATRIXGRAPH_TESTS_LAGRANGE(Hexahedron, RCB)

int main(int argc, char **argv) {
  MppTest mppTest = MppTestBuilder(argc, argv).WithPPM();
  return mppTest.RUN_ALL_MPP_TESTS();
}
