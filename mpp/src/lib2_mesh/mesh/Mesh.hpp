#ifndef _MESH_H_
#define _MESH_H_

#include <functional>
#include <list>
#include <memory>
#include <set>
#include <string>
#include <utility>
#include <vector>

#include <gtest/gtest_prod.h>

#include "Cell.hpp"
#include "Celltype.hpp"
#include "CoarseGeometry.hpp"
#include "Edge.hpp"
#include "Face.hpp"
#include "ICell.hpp"
#include "Identify.hpp"
#include "LevelPair.hpp"
#include "MeshPart.hpp"
#include "MeshSettings.hpp"
#include "Point.hpp"
#include "ProcSet.hpp"
#include "Rule.hpp"
#include "Vertex.hpp"

class Mesh {
public:
  /**
   * @brief Construct 2^commSplit new independent Mesh objects with the lowest possible level.
   *
   * @param settings
   * @param commSplit
   */
  Mesh(const MeshSettings &settings, int commSplit);

  /**
   * @brief Construct a new Mesh object (copied, refined or/and mapped) from a base object.
   *
   * @param base
   * @param refine
   * @param refineInTime whether to refine in time or space
   * @param mapGeometry
   */
  Mesh(const Mesh &base, bool refine = false, bool refineInTime = false, bool mapGeometry = false);

  ~Mesh();

  bool IsSTMesh() const;

  std::string Name() const;

  void PrintBndIdsPerProc() const;

  const LevelPair &Level() const;

  const ProcSets &GetProcSets() const;

  std::pair<double, double> MeshWidth() const;

  const std::vector<int> &IncludedSubdomains() const;

  cell find_cell(const Point &center) const;

  vertex find_vertex(const Point &point) const;

  const std::string &DistributionName() const;

  void PrintInfo() const;

  double slice(int k) const;

  int steps() const;

  int slices() const;

  const TimeSteps &GetTimesteps() const;

  vector<double> GetTimeMidpointValues() const;

  double GetEndTime() const;

  double MaxMeshWidth() const;

  double MinMeshWidth() const;

  int dim() const { return settings.coarseGeometry->Dim(); }

  int BoundaryFacePart(const Point &bndFaceCenter) const;

  int CommSplit() const;

  cell cells() const;

  cell cells_end() const;

  void InitTValues() const;

  int GetBNDFaceCount(int i) const;

  int VertexCount() const;

  int VertexCountGeometry() const;

  int SpaceVertexCountGeometry() const;

  int EdgeCount() const;

  int EdgeCountGeometry() const;

  int FaceCount() const;

  int FaceCountGeometry() const;

  int CellCount() const;

  int CellCountGeometry() const;

  int SpaceCellCountGeometry(int slice_index) const;

  int ProcSetsCount() const;

  int ProcSetsCountWithoutInfty() const;

  int ProcSetsCountGeometry() const;

  int ProcSetsCountGeometryWithoutInfty() const;

  int BoundaryFaceCount() const;
protected:
  friend class IMatrixGraph;
  friend class EGMatrixGraph;
  friend class DGMatrixGraph;
  friend class TestRefineMesh;
  friend class Distribution;
  friend class TestMesh;
  friend class TestProcSets;
  friend class MeshInfoOnProc;
  friend class TestMeshesCreator;
  friend class CellBoundaryFaces;
  friend void match(const Mesh &mesh, const Mesh &other);
  friend void match(Mesh &mesh, CoarseGeometry &geo);
  FRIEND_TEST(TestMesh, GenerateFromMesh);

  template<typename S>
  friend LogTextStream<S> &operator<<(LogTextStream<S> &s, const Mesh &M);

  friend Buffer &operator>>(Buffer &b, Mesh &M);

  int verbose = 1;

  const MeshSettings &settings;

  LevelPair level;

  bool Parallel;

  bool identifyBnd = false;

  double maxMeshWidth = 0.0;

  double minMeshWidth = 0.0;

  Cells meshCells{};

  Edges meshEdges{};

  Faces meshFaces{};

  TimeSteps timesteps;

  Vertices meshVertices{};

  BoundaryFaces meshBndFaces{}; // Todo are they really needed?

  IdentifySets identifySets{};

  ProcSets procSets{};

  mutable std::vector<double> timeMidPointValues;

  std::function<void(Mesh &, const Mesh &base)> refineInSpaceWrapper;

  std::function<void(Mesh &, const Mesh &base)> refineInTimeWrapper = std::mem_fn(&Mesh::fillST);

  void insertVertex(const Point &p);

  bool removeVertex(const Point &p);

  void insertEdge(const Point &left, const Point &right);

  bool removeEdge(const Point &edgeCenter);

  bool refineEdge(const procset &p, const Mesh &base);

  void removeFace(const Point &F, const face &f, const Point &C);

  void removeFace(const cell &c, const Point &F);

  void removeBoundaryFace(const Point &x);

  void refineFace(const procset &p, const Mesh &base);

  void findMinMaxMeshWidth();

  bool isIncluded(int subdomain) const;

  cell insertCell(CELLTYPE tp, short subdomain, const std::vector<Point> &x, double a = 0.0,
                  double b = 0.0);

  void insertFace(const Point &faceCenter, const Point &cellCenter);

  void insertFace(const Point &faceCenter, const Face &faceValue);

  void finishParallel();

  void removeCell(cell c);

  void finish();

  void fill();

  void fillST(const Mesh &base);

  using CellRefinementFunction =
      std::function<const std::vector<Rule> &(const Cell &, std::vector<Point> &)>;
  void refineInSpaceInternal(const Mesh &base, const Cells &baseCells,
                             const CellRefinementFunction &cellRefinement);

  void refineInSpace(const Mesh &base);

  void refineInSpaceBaryCentric(const Mesh &base);

  using STCellsRefinementFunction = std::function<void(const Mesh &)>;

  void refineSTInternal(const Mesh &base, const STCellsRefinementFunction &cellsRefinement);

  void refineSTMeshInSpace(const Mesh &base);

  void refineSTMeshInTime(const Mesh &base);

  void mapGeometry(const Mesh &base);

  auto &bnd_ref() { return meshBndFaces.ref(); }

  void insertBoundaryFace(const Point &z, int part);

  void insertCellsWithTime(cell &);

  edge edges() const;

  edge edges_end() const;

  vertex vertices() const;

  vertex vertices_end() const;

  edge find_edge(const Point &center) const;

  face faces() const;

  face faces_end() const;

  face find_face(const Point &center) const;

  bnd_face bnd_faces() const;

  bnd_face bnd_faces_end() const;

  bnd_face find_bnd_face(const Point &z) const;

  cell find_cell_on_bndface(face f) const;

  procset procsets() const;

  procset procsets_end() const;

  procset find_procset(const Point &z) const;

  identifyset identifysets() const;

  identifyset identifysets_end() const;

  identifyset find_identifyset(const Point &z) const;

  bool master(const Point &z) const;

  int find_neighbour_face_id(const Point &f_c, const cell &cf) const;

  bool onBoundary(const Cell &c, int f) const;

  bool onBoundary(const Cell &c) const;

  bool identify() const;

  bool parallel() const;

  bool onBnd(const Point &x) const;

  std::string meshWidthStr() const;

  double t(int k) const;
};

template<typename S>
LogTextStream<S> &operator<<(LogTextStream<S> &s, const Mesh &M) {
  std::pair<const Vertices &, const ProcSets &> vertices_pair =
      std::make_pair(std::ref(M.meshVertices), std::ref(M.procSets));
  std::pair<const Edges &, const ProcSets &> edges_pair =
      std::make_pair(std::ref(M.meshEdges), std::ref(M.procSets));
  std::pair<const Faces &, const ProcSets &> faces_pair =
      std::make_pair(std::ref(M.meshFaces), std::ref(M.procSets));
  std::pair<const Cells &, const ProcSets &> cells_pair =
      std::make_pair(std::ref(M.meshCells), std::ref(M.procSets));
  std::pair<const BoundaryFaces &, const ProcSets &> bndfaces_pair =
      std::make_pair(std::ref(M.meshBndFaces), std::ref(M.procSets));

  s << "Vertices: " << M.VertexCount() << "\n"
    << vertices_pair << "Edges: " << M.EdgeCount() << "\n"
    << edges_pair << "Faces: " << M.FaceCount() << "\n"
    << faces_pair << "Cells: " << M.CellCount() << "\n"
    << cells_pair << "BoundaryFaces: " << M.BoundaryFaceCount() << "\n"
    << bndfaces_pair << "ProcSets: " << M.procSets.size() << "\n";
  s << M.procSets << "IdentifySets: " << M.identifySets.size() << "\n";
  s << M.identifySets;
  return s;
}

Buffer &operator>>(Buffer &b, Mesh &M);

template<class T>
std::ostream &operator<<(std::ostream &s,
                         std::pair<const MeshPart<Point, T> &, const ProcSets &> const &meshData) {
  const MeshPart<Point, T> &M = meshData.first;
  const ProcSets &PS = meshData.second;
  for (auto &[point, data] : M) {
    if constexpr (std::is_pointer_v<T>) {
      s << point << ": " << *data << " p ";
    } else {
      s << point << ": " << data << " p ";
    }

    procset ps = PS.find_procset(point);
    if (ps == PS.procsets_end()) {
      s << "[]";
    } else {
      s << (*ps).ToString();
    }
    s << "\n";
  }
  s << endl;
  return s;
}

template std::ostream &operator<<(std::ostream &s,
                                  std::pair<const Vertices &, const ProcSets &> const &meshData);

template std::ostream &operator<<(std::ostream &s,
                                  std::pair<const Cells &, const ProcSets &> const &meshData);

template std::ostream &operator<<(std::ostream &s,
                                  std::pair<const Faces &, const ProcSets &> const &meshData);

template std::ostream &
operator<<(std::ostream &s, std::pair<const BoundaryFaces &, const ProcSets &> const &meshData);

template std::ostream &operator<<(std::ostream &s,
                                  std::pair<const Edges &, const ProcSets &> const &meshData);

class BFParts {
  int bnd[6];
  int n;
  bool onbnd;
public:
  BFParts(const Mesh &M, const Cell &c);

  BFParts(const Mesh &M, const cell &c) : BFParts(M, *c){};

  int operator[](int i) const { return bnd[i]; }

  const int *operator()() const {
    if (onbnd) return bnd;
    return 0;
  }

  bool onBnd() const { return onbnd; }

  int size() const { return n; }

  template<typename S>
  friend LogTextStream<S> &operator<<(LogTextStream<S> &s, const BFParts &BF);
};

#endif // of #ifndef _MESH_H_
