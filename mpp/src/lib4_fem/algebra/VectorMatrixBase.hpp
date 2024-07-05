#ifndef VECTORMATRIXBASE_HPP
#define VECTORMATRIXBASE_HPP

#include "IDiscretization.hpp"
#include "MixedDoF.hpp"

class VectorMatrixBase {
protected:
  std::shared_ptr<const IDiscretization> disc;

#ifdef BUILD_IA
  std::shared_ptr<const IAIDiscretization> iadisc;
#endif

  // reference to avoid multiple access to matrix graph map in discretization
  IMatrixGraph &graph;

  VectorMatrixBase(std::shared_ptr<const IDiscretization> disc, LevelPair levels) :
      disc(std::move(disc)), graph((*this->disc)(levels)) {}

  VectorMatrixBase(const VectorMatrixBase &g) : disc(g.disc), graph(g.graph) {
#ifdef BUILD_IA
    iadisc = g.iadisc;
#endif
  }
public:
  const IDiscretization &GetDisc() const { return *disc; }

  std::shared_ptr<const IDiscretization> GetSharedDisc() const { return disc; }

  const Meshes &GetMeshes() const { return GetDisc().GetMeshes(); }

  template<typename T, int sDim, int tDim>
  const IDiscretizationT<T, sDim, tDim> &GetDiscT() const;

#ifdef BUILD_IA

  const IAIDiscretization &GetIADisc() const { return *iadisc; }

  std::shared_ptr<const IAIDiscretization> GetSharedIADisc() const { return iadisc; }

#endif

  const ProcSets &GetProcSets() const { return graph.GetProcSets(); }

  const IMatrixGraph &GetMatrixGraph() const { return graph; };

  int CommSplit() const { return graph.CommSplit(); }

  LevelPair Level() const { return graph.GetLevel(); }

  int SpaceLevel() const { return Level().space; }

  int TimeLevel() const { return Level().time; }

  void PrintInfo() const { graph.PrintInfo(); }

  row rows() const { return graph.rows(); }

  row rows_end() const { return graph.rows_end(); }

  row find_row(const Point &z) const { return graph.find_row(z); }

  int pMatrixSize() const { return (int)graph.pMatrixSize(); }

  int Id(const Point &z) const { return graph.GetRows().RowId(z); }

  int Idx(const Point &z) const { return graph.GetRows().RowIdChecked(z); }

  int size() const { return graph.size(); }

  int Size() const { return (int)graph.Size(); }

  int pSize() const { return graph.pSize(); }

  int Index(int i) const { return graph.Index(i); }

  int Column(int i) const { return graph.Column(i); }

  int Diag(int i) const { return graph.Diag(i); }

  int Entry(int i) const { return graph.Entry(i); }

  bool SingleEntry(int i) const { return graph.SingleEntry(i); }

  int Dof(int i) const { return graph.Dof(i); }

  int GetEntry(const row &r0, const row &r1) const { return graph.GetRows().Index(r0, r1); }

  int GetEntryX(const row &r0, const row &r1) const { return graph.GetRows().IndexChecked(r0, r1); }

  int GetDoubleEntryX(const row &r0, const row &r1) const {
    return graph.GetRows().DoubleIndexChecked(r0, r1);
  }

  int nR() const { return graph.GetRows().size(); }

  template<typename S>
  friend LogTextStream<S> &operator<<(LogTextStream<S> &s, const VectorMatrixBase &G) {
    return s << G.disc;
  }

  int rowsize() const { return graph.rowsize(); }

  const Mesh &GetMesh() const { return graph.GetMesh(); }

  cell cells() const { return graph.cells(); }

  cell cells_end() const { return graph.cells_end(); }

  cell find_cell(const Point &z) const { return graph.find_cell(z); }

  cell overlap() const { return graph.overlap(); }

  cell overlap_end() const { return graph.overlap_end(); }

  face faces() const { return graph.faces(); }

  face find_face(const Point &center) const { return graph.find_face(center); }

  vertex vertices() const { return graph.vertices(); }

  vertex vertices_end() const { return graph.vertices_end(); }

  vertex find_vertex(const Point &p) const { return graph.find_vertex(p); }

  edge edges() const { return graph.edges(); }

  edge edges_end() const { return graph.edges_end(); }

  edge find_edge(const Point &center) const { return graph.find_edge(center); }

  cell find_neighbour_cell(const cell &c, int f) const { return graph.find_neighbour_cell(c, f); }

  cell find_neighbour_cell(const Cell &c, int f) const { return graph.find_neighbour_cell(c, f); }

  int find_neighbour_face_id(const Point &f_c, const cell &cf) const {
    return graph.find_neighbour_face_id(f_c, cf);
  }

  cell find_cell_or_overlap_cell(const Point &z) const {
    return graph.find_cell_or_overlap_cell(z);
  }

  Point neighbour_center(const Cell &c, int i) const { return graph.neighbour_center(c, i); }

  bool has_cell_or_overlap_cell(const Point &z) const { return graph.has_cell_or_overlap_cell(z); }

  bool OnBoundary(const Cell &c) const { return graph.OnBoundary(c); }

  bool has_previous_cell(const cell &c) const { return graph.has_previous_cell(c); }

  bool OnBoundary(const Cell &c, int f) const { return graph.OnBoundary(c, f); }

  cell find_previous_cell(const cell &c) const { return graph.find_previous_cell(c); }

  face faces_end() const { return graph.faces_end(); }

  cell find_overlap_cell(const Point &z) const { return graph.find_overlap_cell(z); }

  double t(int k) const { return graph.t(k); }

  bool onBnd(const Point &x) const { return graph.onBnd(x); }

  int BndPart(const Point &x) const { return graph.BndPart(x); }

  bnd_face bnd_faces() const { return graph.bnd_faces(); }

  bnd_face bnd_faces_end() const { return graph.bnd_faces_end(); }

  bnd_face find_bnd_face(const Point &z) const { return graph.find_bnd_face(z); }

  identifyset identifysets() const { return graph.GetIdentifySets().identifysets(); }

  identifyset identifysets_end() const { return graph.GetIdentifySets().identifysets_end(); }

  identifyset find_identifyset(const Point &z) const { return graph.find_identifyset(z); }

  bool identify() const { return graph.identify(); }

  bool parallel() const { return graph.parallel(); }

  int dim() const { return graph.dim(); }

  procset procsets() const { return graph.procsets(); }

  procset procsets_end() const { return graph.procsets_end(); }

  procset find_procset(const Point &z) const { return graph.find_procset(z); }

  bool master(const Point &z) const { return graph.GetProcSets().master(z); }

  const SubVectorMask &Mask(const char *name) const { return graph.Mask(name); }

  MatrixGraphBuffers &Buffers() { return graph.Buffers(); }

  int OverlapCount() const { return graph.OverlapCount(); }

  int OverlapCountGeometry() const { return graph.OverlapCountGeometry(); }

  /// functions from IDoF to provide functions for Vector, Matrix, etc.
  ///===============================================================================================
  IDoF &GetDoF() const { return graph.GetDoF(); }

  // TODO: Remove
  const IDoF &GetShapeDoF() const { return graph.GetShapeDoF(); }

  std::string DoFsName() const;

  short NumberOfDoFs() const;

  short NumberOfComponents() const;

  short NumberOfNodalPoints(const Cell &c) const;

  short NumberOfNodalPoints(int n, const Cell &c) const;

  short NumberOfNodalPoints(const int sDeg, const int tDeg) const;

  std::vector<Point> GetNodalPoints(const Cell &c) const;

  std::vector<Point> GetNodalPoints(int n, const Cell &c) const;

  std::vector<short> DoFSizesAtNodalPoints(const Cell &c) const;

  std::vector<short> DoFSizesAtNodalPoints(int n, const Cell &c) const;

  short NumberOfNodalPointsOnFace(const Cell &c, int faceId) const;

  short NumberOfNodalPointsOnFace(int n, const Cell &c, int faceId) const;

  short IdOfNodalPointOnFace(const Cell &c, int faceId, int k) const;

  short IdOfNodalPointOnFace(int n, const Cell &c, int faceId, int k) const;

  std::vector<Point> GetNodalPointsOnFace(const Cell &c, int faceId) const;

  std::vector<Point> GetNodalPointsOnFace(int n, const Cell &c, int faceId) const;

  short NumberOfNodalPointsOnEdge(const Cell &c, int edgeId) const;

  short NumberOfNodalPointsOnEdge(int n, const Cell &c, int edgeId) const;

  short IdOfNodalPointOnEdge(const Cell &c, int edgeId, int k) const;

  short IdOfNodalPointOnEdge(int n, const Cell &c, int edgeId, int k) const;

  std::vector<Point> GetNodalPointsOnEdge(const Cell &c, int edgeId) const;

  std::vector<Point> GetNodalPointsOnEdge(int n, const Cell &c, int edgeId) const;
};


#endif // of #ifndef VECTORMATRIXBASE_HPP
