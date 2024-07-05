#ifndef IMATRIXGRAPH_HPP
#define IMATRIXGRAPH_HPP

#include "ExchangeBuffer.hpp"
#include "IDoF.hpp"
#include "Mesh.hpp"
#include "RowIterators.hpp" //TODO: This should be Row.hpp

class MatrixGraphBuffers {
protected:
  ExchangeBuffer APE;
  ExchangeBuffer CPE;
  ExchangeBuffer AIE;
  ExchangeBuffer CIE;
  ExchangeBuffer AVE;
  ExchangeBuffer MAE;
  ExchangeBuffer MIE;
  ExchangeBuffer DPE;
public:
  explicit MatrixGraphBuffers(int commSplit = 0) :
      APE(commSplit), CPE(commSplit), AIE(commSplit), CIE(commSplit), AVE(commSplit),
      MAE(commSplit), MIE(commSplit), DPE(commSplit) {}

  ExchangeBuffer &AccumulateParallelBuffer() { return APE.Rewind(); }

  ExchangeBuffer &CollectParallelBuffer() { return CPE.Rewind(); }

  ExchangeBuffer &AccumulateIdentifyBuffer() { return AIE.Rewind(); }

  ExchangeBuffer &CollectIdentifyBuffer() { return CIE.Rewind(); }

  ExchangeBuffer &DirichletParallelBuffer() { return DPE.Rewind(); }

  ExchangeBuffer &AverageParallelBuffer() { return AVE.Rewind(); }

  ExchangeBuffer &AccumulateMatrixIdentifyBuffer() { return MIE.Rewind(); }

  ExchangeBuffer &AccumulateMatrixBuffer() { return MAE.Rewind(); }
};

struct RowIDCommunicationPair {
  int ReceivingProc{};
  int ReceivingRowIndex{};
  int DofSize{};
  int SendingRowIndex{};
  bool IsReceivingMaster{};
  bool IsSendingMaster{};
};

struct RowIDCommunicationPairHelper {
  Point ReceivingCellMid{};
  int ReceivingRowIndex{};
  int DofSize{};
  bool IsReceivingMaster{};
  intproc master{};
};

using RowIdCommunicationVector = std::vector<RowIDCommunicationPair>;

template<typename T>
struct MatrixGraphMemoryInfo {
  std::string name;
  size_t totalSize;
  size_t sizeOnProc;
  size_t memOnProc_MB;
  size_t memOnProc_GB;
  size_t sizeMaxOnProc;
  size_t memMaxOnProc_MB;
  size_t memMaxOnProc_GB;
  size_t memTotal_MB;
  size_t memTotal_GB;
  double memBalance;

  MatrixGraphMemoryInfo(size_t totalSize, size_t sizeOnProc, int commSplit,
                        std::string name = "MemoryInfo") :
      name(name), sizeOnProc(sizeOnProc), totalSize(totalSize) {
    sizeMaxOnProc = PPM->Max(sizeOnProc);
    memOnProc_MB = MAX((sizeOnProc * sizeof(T)) / (1024 * 1024), 1);
    memOnProc_GB = memOnProc_MB / 1024;
    memMaxOnProc_MB = PPM->Max(memOnProc_MB, commSplit);
    memMaxOnProc_GB = memMaxOnProc_MB / 1024;
    memTotal_MB = PPM->SumOnCommSplit(memOnProc_MB, commSplit);
    memTotal_GB = memTotal_MB / 1024;
    memBalance = memMaxOnProc_MB / double(PPM->Min(memOnProc_MB, commSplit));
  }

  void PrintInfo() {
    mout.PrintInfo(name, 1, PrintInfoEntry("Problem size", totalSize),
                   PrintInfoEntry("Max Problem size on Proc", sizeMaxOnProc),
                   PrintInfoEntry("Mem max on proc in MB", memMaxOnProc_MB),
                   PrintInfoEntry("Mem max on proc in GB", memMaxOnProc_GB),
                   PrintInfoEntry("Total Memory in MB", memTotal_MB),
                   PrintInfoEntry("Total Memory in GB", memTotal_GB),
                   PrintInfoEntry("Memory balancing factor", memBalance));
  }
};

class IMatrixGraph {
protected:
  int verbose = 0;

  mutable MatrixGraphBuffers buffers;

  std::unique_ptr<IDoF> dof = nullptr;
  const Mesh &mesh;
  const LevelPair level;


  Rows _rows;
  IdentifySets identifySets;
  ProcSets procSets;
  Cells overlapCells{};
  ProcSets meshAndOverlapProcSets{};
  Faces meshAndOverlapFaces{};


  bool singleEntries = false; // TODO: Move to Rows?
  int N = -1;
  // index for accessing the subentries of Vector
  int *index = nullptr;
  int *diag = nullptr;

  int *column = nullptr;
  int *matentry = nullptr;
  int *n = nullptr;
  bool *singleentry = nullptr;
  mutable RowIdCommunicationVector Comm;

  void InitRowIDCommunicationList() const;

  void InitIndexData();

  const RowIdCommunicationVector &GetRowIdCommunicationList() const {
    InitRowIDCommunicationList();
    return Comm;
  };

  void AddCell(const Cell &c, int depth = 1);

  void AddCells(int depth = 1);

  void AddOverlapCells(int depth = 1);

  virtual void Init();

  void InitParallel();

  void InitIdentify();

  virtual void SetProcSetsCell(const Cell &);

  virtual void SetProcSetsOverlapCell(const Cell &);

  virtual void IdentifyCell(const Cell &);

  virtual void IdentifyOverlapCell(const Cell &);

  IMatrixGraph(const Mesh &mesh, std::unique_ptr<IDoF> _dof, LevelPair level, bool single = false) :
      dof(std::move(_dof)), mesh(mesh), singleEntries(single), buffers(mesh.CommSplit()),
      level(level) {

    // This information is copied from mesh, since overlap manipulates both of these containers
    // Would be nice if after some major refactoring of MatrixGraph classes this disappears again
    meshAndOverlapFaces = mesh.meshFaces;
    meshAndOverlapProcSets.SetCommSplit(CommSplit());
    for (procset p = mesh.procsets(); p != mesh.procsets_end(); p++) {
      meshAndOverlapProcSets.InsertFront(p(), *p);
    }

    procSets.SetCommSplit(CommSplit());
    Config::Get("MatrixGraphVerbose", verbose);
  }

  IMatrixGraph(const Mesh &mesh, std::unique_ptr<IDoF> _dof, bool single = false) :
      IMatrixGraph(mesh, std::move(_dof), mesh.Level(), single) {}
public:
  virtual ~IMatrixGraph();

  virtual void Destruct();

  LevelPair GetLevel() const { return level; }

  const IDoF &GetDoF() const { return *dof; }

  IDoF &GetDoF() { return *dof; }

  /// usually vector and shape DoFs coincide, if not override this function (cf. HybridMatrixGraph)
  virtual const IDoF &GetShapeDoF() const { return *dof; }

  const Mesh &GetMesh() const { return mesh; }

  int CommSplit() const { return GetMesh().CommSplit(); }

  MatrixGraphBuffers &Buffers() const { return buffers; }

  Rows &GetRows() { return _rows; }

  const Rows &GetRows() const { return _rows; }

  IdentifySets &GetIdentifySets() { return identifySets; }

  const IdentifySets &GetIdentifySets() const { return identifySets; }

  ProcSets &GetProcSets() { return procSets; }

  const ProcSets &GetProcSets() const { return procSets; }

  size_t Size() const { return _rows.NumberOfDofsTotal(); }

  int size() const { return N; }

  int pSize() const;

  size_t pMatrixSize() const;

  int rowsize() const { return _rows.size(); }

  int Index(int i) const { return index[i]; }

  int Column(int i) const { return column[i]; }

  int Diag(int i) const { return diag[i]; }

  int Entry(int i) const { return matentry[i]; }

  int Dof(int i) const { return n[i]; }

  bool SingleEntry(int i) const;

  const SubVectorMask &Mask(const char *name) const;

  void PrintInfo() const;

  MatrixGraphMemoryInfo<double> VectorMemoryInfo() const;

  void PrintVectorMemoryInfo() const;

  MatrixGraphMemoryInfo<double> MatrixMemoryInfo() const;

  void PrintMatrixMemoryInfo() const;

  virtual void update(){};

  double t(int k) const { return GetMesh().t(k); }

  bool onBnd(const Point &x) const { return GetMesh().onBnd(x); }

  int BndPart(const Point &x) const { return GetMesh().BoundaryFacePart(x); }

  bnd_face bnd_faces() const { return GetMesh().bnd_faces(); }

  bnd_face bnd_faces_end() const { return GetMesh().bnd_faces_end(); }

  bnd_face find_bnd_face(const Point &z) const { return GetMesh().find_bnd_face(z); }

  face find_face(const Point &center) const { return meshAndOverlapFaces.Find(center); }

  void InsertOverlapCell(const Cell *C) {
    if (!GetMesh().IsSTMesh()) {
      auto *oCell = CreateCell(C->Type(), C->Subdomain(), C->AsVector());
#ifdef USE_DATAMESH
      oCell->SetData(C->GetData());
#endif
      overlapCells.Insert(C->Center(), oCell);
    } else {
      auto *TC =
          CreateCell(C->Type(), C->Subdomain(), C->SpaceCell().AsVector(), C->min(), C->max());
      overlapCells.Insert((*TC)(), TC);
    }
  }

  void InsertOverlapFace(const Point &faceCenter, const Point &cellCenter) {
    auto f = meshAndOverlapFaces.Find(faceCenter);
    if (f == meshAndOverlapFaces.End()) {
      meshAndOverlapFaces.Insert(faceCenter, Face(cellCenter));
    } else if (f->second.Right() == Infty) {
      if (f->second.Left() != cellCenter) {
        meshAndOverlapFaces.Replace(faceCenter, Face(f->second.Left(), cellCenter));
      }
    }
  }

  face faces_end() const { return {meshAndOverlapFaces.End()}; }

  face faces() const { return {meshAndOverlapFaces.Begin()}; }

  vertex vertices() const { return GetMesh().vertices(); }

  vertex vertices_end() const { return GetMesh().vertices_end(); }

  vertex find_vertex(const Point &p) const { return GetMesh().find_vertex(p); }

  cell find_neighbour_cell(const cell &c, int f) const { return find_neighbour_cell(*c, f); }

  cell find_neighbour_cell(const Cell &c, int f) const {
    face ff = find_face(c.Face(f));
    if constexpr (DebugLevel > 0) {
      if (ff == faces_end()) { THROW("find_neighbour_cell, Face not found!") }
    }
    if (c() != ff.Right()) return find_cell_or_overlap_cell(ff.Right());
    else return find_cell_or_overlap_cell(ff.Left());
  }

  int find_neighbour_face_id(const Point &f_c, const cell &cf) const {
    for (int f1 = 0; f1 < cf.Faces(); ++f1) {
      if (f_c == cf.Face(f1)) return f1;
    }
    THROW("find_neighbour_face_id Face not found!")
  }

  Point neighbour_center(const Cell &c, int i) const {
    face f = find_face(c.Face(i));
    if (f.Right() == c()) return f.Left();
    else return f.Right();
  }

  edge edges() const { return GetMesh().edges(); }

  edge edges_end() const { return GetMesh().edges_end(); }

  edge find_edge(const Point &center) const { return GetMesh().find_edge(center); }

  identifyset find_identifyset(const Point &z) const { return {identifySets.Find(z)}; }

  identifyset identifysets() const { return {GetMesh().identifysets()}; }

  identifyset identifysets_end() const { return {GetMesh().identifysets_end()}; }

  procset procsets() const { return {procSets.Begin()}; }

  procset procsets_end() const { return {procSets.End()}; }

  procset find_procset(const Point &z) const { return {procSets.Find(z)}; }

  cell find_cell_or_overlap_cell(const Point &z) const {
    cell c = find_cell(z);
    if (c != cells_end()) return c;
    c = find_overlap_cell(z);
    if (c != overlap_end()) return c;

    std::cout << "No cell found for Point:" << z << "  Proc: " << PPM->Proc(CommSplit()) << endl;
    THROW("Error: no cell found")
  }

  bool has_cell_or_overlap_cell(const Point &z) const {
    cell c = find_cell(z);
    if (c != cells_end()) return true;
    c = find_overlap_cell(z);
    if (c != overlap_end()) return true;
    return false;
  }

  cell cells_end() const { return {GetMesh().cells_end()}; }

  cell overlap_end() const { return {overlapCells.End()}; }

  cell find_cell(const Point &center) const { return {GetMesh().find_cell(center)}; }

  bool OnBoundary(const Cell &c) const {
    for (int i = 0; i < c.Faces(); ++i) {
      int part = GetMesh().BoundaryFacePart(c.Face(i));
      if (part != -1) return true;
    }
    return false;
  }

  int OverlapCount() const { return (int)overlapCells.size(); }

  int OverlapCountGeometry() const { return PPM->SumOnCommSplit(OverlapCount(), level.commSplit); }

  bool has_previous_cell(const cell &c) const { return find_previous_cell(c) != c; }

  bool OnBoundary(const Cell &c, int f) const {
#ifdef USE_SPACETIME
    face ff = find_face(c.Face(f));
    return ((ff.Right() == Infty) || (ff.Left() == Infty));
#else
    return onBnd(c.Face(f));
#endif
  }

  cell overlap() const {
    return overlapCells.Begin();
    ;
  }

  cell cells() const { return GetMesh().cells(); }

  cell find_previous_cell(const cell &c) const {
    face f = find_face(c.Face(c.Faces() - 2));
    if (f == faces_end()) { return c; }
    cell c_prev = find_cell_or_overlap_cell(f.Left());
    if (c_prev == c && has_cell_or_overlap_cell(f.Right())) {
      c_prev = find_cell_or_overlap_cell(f.Right());
    }
    if (c_prev == cells_end() || c_prev == overlap_end()) { return c; }

    return c_prev;
  }

  cell find_overlap_cell(const Point &center) const { return {overlapCells.Find(center)}; }

  bool identify() const { return GetMesh().identify(); }

  bool parallel() const { return GetMesh().parallel(); }

  int dim() const { return GetMesh().dim(); }

  row rows() const { return row(_rows.begin()); }

  row rows_end() const { return row(_rows.end()); }

  row find_row(const Point &z) const { return row(_rows.find(z)); }

  virtual bool IsSpaceTime() const { return false; }

  template<typename S>
  friend LogTextStream<S> &operator<<(LogTextStream<S> &s, const IMatrixGraph &G);
};

class rows : public vector<row> {
public:
  rows(const IMatrixGraph &graph, const Cell &c);

  rows(int n, const IMatrixGraph &graph, const Cell &c);

  rows(const IMatrixGraph &graph, const Cell &c, int face);

  rows(const IMatrixGraph &graph, const Cell &c, int face, int n);

  int n(int i) const { return (*this)[i].n(); }
};

#endif // IMATRIXGRAPH_HPP
