#include "IMatrixGraph.hpp"
#include "MixedDoF.hpp"

void IMatrixGraph::InitIndexData() {
  _rows.InitIndices(singleEntries);
  int nR = _rows.size();
  index = new int[nR + 1];
  diag = new int[nR + 1];
  n = new int[nR];
  column = new int[_rows.NumberOfEntriesTotal() + 1];
  matentry = new int[_rows.NumberOfEntriesTotal() + 1];
  if (singleEntries) singleentry = new bool[_rows.NumberOfEntriesTotal() + 1];
  vector<row> r(nR);
  for (row ri = rows(); ri != rows_end(); ++ri) {
    r[ri.Id()] = ri;
    n[ri.Id()] = ri.n();
  }
  int d = 0;
  index[0] = 0;
  diag[0] = d;
  matentry[d] = 0;
  for (int i = 0; i < nR; ++i) {
    index[i + 1] = index[i] + n[i];
    column[d] = i;
    matentry[d] = r[i].GetEntry();
    ++d;
    double t;
    if (singleEntries && r[i]().isTimeDep()) t = r[i]().t();
    for (entry e = r[i].entries(); e != r[i].entries_end(); ++e) {
      column[d] = e.Id();
      matentry[d] = e.GetEntry();
      if (singleEntries) singleentry[d] = singleEntries && (e().t() != t);
      ++d;
    }
    diag[i + 1] = d;
  }
  N = index[nR];
}

void IMatrixGraph::AddCell(const Cell &c, int depth) {
  vector<Point> sp = dof->GetStoragePoints(c);
  vector<short> allocSizes = dof->AllocationSizesAtStoragePoints(c);
  for (int i = 0; i < sp.size(); ++i)
    _rows.Insert(sp[i], allocSizes[i]);

  switch (depth) {
  case 0:
    return;
  case 1: {
    for (int i = 1; i < sp.size(); ++i)
      for (int j = 0; j < i; ++j)
        _rows.AddEntry(sp[i], sp[j]);
    break;
  }
  default:
    THROW("Only depth 0 and 1 are implemented!")
  }
}

void IMatrixGraph::AddCells(int depth) {
  for (cell c = cells(); c != cells_end(); ++c)
    AddCell(*c, depth);
}

void IMatrixGraph::AddOverlapCells(int depth) {
  for (cell c = overlap(); c != overlap_end(); ++c)
    AddCell(*c, depth);
}

void IMatrixGraph::Init() {
  _rows.InsertInfty(dof->n_infty());
  if (parallel()) { InitParallel(); }
  if (identify()) { InitIdentify(); }
  InitIndexData();
}

Buffer &operator>>(Buffer &buffer, RowIDCommunicationPairHelper &pair) {
  buffer >> pair.ReceivingCellMid >> pair.ReceivingRowIndex >> pair.DofSize
      >> pair.IsReceivingMaster >> pair.master;
  return buffer;
}

Buffer &operator<<(Buffer &buffer, const RowIDCommunicationPairHelper &pair) {
  buffer << pair.ReceivingCellMid << pair.ReceivingRowIndex << pair.DofSize
         << pair.IsReceivingMaster << pair.master;
  return buffer;
}

void IMatrixGraph::InitRowIDCommunicationList() const {
  Comm.clear();
  ExchangeBuffer exBuf(CommSplit());
  size_t nrOfMaster = 0;
  for (procset p = GetProcSets().procsets(); p != GetProcSets().procsets_end(); ++p) {
    bool isMaster = p.master() == PPM->Proc();
    auto R = (*this).find_row(p());
    int id = R.Id();
    for (int k = 0; k < p.size(); ++k) {
      exBuf.Send(p[k]) << RowIDCommunicationPairHelper{p(), Index(id), n[id], isMaster, p.master()};
    }
  }
  exBuf.Communicate();
  for (int q = 0; q < PPM->Size(CommSplit()); ++q) {
    if (q == PPM->Proc(CommSplit())) continue;
    while (exBuf.Receive(q).size() < exBuf.Receive(q).Size()) {
      RowIDCommunicationPairHelper pairhelp;
      exBuf.Receive(q) >> pairhelp;
      Comm.push_back(RowIDCommunicationPair{q, pairhelp.ReceivingRowIndex, pairhelp.DofSize,
                                            index[(*this).find_row(pairhelp.ReceivingCellMid).Id()],
                                            pairhelp.IsReceivingMaster,
                                            pairhelp.master == PPM->Proc()});
    }
  }
}

void IMatrixGraph::InitParallel() {
  for (row r = rows(); r != rows_end(); ++r) {
    procset p = meshAndOverlapProcSets.Find(r());
    if (p != meshAndOverlapProcSets.End()) procSets.Copy(p);
  }
  for (cell c = cells(); c != cells_end(); ++c)
    SetProcSetsCell(*c);
  for (cell c = overlap(); c != overlap_end(); ++c)
    SetProcSetsOverlapCell(*c);
  if (dof->n_infty()) procSets.AddInfty();

  ExchangeBuffer E(CommSplit());
  for (procset p = procSets.procsets(); p != procSets.procsets_end(); ++p) {
    row r = find_row(p());
    if (r == row(_rows.end())) continue;
    if (r.n() == 1) continue;
    for (unsigned int j = 0; j < p.size(); ++j) {
      int q = p[j];
      if (q == PPM->Proc(CommSplit())) continue;
      E.Send(q) << p() << short(r.n());
    }
  }
  E.Communicate();
  for (short q = 0; q < PPM->Size(CommSplit()); ++q) {
    while (E.Receive(q).size() < E.ReceiveSize(q)) {
      Point z;
      short n;
      E.Receive(q) >> z >> n;
      Rows::RowIterator r = _rows.find(z);
      if (r != _rows.end()) r->second.NumberOfDofs(n);
    }
  }
}

void IMatrixGraph::InitIdentify() {
  for (identifyset is = mesh.identifysets(); is != mesh.identifysets_end(); ++is) {
    row r = find_row(is());
    if (r != rows_end()) identifySets.Insert(is);
  }
  for (cell c = cells(); c != cells_end(); ++c)
    IdentifyCell(*c);
  for (cell c = overlap(); c != overlap_end(); ++c)
    IdentifyOverlapCell(*c);
  for (identifyset is = mesh.identifysets(); is != mesh.identifysets_end(); ++is) {
    row r = find_row(is());
    if (r != rows_end()) identifySets.Insert(is);
  }
  for (identifyset is = identifySets.identifysets(); is != identifySets.identifysets_end(); ++is) {
    for (unsigned int i = 0; i < is.size(); ++i) {
      procset p = procSets.find_procset(is[i]);
      if (p == procSets.procsets_end()) continue;
      for (unsigned int j = 0; j < is.size(); ++j)
        procSets.Add(is[j], p);
    }
  }
  for (identifyset is = identifySets.identifysets(); is != identifySets.identifysets_end(); ++is) {
    int n = _rows.RowNumberOfDofs(is());
    for (unsigned int i = 0; i < is.size(); ++i)
      _rows.Insert(is[i], n);
    procset p = procSets.find_procset(is());
    if (p == procSets.procsets_end()) continue;
    for (unsigned int i = 0; i < is.size(); ++i)
      procSets.Copy(p, is[i]);
  }
  for (identifyset is = mesh.identifysets(); is != mesh.identifysets_end(); ++is) {
    row r = find_row(is());
    if (r == rows_end()) continue;
    identifySets.Insert(is);
    if (is.master()) continue;
    for (int i = 0; i < is.size(); ++i) {
      row rr = find_row(is[i]);
      if (rr == rows_end()) continue;
      for (entry e = rr.entries(); e != rr.entries_end(); ++e)
        _rows.AddEntry(r(), e());
    }
  }
  for (row r = rows(); r != rows_end(); ++r) {
    for (entry e = r.entries(); e != r.entries_end(); ++e) {
      identifyset is = identifySets.find_identifyset(e());
      if (is == identifySets.identifysets_end()) continue;
      if (is.master()) continue;
      row r0 = find_row(is[0]);
      _rows.AddEntry(r(), r0());
    }
  }
  for (identifyset is = identifySets.identifysets(); is != identifySets.identifysets_end(); ++is) {
    if (is.master()) continue;
    row r2 = find_row(is());
    row r0 = find_row(is[0]);
    for (entry e = r2.entries(); e != r2.entries_end(); ++e)
      _rows.AddEntry(e(), r0());
    _rows.AddEntry(r2(), r0());
  }
}

void IMatrixGraph::SetProcSetsCell(const Cell &c) {
  procset pc = meshAndOverlapProcSets.Find(c());
  // set procsets on faces
  for (int i = 0; i < c.Faces(); ++i) {
    procset pf = meshAndOverlapProcSets.Find(c.Face(i));
    if (pf == meshAndOverlapProcSets.End()) continue;
    vector<Point> z = dof->GetStoragePoints(c);
    for (int l = 0; l < c.FaceEdges(i); ++l) {
      procset pfe = meshAndOverlapProcSets.Find(c.FaceEdge(i, l));
      if (pfe == meshAndOverlapProcSets.End()) continue;
      for (int k = 0; k < dof->NumberOfStoragePointsOnEdge(c, c.faceedge(i, l)); ++k) {
        int nk = dof->IdOfStoragePointOnEdge(c, c.faceedge(i, l), k);
        procSets.Add(z[nk], pfe);
      }
    }
    for (int k = 0; k < dof->NumberOfStoragePointsOnFace(c, i); ++k) {
      int j = dof->IdOfStoragePointOnFace(c, i, k);
      procSets.Add(z[j], pf);
    }
  }
  // set procsets on edges
  for (int i = 0; i < c.Edges(); ++i) {
    procset p = meshAndOverlapProcSets.Find(c.Edge(i));
    if (p == meshAndOverlapProcSets.End()) continue;
    vector<Point> z = dof->GetStoragePoints(c);
    for (int k = 0; k < dof->NumberOfStoragePointsOnEdge(c, i); ++k) {
      int nk = dof->IdOfStoragePointOnEdge(c, i, k);
      procSets.Add(z[nk], p);
    }
  }
  // set procsets in cell
  if (pc != meshAndOverlapProcSets.End()) {
    vector<Point> z = dof->GetStoragePoints(c);
    for (int i = 0; i < z.size(); ++i) {
      procSets.Append(z[i], pc);
    }
  }
}

void IMatrixGraph::SetProcSetsOverlapCell(const Cell &c) {
  procset p = meshAndOverlapProcSets.Find(c());
  Assert(p != meshAndOverlapProcSets.End());

  // set procsets on faces
  for (int i = 0; i < c.Faces(); ++i) {
    procset pf = meshAndOverlapProcSets.Find(c.Face(i));
    if (pf == meshAndOverlapProcSets.End()) continue;
    vector<Point> z = dof->GetStoragePoints(c);
    for (int l = 0; l < c.FaceEdges(i); ++l) {
      procset pfe = meshAndOverlapProcSets.Find(c.FaceEdge(i, l));
      if (pfe == meshAndOverlapProcSets.End()) continue;
      for (int k = 0; k < dof->NumberOfStoragePointsOnEdge(c, c.faceedge(i, l)); ++k) {
        int nk = dof->IdOfStoragePointOnEdge(c, c.faceedge(i, l), k);
        procSets.Append(z[nk], pfe);
      }
    }
    for (int k = 0; k < dof->NumberOfStoragePointsOnFace(c, i); ++k) {
      int j = dof->IdOfStoragePointOnFace(c, i, k);
      procSets.Append(z[j], pf);
    }
  }
  // set procsets in cell
  vector<Point> z = dof->GetStoragePoints(c);
  for (int i = 0; i < z.size(); ++i) {
    procSets.Append(z[i], p);
  }
}

void IMatrixGraph::IdentifyCell(const Cell &c) {
  for (int i = 0; i < c.Faces(); ++i) {
    identifyset is = mesh.find_identifyset(c.Face(i));
    if (is == mesh.identifysets_end()) continue;
    int mode = mesh.BoundaryFacePart(c.Face(i));
    vector<Point> z = dof->GetStoragePoints(c);
    for (int k = 0; k < dof->NumberOfStoragePointsOnFace(c, i); ++k) {
      int j = dof->IdOfStoragePointOnFace(c, i, k);
      identifySets.Identify(z[j], mode);
      procset p = meshAndOverlapProcSets.Find(z[j]);
      if (p == meshAndOverlapProcSets.End()) continue;
      procSets.Add(z[j], p);
    }
    procset p = meshAndOverlapProcSets.Find(is());
    if (p == meshAndOverlapProcSets.End()) continue;
    for (int k = 0; k < dof->NumberOfStoragePointsOnFace(c, i); ++k) {
      int j = dof->IdOfStoragePointOnFace(c, i, k);
      procSets.Add(z[j], p);
    }
  }
}

void IMatrixGraph::IdentifyOverlapCell(const Cell &c) {
  for (int i = 0; i < c.Faces(); ++i) {
    identifyset is = mesh.find_identifyset(c.Face(i));
    if (is == mesh.identifysets_end()) continue;
    int mode = mesh.BoundaryFacePart(c.Face(i));
    vector<Point> z = dof->GetStoragePoints(c);
    for (int k = 0; k < dof->NumberOfStoragePointsOnFace(c, i); ++k) {
      int j = dof->IdOfStoragePointOnFace(c, i, k);
      identifySets.Identify(z[j], mode);
      procset p = meshAndOverlapProcSets.Find(z[j]);
      if (p == meshAndOverlapProcSets.End()) continue;
      procSets.Append(z[j], p);
    }
    procset p = meshAndOverlapProcSets.Find(is());
    if (p == meshAndOverlapProcSets.End()) continue;
    for (int k = 0; k < dof->NumberOfStoragePointsOnFace(c, i); ++k) {
      int j = dof->IdOfStoragePointOnFace(c, i, k);
      procSets.Append(z[j], p);
    }
  }
}

IMatrixGraph::~IMatrixGraph() { Destruct(); }

void IMatrixGraph::Destruct() {
  if (n) delete[] n;
  if (column) delete[] column;
  if (diag) delete[] diag;
  if (index) delete[] index;
  if (matentry) delete[] matentry;
  if (singleentry) delete[] singleentry;
}

int IMatrixGraph::pSize() const {
  int cnt = 0;
  for (const auto &[p, r] : _rows)
    if (procSets.master(p)) cnt += r.NumberOfDofs();
  return PPM->SumOnCommSplit(cnt, CommSplit());
}

size_t IMatrixGraph::pMatrixSize() const {
  size_t cnt = 0;
  for (const auto &[p_r, r] : _rows) {
    if (procSets.master(p_r)) {
      cnt += r.NumberOfDofs() * r.NumberOfDofs();
      for (auto const &[p_e, e] : r.Entries())
        cnt += 2 * r.NumberOfDofs() * Dof(e.Id());
    } else {
      procset p = meshAndOverlapProcSets.Find(p_r);
      for (auto const &[p_e, e] : r.Entries()) {
        if (procSets.master(p_e)) cnt += 2 * r.NumberOfDofs() * Dof(e.Id());
        else {
          procset pe = meshAndOverlapProcSets.Find(p_e);
          if (p[0] != pe[0]) cnt += 2 * r.NumberOfDofs() * Dof(e.Id());
        }
      }
    }
  }
  return PPM->SumOnCommSplit(cnt, CommSplit());
}

bool IMatrixGraph::SingleEntry(int i) const {
  if (singleEntries) return singleentry[i];
  return false;
}

const SubVectorMask &IMatrixGraph::Mask(const char *name) const { return dof->GetSubVector(name); }

void IMatrixGraph::PrintInfo() const {
  MatrixGraphMemoryInfo<double> vectorMemInfo = VectorMemoryInfo();
  MatrixGraphMemoryInfo<double> matrixMemInfo = MatrixMemoryInfo();
  mout.PrintInfo("Matrix Graph", verbose, PrintInfoEntry("DoF", dof->Name()),
#ifndef USE_SPACETIME
                 PrintInfoEntry("Steps", GetMesh().steps()),
#endif
                 PrintInfoEntry("Problem size", vectorMemInfo.totalSize),
                 PrintInfoEntry("Vector Mem in MB", vectorMemInfo.memTotal_MB, 3),
                 PrintInfoEntry("Vector Mem max on proc in MB", vectorMemInfo.memMaxOnProc_MB, 3),
                 PrintInfoEntry("Vector Mem max on proc in GB", vectorMemInfo.memMaxOnProc_GB, 3),
                 PrintInfoEntry("Vector Total Memory in MB", vectorMemInfo.memTotal_MB, 3),
                 PrintInfoEntry("Vector Total Memory in GB", vectorMemInfo.memTotal_GB, 3),
                 PrintInfoEntry("Vector Memory balancing factor", vectorMemInfo.memBalance, 3),
                 PrintInfoEntry("Matrix size", matrixMemInfo.totalSize, 2),
                 PrintInfoEntry("Matrix Mem in MB", matrixMemInfo.memTotal_MB, 2),
                 PrintInfoEntry("Matrix Mem max on proc in MB", matrixMemInfo.memMaxOnProc_MB, 2),
                 PrintInfoEntry("Matrix Mem max on proc in GB", matrixMemInfo.memMaxOnProc_GB, 2),
                 PrintInfoEntry("Matrix Total Memory in MB", matrixMemInfo.memTotal_MB, 2),
                 PrintInfoEntry("Matrix Total Memory in GB", matrixMemInfo.memTotal_GB, 2),
                 PrintInfoEntry("Matrix Memory balancing factor", matrixMemInfo.memBalance, 2));
}

MatrixGraphMemoryInfo<double> IMatrixGraph::VectorMemoryInfo() const {
  return MatrixGraphMemoryInfo<double>(pSize(), size(), CommSplit(), "VectorMemoryInfo");
}

void IMatrixGraph::PrintVectorMemoryInfo() const { VectorMemoryInfo().PrintInfo(); }

MatrixGraphMemoryInfo<double> IMatrixGraph::MatrixMemoryInfo() const {
  return MatrixGraphMemoryInfo<double>(pMatrixSize(), Size(), CommSplit(), "MatrixMemoryInfo");
}

void IMatrixGraph::PrintMatrixMemoryInfo() const { MatrixMemoryInfo().PrintInfo(); }

template<typename S>
LogTextStream<S> &operator<<(LogTextStream<S> &s, const IMatrixGraph &G) {
  s << "DoFs" << endl;
  for (int i = 0; i < G._rows.size(); ++i)
    s << G.Dof(i) << " ";
  s << endl << "index" << endl;
  for (int i = 0; i <= G._rows.size(); ++i)
    s << G.Index(i) << " ";
  s << endl << "diag" << endl;
  for (int i = 0; i <= G._rows.size(); ++i)
    s << G.Diag(i) << " ";
  s << endl << "column" << endl;
  for (int i = 0; i < G._rows.NumberOfEntriesTotal(); ++i)
    s << G.Column(i) << " ";
  s << endl << "entry" << endl;
  for (int i = 0; i < G._rows.NumberOfEntriesTotal(); ++i)
    s << G.Entry(i) << " ";
  return s << endl
           << "Rows: " << endl
           << G._rows << "ProcSets: " << G.procSets.size() << endl
           << G.procSets.ref() << "IdentifySets: " << endl
           << G.identifySets.size() << endl
           << G.identifySets.ref();
}

rows::rows(const IMatrixGraph &g, const Cell &c) {
  if (g.IsSpaceTime()) {
    resize(1);
    (*this)[0] = g.find_row(c());
    return;
  }
  vector<Point> z = g.GetDoF().GetStoragePoints(c);
  resize(z.size());
  for (int i = 0; i < size(); ++i) {
    (*this)[i] = g.find_row(z[i]);
  }
}

rows::rows(int n, const IMatrixGraph &g, const Cell &c) {
  if (g.IsSpaceTime()) {
    resize(1);
    (*this)[0] = g.find_row(c());
    return;
  }
  vector<Point> z = MixedDoF::Cast(g.GetShapeDoF()).GetStoragePoints(n, c);
  resize(z.size());
  for (int i = 0; i < size(); ++i) {
    (*this)[i] = g.find_row(z[i]);
  }
}

rows::rows(const IMatrixGraph &g, const Cell &c, int face) {
  vector<Point> z = g.GetDoF().GetStoragePoints(c);
  resize(g.GetDoF().NumberOfStoragePointsOnFace(c, face));
  for (int i = 0; i < size(); ++i) {
    (*this)[i] = g.find_row(z[g.GetDoF().IdOfStoragePointOnFace(c, face, i)]);
  }
}

rows::rows(const IMatrixGraph &g, const Cell &c, int face, int n) {
  const MixedDoF &mixedDoF = MixedDoF::Cast(g.GetShapeDoF());
  vector<Point> z = mixedDoF.GetStoragePoints(n, c);
  resize(mixedDoF.NumberOfStoragePointsOnFace(n, c, face));
  for (int i = 0; i < size(); ++i) {
    (*this)[i] = g.find_row(z[mixedDoF.IdOfStoragePointOnFace(n, c, face, i)]);
  }
}
