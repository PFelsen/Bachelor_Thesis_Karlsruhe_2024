#include "SpaceTimeMatrixGraph.hpp"
#include "AdaptiveSpaceTimeDof.hpp"

void SpaceTimeMatrixGraph::AddNBCells() {
  if (isBlockDiagonal) { return; }
  for (cell c = cells(); c != cells_end(); c++) {
    for (int i = 0; i < c.Faces(); ++i) {
      face f = find_face(c.Face(i));
      if (f.Right() != Infty && f.Left() != Infty) { _rows.AddEntry(f.Left(), f.Right()); }
    }
    if (isDgInTime) continue;

    face f_up = find_face(c.Face(c.Faces() - 1));
    if (f_up.Right() != Infty) {
      cell c_up = find_cell_or_overlap_cell(f_up.Right());
      for (int i = 0; i < c_up.Faces() - 2; ++i) {
        face f = find_face(c_up.Face(i));
        if (f == faces_end()) continue;
        if (f.Right() != Infty) {
          if (f.Left() == c_up()) _rows.AddEntry(c(), f.Right());
          else _rows.AddEntry(c(), f.Left());
        }
      }
    }

    face f_lo = find_face(c.Face(c.Faces() - 2));
    if (f_lo == faces_end()) continue;
    if (f_lo.Right() != Infty) {
      cell c_lo = find_cell_or_overlap_cell(f_lo.Left());
      for (int i = 0; i < c_lo.Faces() - 2; ++i) {
        face f = find_face(c_lo.Face(i));
        if (f == faces_end()) continue;
        if (f.Right() != Infty) {
          if (f.Left() == c_lo()) _rows.AddEntry(c(), f.Right());
          else _rows.AddEntry(c(), f.Left());
        }
      }
    }
  }
}

void SpaceTimeMatrixGraph::InitProcsets() {
  _rows.InsertInfty(this->dof->n_infty());

  if (parallel()) {
    for (row r = rows(); r != rows_end(); ++r) {
      procset p = meshAndOverlapProcSets.Find(r());
      if (p != meshAndOverlapProcSets.End()) procSets.Copy(p);
    }
  }
  InitIndexData();
}

void SpaceTimeMatrixGraph::Init() {

  auto &dof2 = dynamic_cast<AdaptiveSpaceTimeDof &>(*dof);
  dof2.SetMatrixGraph(this);
  LOWEST_POSSIBLE_TIME_DEG = dof2.GetLOWEST_POSSIBLE_TIME_DEG();

  // Todo: This shouldn't be set via Config. Discretization should pick the right Overlap
  std::string overlapName = "STCellsWithCorners";
  Config::Get("Overlap", overlapName);
  if (overlapName == "STCellsWithFaces") {
    STOverlapCellsWithFaces();
  } else if (overlapName == "STCellsWithCorners") {
    STOverlapCellsWithCorners();
  }


  if (cell_deg.empty()) {
    for (cell c = cells(); c != cells_end(); ++c) {
      cell_deg[c()] = defaultDegree;
    }
  }
  max_deg = {0, 0};
  for (cell c = cells(); c != cells_end(); ++c) {
    DegreePair degs = cell_deg[c()];
    max_deg = {max(max_deg.space, degs.space), max(max_deg.time, degs.time)};
  }
  communicate();

  AddCells(0);
  AddOverlapCells(0);
  AddNBCells();
  InitProcsets();
}

SpaceTimeMatrixGraph::SpaceTimeMatrixGraph(const Mesh &stMesh, std::unique_ptr<IDoF> _dofs,
                                           DegreePair defaultDegreePair, bool isDgInTime,
                                           bool single, LevelPair level,
                                           const std::unordered_map<Point, DegreePair> &polyDist,
                                           bool blockDiagonal) :
    IMatrixGraph(stMesh, std::move(_dofs), level, single), isBlockDiagonal(blockDiagonal),
    isDgInTime(isDgInTime), cell_deg(polyDist), defaultDegree(defaultDegreePair) {
  Init();
}

void SpaceTimeMatrixGraph::update() {
  _rows.clear();
  procSets.clear();
  identifySets.clear();
  Destruct();
  buffers = MatrixGraphBuffers();
  Init();
}

void SpaceTimeMatrixGraph::STOverlapCellsWithCorners() {
  int commSplit = CommSplit();
  if (PPM->Size(commSplit) == 1) return;
  std::map<Point, bool> bnd_vert;
  ExchangeBuffer E_cells(CommSplit());
  for (cell c = cells(); c != cells_end(); ++c) {
    for (int fi = 0; fi < c.Faces(); fi++) {
      face f_tmp = find_face(c.Face(fi));
      cell c1 = find_cell(f_tmp.Left());
      cell c2 = find_cell(f_tmp.Right());
      if (c1 == cells_end() || c2 == cells_end()) {
        for (int fci = 0; fci < c.FaceCorners(fi); ++fci) {
          Point facecorner = c.FaceCorner(fi, fci);
          if (bnd_vert.find(facecorner) == bnd_vert.end()) {
            bnd_vert[facecorner] = true;
            ;
          }
        }
      }
    }
  }

  for (cell c = cells(); c != cells_end(); ++c) {
    bool bnd_proc_cell = false;

    for (int i = 0; i < c.Corners(); ++i) {
      if (bnd_vert.find(c.Corner(i)) != bnd_vert.end()) {
        bnd_proc_cell = true;
        break;
      }
    }

    if (bnd_proc_cell) {
      meshAndOverlapProcSets.Add(c(), PPM->Proc(commSplit));
      for (int vi = 0; vi < c.Corners(); ++vi) {
        Point cp = c.Corner(vi);
        procset ps2 = meshAndOverlapProcSets.Find(cp);
        if (ps2 != meshAndOverlapProcSets.End()) { meshAndOverlapProcSets.Append(c(), ps2); }
      }
      ProcSet ps = *meshAndOverlapProcSets.find_procset(c());
      for (short q : ps) {
        if (q == PPM->Proc(commSplit)) continue;
        E_cells.Send(q) << *c;
        E_cells.Send(q) << ps;
      }
    }
  }
  E_cells.Communicate();

  for (int q = 0; q < PPM->Size(commSplit); ++q) {
    if (q == PPM->Proc(commSplit)) continue;
    while (E_cells.Receive(q).size() < E_cells.Receive(q).Size()) {
      Cell *C = nullptr;
      E_cells.Receive(q) >> C;
      ProcSet ps;
      E_cells.Receive(q) >> ps;
      for (int i = 0; i < C->Faces(); ++i) {
        InsertOverlapFace(C->Face(i), (*C)());
      }
      InsertOverlapCell(C);
      meshAndOverlapProcSets.Add((*C)(), q);
      meshAndOverlapProcSets.Insert((*C)(), ps);
      delete C;
    }
  }
  meshAndOverlapProcSets.RemoveSingle();
}

void SpaceTimeMatrixGraph::communicate() {
  ExchangeBuffer E(CommSplit());
  for (procset p = meshAndOverlapProcSets.Begin(); p != meshAndOverlapProcSets.End(); ++p)
    if (find_cell(p()) != cells_end()) {
      for (int j = 0; j < p.size(); ++j) {
        int q = p[j];
        if (q == PPM->Proc(CommSplit())) continue;
        DegreePair degs = cell_deg[p()];
        E.Send(q) << p();
        E.Send(q) << degs.space;
        E.Send(q) << degs.time;
      }
    }
  E.Communicate();
  for (int q = 0; q < PPM->Size(CommSplit()); ++q) {
    while (E.Receive(q).size() < E.ReceiveSize(q)) {
      Point z;
      DegreePair degs{-1, -1};
      E.Receive(q) >> z;
      E.Receive(q) >> degs.space;
      E.Receive(q) >> degs.time;
      SetDegree(z, degs);
    }
  }
  E.ClearBuffers();
  max_deg = {PPM->Max(max_deg.space, CommSplit()), PPM->Max(max_deg.space, CommSplit())};
}

void SpaceTimeMatrixGraph::STOverlapCellsWithFaces() {
  int commSplit = CommSplit();
  if (PPM->Size(commSplit) == 1) return;

  ExchangeBuffer E_cells(CommSplit());
  for (cell c = cells(); c != cells_end(); ++c) {
    for (int fi = 0; fi < c.Faces(); fi++) {
      face f_tmp = find_face(c.Face(fi));
      cell c1 = find_cell(f_tmp.Left());
      procset ps = meshAndOverlapProcSets.Find(f_tmp());
      if (ps == meshAndOverlapProcSets.End()) continue;
      if (c1 != cells_end()) { c1 = find_cell(f_tmp.Right()); }
//      if (c1 != cells_end()) { THROW("Error should find at least one cell of face") }
      meshAndOverlapProcSets.Add(c(), PPM->Proc(commSplit));
      meshAndOverlapProcSets.Append(c(), ps);
      for (short q : *ps) {
        if (q == PPM->Proc(commSplit)) continue;
        E_cells.Send(q) << *c;
        E_cells.Send(q) << *ps;
      }
    }
  }
  E_cells.Communicate();

  for (int q = 0; q < PPM->Size(commSplit); ++q) {
    if (q == PPM->Proc(commSplit)) continue;
    while (E_cells.Receive(q).size() < E_cells.Receive(q).Size()) {
      Cell *C = nullptr;
      E_cells.Receive(q) >> C;
      ProcSet ps;
      E_cells.Receive(q) >> ps;
      for (int i = 0; i < C->Faces(); ++i) {
        InsertOverlapFace(C->Face(i), (*C)());
      }
      InsertOverlapCell(C);
      meshAndOverlapProcSets.Add((*C)(), q);
      meshAndOverlapProcSets.Insert((*C)(), ps);
      delete C;
    }
  }

  for (int i = 0; i < 1; i++) {
    ExchangeBuffer CB(commSplit);
    for (procset p = meshAndOverlapProcSets.Begin(); p != meshAndOverlapProcSets.End(); ++p) {
      for (short q = 0; q < p.size(); ++q) {
        CB.Send(p[q]) << p;
      }
    }
    CB.Communicate();
    for (short q = 0; q < PPM->Size(commSplit); ++q) {
      if (q == PPM->Proc(commSplit)) continue;
      while (CB.Receive(q).size() < CB.Receive(q).Size()) {
        std::pair<Point, ProcSet> receivedProcSetData;
        CB.Receive(q) >> receivedProcSetData;
        meshAndOverlapProcSets.Insert(receivedProcSetData.first, receivedProcSetData.second);
        Assert(receivedProcSetData.second.existselementof(q));
      }
    }
  }
  meshAndOverlapProcSets.RemoveSingle();
}