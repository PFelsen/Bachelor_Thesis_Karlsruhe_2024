#include "DGMatrixGraph.hpp"

#include <vector>

#include "Point.hpp"

void DGMatrixGraph::AddNBCells() {
  for (face f = faces(); f != faces_end(); ++f) {
    const Point &L = f.Left();
    const Point &R = f.Right();
    if (R != Infty) {
      if (L != Infty) _rows.AddEntry(L, R);
    } else {
      identifyset is = mesh.find_identifyset(f());
      if (is == mesh.identifysets_end()) continue;
      face ff = find_face((*is)[0]);
      _rows.AddEntry(L, ff.Left());
    }
  }
}

void DGMatrixGraph::AddCrossPoints() {
  ExchangeBuffer E(CommSplit());
  for (cell c = cells(); c != cells_end(); ++c) {
    for (int i = 0; i < c.Corners(); ++i) {
      procset ps = meshAndOverlapProcSets.Find(c[i]);
      if (ps == meshAndOverlapProcSets.End()) continue;
      if (ps.size() < 3) continue;
      for (int k = 0; k < ps.size(); ++k)
        if (ps[k] != PPM->Proc(CommSplit())) {
          vector<short> z = dof->AllocationSizesAtStoragePoints(*c);
          E.Send(ps[k]) << c() << z[0];
        }
    }
  }
  E.Communicate();
  for (int q = 0; q < PPM->Size(CommSplit()); ++q) {
    if (q == PPM->Proc(CommSplit())) continue;
    while (E.Receive(q).size() < E.Receive(q).Size()) {
      Point z;
      short m;
      E.Receive(q) >> z >> m;
      _rows.Insert(z, m);
    }
  }
}

void DGMatrixGraph::SendProcSets() {
  ExchangeBuffer E(CommSplit());
  for (cell c = cells(); c != cells_end(); ++c) {
    for (int i = 0; i < c.Corners(); ++i) {
      procset ps = meshAndOverlapProcSets.Find(c[i]);
      if (ps == meshAndOverlapProcSets.End()) continue;
      if (ps.size() < 3) continue;
      for (int k = 0; k < ps.size(); ++k)
        procSets.Append(c(), ps[k]);
    }
  }
  for (cell c = cells(); c != cells_end(); ++c) {
    procset ps = procSets.find_procset(c());
    if (ps == procSets.procsets_end()) continue;
    if (ps.size() < 3) continue;
    for (int k = 0; k < ps.size(); ++k) {
      if (ps[k] != PPM->Proc(CommSplit())) {
        E.Send(ps[k]) << c() << short(ps.size());
        for (int i = 0; i < ps.size(); ++i)
          E.Send(ps[k]) << short(ps[i]);
      }
    }
  }
  E.Communicate();
  for (int q = 0; q < PPM->Size(CommSplit()); ++q) {
    if (q == PPM->Proc(CommSplit())) continue;
    while (E.Receive(q).size() < E.Receive(q).Size()) {
      Point z;
      short m;
      E.Receive(q) >> z >> m;
      procSets.Add(z, q);
      for (int i = 0; i < m; ++i) {
        short p;
        E.Receive(q) >> p;
        procSets.Append(z, p);
      }
    }
  }
}

DGMatrixGraph::DGMatrixGraph(const Mesh &mesh, std::unique_ptr<IDoF> _dof) :
    IMatrixGraph(mesh, std::move(_dof)) {
  CellsOverlap_dG1();
  AddCells(0);
  AddOverlapCells(0);
  AddCrossPoints();
  AddNBCells();
  Init();
  SendProcSets();
}

void DGMatrixGraph::CellsOverlap_dG1() {
  int commSplit = CommSplit();
  std::vector<Point> clean_list;

  for (procset p = meshAndOverlapProcSets.Begin(); p != meshAndOverlapProcSets.End(); ++p) {
    if (find_vertex(p()) != vertices_end()) continue;
    if (find_face(p()) != faces_end()) continue;
    if (find_edge(p()) != edges_end()) continue;
    if (p() == Infty) continue;
    clean_list.push_back(p());
  }

  for (const auto &p : clean_list) {
    meshAndOverlapProcSets.Remove(p);
  }

  ExchangeBuffer cellBuffer(commSplit);
  ExchangeBuffer faceBuffer(commSplit);
  ExchangeBuffer edgeBuffer(commSplit);

  for (cell c = cells(); c != cells_end(); ++c) {
    for (int fi = 0; fi < c.Faces(); fi++) {
      face f_tmp = find_face(c.Face(fi));
      procset pf = meshAndOverlapProcSets.Find(f_tmp());
      if (pf != meshAndOverlapProcSets.End()) {
        if (meshAndOverlapProcSets.Find(c()) == meshAndOverlapProcSets.End()) {
          meshAndOverlapProcSets.Add(c(), PPM->Proc(commSplit));
        }

        for (int i = 0; i < pf.size(); ++i)
          meshAndOverlapProcSets.Append(c(), pf[i]);
      }
    }

    procset pc = meshAndOverlapProcSets.Find(c());
    if (pc != meshAndOverlapProcSets.End())
      for (int q = 0; q < pc.size(); ++q) {
        if (pc[q] == PPM->Proc(commSplit)) continue;
        cellBuffer.Send(pc[q]) << c;
        cellBuffer.Send(pc[q]) << short(pc.size() - 1); // send size of procset
        for (int ql = 1; ql < pc.size(); ++ql)          // send procset without master
          cellBuffer.Send(pc[q]) << short(pc[ql]);
      }
  }

  cellBuffer.Communicate();

  for (int q = 0; q < PPM->Size(commSplit); ++q) {
    if (q == PPM->Proc(commSplit)) continue;
    while (cellBuffer.Receive(q).size() < cellBuffer.Receive(q).Size()) {
      int tp;
      short n, sd, q_list_size;
      cellBuffer.Receive(q) >> tp >> sd >> n;
      vector<Point> x(n);
      for (int i = 0; i < n; ++i)
        cellBuffer.Receive(q) >> x[i];

      Cell *C = CreateCell(CELLTYPE(tp), sd, x);
#ifdef USE_DATAMESH
      DataContainer d;
      cellBuffer.Receive(q) >> d;
      C->SetData(d);
#endif

      InsertOverlapCell(C);
      meshAndOverlapProcSets.Add(C->Center(), q);

      for (int i = 0; i < C->Faces(); ++i)
        InsertOverlapFace(C->Face(i), C->Center());

      cellBuffer.Receive(q) >> q_list_size;
      for (int ql = 0; ql < q_list_size; ++ql) {
        short proc;
        cellBuffer.Receive(q) >> proc;
        meshAndOverlapProcSets.Append(C->Center(), proc);
      }
      delete C;
    }
  }
}