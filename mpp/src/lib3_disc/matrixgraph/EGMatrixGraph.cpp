#include "EGMatrixGraph.hpp"

#include <vector>

void EGMatrixGraph::AddNBCells() {
  for (face f = faces(); f != faces_end(); ++f) {
    if (f.Right() != Infty) {
      if (f.Left() != Infty) {
        cell c = find_cell_or_overlap_cell(f.Left());
        cell cf = find_cell_or_overlap_cell(f.Right());
        vector<Point> z = dof->GetStoragePoints(*c);
        vector<Point> zf = dof->GetStoragePoints(*cf);
        for (int i = 0; i < z.size(); ++i)
          for (int j = 0; j < zf.size(); ++j)
            if (z[i] != zf[j]) _rows.AddEntry(z[i], zf[j]);
      }
    }
  }
}

EGMatrixGraph::EGMatrixGraph(const Mesh &mesh, std::unique_ptr<IDoF> _dof, int depth) :
    IMatrixGraph(mesh, std::move(_dof)) {
  CellsOverlapEG();
  AddCells(depth);
  AddOverlapCells(depth);
  AddNBCells();
  Init();
}

void EGMatrixGraph::CellsOverlapEG() {
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

  std::vector<std::reference_wrapper<const Cell>> cellsToSend;
  for (cell c = cells(); c != cells_end(); ++c) {
    for (int fi = 0; fi < c.Faces(); fi++) {
      face f_tmp = find_face(c.Face(fi));
      procset pf = meshAndOverlapProcSets.Find(f_tmp());
      if (pf != meshAndOverlapProcSets.End()) {
        if (meshAndOverlapProcSets.Find(c()) == meshAndOverlapProcSets.End()) {
          cellsToSend.emplace_back(*c);
          break;
        }
      }
    }
  }

  for (const Cell &c : cellsToSend) {
    InsertOverlapCell(&c);
    for (int fi = 0; fi < c.Faces(); fi++) {
      face f_tmp = find_face(c.Face(fi));
      procset pf = meshAndOverlapProcSets.Find(f_tmp());
      if (pf != meshAndOverlapProcSets.End()) {
        if (meshAndOverlapProcSets.Find(c()) == meshAndOverlapProcSets.End()) {
          // setting current proc as master of procset for cell
          meshAndOverlapProcSets.Add(c(), PPM->Proc(commSplit));
        }
        meshAndOverlapProcSets.Append(c(), pf);
      }
    }
    const procset pc = meshAndOverlapProcSets.Find(c());
    Assert(pc != meshAndOverlapProcSets.End());
    for (int j = 0; j < c.Corners(); ++j) {
      if (meshAndOverlapProcSets.Find(c[j]) == meshAndOverlapProcSets.End()) {
        // setting current proc as potential master of procset for corner
        meshAndOverlapProcSets.Add(c.Corner(j), PPM->Proc(commSplit));
      }
      meshAndOverlapProcSets.Append(c.Corner(j), pc);
    }
    for (int f = 0; f < c.Faces(); ++f) {
      Point face_p = c.Face(f);
      procset pf = meshAndOverlapProcSets.Find(face_p);
      if (pf == meshAndOverlapProcSets.End()) {
        meshAndOverlapProcSets.Add(face_p, PPM->Proc(commSplit));
      }
      Assert(meshAndOverlapProcSets.Find(face_p)->second.existselementof(PPM->Proc(commSplit)));
    }
    for (int e = 0; e < c.Edges(); ++e) {
      Point edge_p = c.Edge(e);
      procset pe = meshAndOverlapProcSets.Find(edge_p);
      if (pe == meshAndOverlapProcSets.End()) {
        meshAndOverlapProcSets.Add(edge_p, PPM->Proc(commSplit));
      }
      Assert(meshAndOverlapProcSets.Find(edge_p)->second.existselementof(PPM->Proc(commSplit)));
    }

    for (int q = 0; q < pc.size(); ++q) {
      if (pc[q] == PPM->Proc(commSplit)) continue;
      cellBuffer.Send(pc[q]) << c;
      cellBuffer.Send(pc[q]) << *pc;

      for (int j = 0; j < c.Corners(); ++j) {
        procset corner_pc = meshAndOverlapProcSets.Find(c[j]);
        cellBuffer.Send(pc[q]) << corner_pc;
      }

      for (int f = 0; f < c.Faces(); ++f) {
        Point face_p = c.Face(f);
        procset pf = meshAndOverlapProcSets.Find(face_p);
        faceBuffer.Send(pc[q]) << pf;
      }

      for (int e = 0; e < c.Edges(); ++e) {
        Point edge_p = c.Edge(e);
        procset pe = meshAndOverlapProcSets.Find(edge_p);
        edgeBuffer.Send(pc[q]) << pe;
      }
    }
  }

  cellBuffer.Communicate();
  faceBuffer.Communicate();
  edgeBuffer.Communicate();

  for (int q = 0; q < PPM->Size(commSplit); ++q) {
    if (q == PPM->Proc(commSplit)) continue;
    while (cellBuffer.Receive(q).size() < cellBuffer.Receive(q).Size()) {
      Cell *receivedCell = nullptr;
      cellBuffer.Receive(q) >> receivedCell;
      InsertOverlapCell(receivedCell);
      for (int i = 0; i < receivedCell->Faces(); ++i) {
        InsertOverlapFace(receivedCell->Face(i), receivedCell->Center());
        meshAndOverlapProcSets.Append(receivedCell->Face(i), q);
        meshAndOverlapProcSets.Append(receivedCell->Face(i), PPM->Proc(commSplit));
        procset pf = meshAndOverlapProcSets.Find(receivedCell->Face(i));
        for (int fe = 0; fe < receivedCell->FaceEdges(i); fe++) {
          meshAndOverlapProcSets.Append(receivedCell->FaceEdge(i, fe), pf);
        }
      }
      ProcSet receivedPS_Cell;
      cellBuffer.Receive(q) >> receivedPS_Cell;

      for (int j = 0; j < receivedCell->Corners(); ++j) {
        std::pair<Point, ProcSet> receivedPC_Corner;
        cellBuffer.Receive(q) >> receivedPC_Corner;
        meshAndOverlapProcSets.InsertFront(receivedPC_Corner.first, receivedPC_Corner.second);
        meshAndOverlapProcSets.Insert((*receivedCell)[j], receivedPS_Cell);
      }

      if (receivedPS_Cell[0] == q) {
        meshAndOverlapProcSets.InsertFront(receivedCell->Center(), receivedPS_Cell);
      } else {
        meshAndOverlapProcSets.Insert(receivedCell->Center(), receivedPS_Cell);
      }
      Assert(meshAndOverlapProcSets.Find(receivedCell->Center())
                 ->second.existselementof(PPM->Proc(commSplit)));
      delete receivedCell;
    }

    while (faceBuffer.Receive(q).size() < faceBuffer.Receive(q).Size()) {
      std::pair<Point, ProcSet> facePS;
      faceBuffer.Receive(q) >> facePS;
      meshAndOverlapProcSets.InsertFront(facePS.first, facePS.second);
    }

    while (edgeBuffer.Receive(q).size() < edgeBuffer.Receive(q).Size()) {
      std::pair<Point, ProcSet> edgePS;
      edgeBuffer.Receive(q) >> edgePS;
      meshAndOverlapProcSets.InsertFront(edgePS.first, edgePS.second);
    }
  }

  for (int i = 0; i < 2; i++) {
    ExchangeBuffer CB(commSplit);
    for (procset p = meshAndOverlapProcSets.Begin(); p != meshAndOverlapProcSets.End(); ++p) {
      for (int q = 0; q < p.size(); ++q) {
        CB.Send(p[q]) << p;
      }
    }
    CB.Communicate();
    for (int q = 0; q < PPM->Size(commSplit); ++q) {
      if (q == PPM->Proc(commSplit)) continue;
      while (CB.Receive(q).size() < CB.Receive(q).Size()) {
        std::pair<Point, ProcSet> receivedProcSetData;
        CB.Receive(q) >> receivedProcSetData;
        meshAndOverlapProcSets.Insert(receivedProcSetData.first, receivedProcSetData.second);
        Assert(receivedProcSetData.second.existselementof(q));
      }
    }
  }
}