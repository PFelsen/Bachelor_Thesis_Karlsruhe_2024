#include "FullMatrixGraph.hpp"

void FullMatrixGraph::FullCellsOverlap() {
  int commSplit = CommSplit();
  ExchangeBuffer exBuffer(commSplit);
  for (cell c = cells(); c != cells_end(); ++c) {
    CellBoundaryFaces bf(mesh, c);

    meshAndOverlapProcSets.Add(c(), PPM->Proc(CommSplit()));
    for (int i = 0; i < c.Corners(); ++i)
      meshAndOverlapProcSets.Add(c.Corner(i), PPM->Proc(CommSplit()));
    for (int i = 0; i < c.Edges(); ++i)
      meshAndOverlapProcSets.Add(c.Edge(i), PPM->Proc(CommSplit()));

    for (int q = 0; q < PPM->Size(commSplit); ++q) {
      if (q == PPM->Proc(commSplit)) continue;
      exBuffer.Send(q) << c << bf;
      meshAndOverlapProcSets.Append(c(), q);
      for (int i = 0; i < c.Corners(); ++i)
        meshAndOverlapProcSets.Append(c.Corner(i), q);
      for (int i = 0; i < c.Edges(); ++i)
        meshAndOverlapProcSets.Append(c.Edge(i), q);
    }
  }
  exBuffer.Communicate();
  for (int q = 0; q < PPM->Size(commSplit); ++q) {
    if (q == PPM->Proc(commSplit)) continue;
    while (exBuffer.Receive(q).size() < exBuffer.Receive(q).Size()) {
      int tp;
      short n, sd;
      exBuffer.Receive(q) >> tp >> sd >> n;
      vector<Point> x(n);
      for (int i = 0; i < n; ++i)
        exBuffer.Receive(q) >> x[i];
      Cell *C = CreateCell(CELLTYPE(tp), sd, x);
      InsertOverlapCell(C);
      meshAndOverlapProcSets.Add(C->Center(), q);
      for (int i = 0; i < C->Faces(); ++i)
        InsertOverlapFace(C->Face(i), C->Center());
      for (int i = 0; i < C->Corners(); ++i)
        meshAndOverlapProcSets.Add(C->Corner(i), q);
      for (int i = 0; i < C->Edges(); ++i)
        meshAndOverlapProcSets.Add(C->Edge(i), q);
      for (int i = 0; i < C->Faces(); ++i)
        meshAndOverlapProcSets.Add(C->Face(i), q);
      for (int qq = 0; qq < PPM->Size(commSplit); ++qq) {
        meshAndOverlapProcSets.Append(C->Center(), qq);
        for (int i = 0; i < C->Corners(); ++i)
          meshAndOverlapProcSets.Append(C->Corner(i), qq);
        for (int i = 0; i < C->Edges(); ++i)
          meshAndOverlapProcSets.Append(C->Edge(i), qq);
      }
      //      exBuffer.Receive(q) >> n;
      //      for (int i = 0; i < n; ++i)
      //        exBuffer.Receive(q) >> mesh.bnd_ref();
    }
  }
  exBuffer.ClearBuffers();
}
