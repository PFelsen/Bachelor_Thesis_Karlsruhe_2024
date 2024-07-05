#include "MatrixTransfer.hpp"

#include "DGDiscretization.hpp"
#include "LagrangeDiscretization.hpp"
#include "MixedEGDiscretization.hpp"

MatrixTransfer::TransferMatrixGraph::TransferMatrixGraph(const VectorMatrixBase &_cg,
                                                         const VectorMatrixBase &_fg) :
    fg(_fg), cg(_cg), simple(false) {
  Config::Get("SimpleMatrixTransfer", simple);
  for (row r = cg.rows(); r != cg.rows_end(); ++r)
    Rows::Insert(r(), r.n());
  for (cell c = cg.cells(); c != cg.cells_end(); ++c)
    for (int i = 0; i < c.Children(); ++i)
      addCell(*c, *(fg.find_cell(c.Child(i))));
  numbering();
  int nR = int(cg.nR());
  int S = cg.size();
  if (simple) S = nR;
  index = new int[S + 1];
  diag = new int[nR + 1];
  column = new int[m];
  matentry = new int[m];
  vector<row> r(nR);
  for (row ri = rows(); ri != rows_end(); ++ri)
    r[ri.Id()] = ri;
  int d = 0;
  index[0] = 0;
  diag[0] = 0;
  matentry[d] = 0;
  for (int i = 0; i < nR; ++i) {
    if (simple) index[i + 1] = index[i] + 1;
    else index[i + 1] = index[i] + cg.Dof(i);
    for (entry e = r[i].entries(); e != r[i].entries_end(); ++e) {
      column[d] = e.Id();
      matentry[d] = e.GetEntry();
      ++d;
    }
    diag[i + 1] = d;
  }
}

void MatrixTransfer::TransferMatrixGraph::addCell(const Cell &cc, const Cell &fc) {
  //  vector<Point> cz = cg.GetNodalPoints(cc);
  //  vector<Point> fz = fg.GetNodalPoints(fc);
  vector<Point> cz = cg.GetDoF().GetStoragePoints(cc);
  vector<Point> fz = fg.GetDoF().GetStoragePoints(fc);
  for (int i = 0; i < cz.size(); ++i) {
    Rows::Insert(cz[i], 1);
    for (int j = 0; j < fz.size(); ++j)
      Rows::AddEntryFixedOrder(cz[i], fz[j]);
  }
}

void MatrixTransfer::TransferMatrixGraph::numbering() {
  int nR = int(cg.nR());
  vector<RowIterator> r(nR);
  for (auto ri = Rows::begin(); ri != Rows::end(); ++ri)
    r[cg.Id(ri->first)] = ri;
  d = 0;
  m = 0;
  for (int i = 0; i < nR; ++i) {
    r[i]->second.Id(i);
    for (auto e = (r[i]->second).Entries().begin(); e != (r[i]->second).Entries().end(); ++e) {
      int j = fg.Id(e->first);
      e->second.Id(j);
      e->second.EntryIndex(m);
      ++d;
      if (simple) ++m;
      else m += r[i]->second.NumberOfDofs() * fg.Dof(j);
    }
  }
}

Scalar *MatrixTransfer::TransferMatrix::operator()(const Point &x, const Point &y) {
  row r = TG.find_row(x);
  return (*this)() + r.GetEntry(y);
}

const Scalar *MatrixTransfer::TransferMatrix::operator()(const Point &x, const Point &y) const {
  row r = TG.find_row(x);
  return (*this)() + r.GetEntry(y);
}

void MatrixTransfer::TransferMatrix::multiply_plus(Vector &fv, const Vector &cv) const {
  const Scalar *a = (*this)();
  if (Simple()) {
    for (int i = 0; i < TG.nR(); ++i) {
      for (int d = TG.Diag(i); d < TG.Diag(i + 1); ++d, ++a) {
        int j = TG.Column(d);
        for (int k = 0; k < fv.Dof(j); ++k)
          fv(j, k) += *a * cv(i, k);
      }
    }
  } else {
    for (int i = 0; i < TG.nR(); ++i) {
      for (int d = TG.Diag(i); d < TG.Diag(i + 1); ++d) {
        int j = TG.Column(d);
        for (int k = 0; k < fv.Dof(j); ++k)
          for (int l = 0; l < cv.Dof(i); ++l, ++a)
            fv(j, k) += *a * cv(i, l);
      }
    }
  }
  fv.ClearDirichletValues();
  fv.MakeAdditive();
  fv.Accumulate();
}

void MatrixTransfer::TransferMatrix::multiply_transpose_plus(Vector &cv, const Vector &fv) const {
  const Scalar *a = (*this)();
  if (Simple()) {
    for (int i = 0; i < TG.nR(); ++i) {
      for (int d = TG.Diag(i); d < TG.Diag(i + 1); ++d, ++a) {
        int j = TG.Column(d);
        for (int k = 0; k < fv.Dof(j); ++k)
          cv(i, k) += *a * fv(j, k);
      }
    }
  } else {
    for (int i = 0; i < TG.nR(); ++i) {
      for (int d = TG.Diag(i); d < TG.Diag(i + 1); ++d) {
        int j = TG.Column(d);
        for (int k = 0; k < fv.Dof(j); ++k)
          for (int l = 0; l < cv.Dof(i); ++l, ++a)
            cv(i, l) += *a * fv(j, k);
      }
    }
  }
  cv.ClearDirichletValues();
  cv.Collect();
}

void MatrixTransfer::constructLagrangeTransferMatrix() {
  tMatrix = 0;
  const auto &cg = tMatrix.CoarseMatrixGraph();
  for (cell c = cg.cells(); c != cg.cells_end(); ++c) {
    for (int cf = 0; cf < c.Corners(); ++cf) {
      tMatrix(c.Corner(cf), c.Corner(cf))[0] = 1;
    }
    for (int ce = 0; ce < c.Edges(); ++ce) {
      for (int cf = 0; cf < 2; ++cf) {
        tMatrix(c.EdgeCorner(ce, cf), c.Edge(ce))[0] = 0.5;
      }
    }
    if (c.Type() == QUADRILATERAL || c.Type() == HEXAHEDRON) {
      for (int cf = 0; cf < c.Corners(); ++cf) {
        tMatrix(c.Corner(cf), c())[0] = 0.25;
      }
    }
  }
}

void MatrixTransfer::constructDGTransferMatrix() {
  tMatrix = 0;
  const auto &cg = tMatrix.CoarseMatrixGraph();
  const auto &fg = tMatrix.FineMatrixGraph();
  for (cell c = cg.cells(); c != cg.cells_end(); ++c) {
    const Shape &S = cg.GetDisc().GetShape(*c);
    vector<Point> z = cg.GetDoF().GetNodalPoints(*c);
    for (int k = 0; k < c.Children(); ++k) {
      cell fc = fg.find_cell(c.Child(k));
      vector<Point> z_fg = fg.GetDoF().GetNodalPoints(*fc);
      for (int i = 0; i < z_fg.size(); ++i) {
        Point lz = c.GlobalToLocal(z_fg[i]);
        for (int j = 0; j < z.size(); ++j) {
          tMatrix(c(), c.Child(k))[i * z.size() + j] = S(lz, j);
        }
      }
    }
  }
}

void MatrixTransfer::constructEGTransferMatrix() {
  tMatrix = 0;
  const auto &cg = tMatrix.CoarseMatrixGraph();
  for (cell c = cg.cells(); c != cg.cells_end(); ++c) {
    for (int cf = 0; cf < c.Corners(); ++cf)
      tMatrix(c.Corner(cf), c.Corner(cf))[0] = 1;
    for (int ce = 0; ce < c.Edges(); ++ce)
      for (int cf = 0; cf < 2; ++cf)
        tMatrix(c.EdgeCorner(ce, cf), c.Edge(ce))[0] = 0.5;
    if (c.Type() == QUADRILATERAL || c.Type() == HEXAHEDRON) {
      for (int cf = 0; cf < c.Corners(); ++cf)
        tMatrix(c.Corner(cf), c())[0] = 0.25;
    }
  }
}

MatrixTransfer::MatrixTransfer(const Vector &coarse, const Vector &fine) :
    tMatrixGraph(coarse, fine), tMatrix(tMatrixGraph) {
  if (typeid(coarse.GetDisc()) == typeid(LagrangeDiscretization)) constructLagrangeTransferMatrix();
  else if (typeid(coarse.GetDisc()) == typeid(DGDiscretization)) constructDGTransferMatrix();
  else if (typeid(coarse.GetDisc()) == typeid(MixedEGDiscretization)) constructEGTransferMatrix();
  else Exit("MatrixTransfer is not implemented")
}

void MatrixTransfer::Prolongate(const Vector &coarse, Vector &fine) const {
  fine = tMatrix * coarse;
  fine.ClearDirichletValues();
}

void MatrixTransfer::ProlongateTransposed(Vector &coarse, const Vector &fine) const {
  coarse = fine * tMatrix;
  coarse.ClearDirichletValues();
  coarse.Collect();
}

void MatrixTransfer::Restrict(Vector &coarse, const Vector &fine) const {
  ProlongateTransposed(coarse, fine);
}

void MatrixTransfer::Project(Vector &coarse, const Vector &fine) const {
  // TODO: Projection missing
  Warning("Projection is not implemented") coarse.Average();
}

RMatrix MatrixTransfer::AsMatrix() const { Warning("AsMatrix is not implemented") return {}; }
