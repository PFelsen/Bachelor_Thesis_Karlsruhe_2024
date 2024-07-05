#include "Matrix.hpp"
#include "Operator.hpp"
#include "Sparse.hpp"
#include "Vector.hpp"

Matrix::Matrix(const Vector &U, bool correctionAddBCArgyris) :
  VectorMatrixBase(U), data(U.Size()), Vec(U), sm(nullptr),
  collectAfterProd(true), addBC(correctionAddBCArgyris) {}

Matrix::Matrix(const Matrix &A) :
    VectorMatrixBase(A), data(A.data), Vec(A.Vec), sm(nullptr),
    collectAfterProd(A.collectAfterProd) {}

Matrix::Matrix(Matrix &&A) :
    VectorMatrixBase(A), data(std::move(A.data)), Vec(A.Vec), sm(nullptr),
    collectAfterProd(A.collectAfterProd) {}

Matrix::~Matrix() {
  if (sm) {
    delete sm;
    sm = nullptr;
  }
}

void Matrix::CreateSparse() {
  //  mout << "mem used: " << PPM->MemoryUsage() << ",  before Matrix::CreateSparse" << endl;
  if (sm) delete sm;
  sm = new SparseMatrix(*this);
  sm->compress();
  //  mout << "mem used: " << PPM->MemoryUsage() << ",  after Matrix::CreateSparse" << endl;
}

Matrix &Matrix::operator=(const Matrix &A) {
  data = A.data;
  //  for (int i = 0; i < Size(); ++i)
  //    (*this)[i] = A[i];
  return *this;
}

Matrix &Matrix::operator=(Matrix &&A) {
  data = std::move(A.data);
  //  for (int i = 0; i < Size(); ++i)
  //    (*this)[i] = A[i];
  return *this;
}

Matrix &Matrix::operator+=(const constAB<Scalar, Matrix> &aA) {
  Scalar s = aA.first();
  for (int i = 0; i < Size(); ++i)
    data[i] += s * aA.second().GetData()[i];
  return *this;
}

Matrix &Matrix::operator*=(const Scalar &a) {
  data *= a;
  return *this;
}

void Matrix::multiply_plus(Vector &b, const Vector &u) const {
  plusMatVec(b, u);
  if (collectAfterProd) {
    b.Collect();
  } else {
    b.SetAccumulateFlag(false);
  }
}

void Matrix::multiply_minus(Vector &b, const Vector &u) const {
  minusMatVec(b, u);
  if (collectAfterProd) {
    b.Collect();
  } else {
    b.SetAccumulateFlag(false);
  }
}

void Matrix::multiply(Vector &b, const Vector &u) const {
  matVec(b, u);
  if (collectAfterProd) {
    b.Collect();
  } else {
    b.SetAccumulateFlag(false);
  }
}

void Matrix::multiply_transpose_plus(Vector &b, const Vector &u) const {
  plusMatTVec(b, u);
  if (collectAfterProd) {
    b.Collect();
  } else {
    b.SetAccumulateFlag(false);
  }
}

void Matrix::multiply_plus(Vectors &b, const Vectors &u) const {
  plusMatVecs(b, u);
  if (collectAfterProd) {
    b.Collect();
  } else {
    b.SetAccumulateFlags(false);
  }
}

void Matrix::multiply(Vectors &b, const Vectors &u) const {
  matVecs(b, u);
  if (collectAfterProd) {
    b.Collect();
  } else {
    b.SetAccumulateFlags(false);
  }
}

void Matrix::multiply_transpose_plus(Vectors &b, const Vectors &u) const {
  plusMatTVecs(b, u);
  if (collectAfterProd) {
    b.Collect();
  } else {
    b.SetAccumulateFlags(false);
  }
}

void Matrix::minusMatVec(Vector &b, const Vector &u) const {
  if (!u.GetAccumulateFlag()) { Warning("Vector not accumulated in Matrix Vector product"); }
  applications++;
  if (sm) {
    sm->minusMatVec(b, u);
    return;
  }
  const Scalar *a = data();
  for (int i = 0; i < nR(); ++i) {
    int d = Diag(i);
    int n = Dof(i);
    for (int k = 0; k < n; ++k)
      for (int l = 0; l < n; ++l, ++a)
        b(i, k) -= *a * u(i, l);
    for (++d; d < Diag(i + 1); ++d) {
      int j = Column(d);
      int m = Dof(j);
      for (int k = 0; k < n; ++k)
        for (int l = 0; l < m; ++l, ++a)
          b(i, k) -= *a * u(j, l);
      if (SingleEntry(d)) continue;
      for (int k = 0; k < m; ++k)
        for (int l = 0; l < n; ++l, ++a)
          b(j, k) -= *a * u(i, l);
    }
  }
}

void Matrix::plusMatVec(Vector &b, const Vector &u) const {
  if (!u.GetAccumulateFlag()) { Warning("Vector not accumulated in Matrix Vector product"); }
  applications++;
  if (sm) {
    //        mout << "use sm " << endl;
    sm->plusMatVec(b, u);
    return;
  }
  const Scalar *a = data();
  for (int i = 0; i < nR(); ++i) {
    int d = Diag(i);
    int n = Dof(i);
    for (int k = 0; k < n; ++k)
      for (int l = 0; l < n; ++l, ++a)
        b(i, k) += *a * u(i, l);
    for (++d; d < Diag(i + 1); ++d) {
      int j = Column(d);
      int m = Dof(j);
      for (int k = 0; k < n; ++k)
        for (int l = 0; l < m; ++l, ++a)
          b(i, k) += *a * u(j, l);
      if (SingleEntry(d)) continue;
      for (int k = 0; k < m; ++k)
        for (int l = 0; l < n; ++l, ++a)
          b(j, k) += *a * u(i, l);
    }
  }
}

void Matrix::matVec(Vector &u, const Vector &v) const {
  u = 0;
  plusMatVec(u, v);
}

void Matrix::plusMatTVec(Vector &b, const Vector &u) const {
  if (!u.GetAccumulateFlag()) { Warning("Vector not accumulated in Matrix Vector product!"); }
  const Scalar *a = data();
  for (int i = 0; i < nR(); ++i) {
    int d = Diag(i);
    int n = Dof(i);
    for (int k = 0; k < n; ++k)
      for (int l = 0; l < n; ++l, ++a)
        b(i, l) += *a * u(i, k);
    for (++d; d < Diag(i + 1); ++d) {
      int j = Column(d);
      int m = Dof(j);
      for (int k = 0; k < n; ++k)
        for (int l = 0; l < m; ++l, ++a)
          b(j, l) += *a * u(i, k);
      if (SingleEntry(d)) continue;
      for (int k = 0; k < m; ++k)
        for (int l = 0; l < n; ++l, ++a)
          b(i, l) += *a * u(j, k);
    }
  }
}

void Matrix::plusMatVecs(Vectors &b, const Vectors &u) const {
  for (size_t i = 0; i < u.size(); ++i) {
    if (!u[i].GetAccumulateFlag()) {
      Warning("At least one Vector in Vectors not accumulated!");
      break;
    }
  }
  if (sm) {
    for (size_t i = 0; i < b.size(); ++i) {
      sm->plusMatVec(b[i], u[i]);
    }
    return;
  }
  const Scalar *a = data();
  for (int i = 0; i < nR(); ++i) {
    int d = Diag(i);
    int n = Dof(i);
    for (int k = 0; k < n; ++k)
      for (int l = 0; l < n; ++l, ++a)
        for (int v = 0; v < b.size(); ++v)
          b(v, i, k) += *a * u(v, i, l);
    for (++d; d < Diag(i + 1); ++d) {
      int j = Column(d);
      int m = Dof(j);
      for (int k = 0; k < n; ++k)
        for (int l = 0; l < m; ++l, ++a)
          for (int v = 0; v < b.size(); ++v)
            b(v, i, k) += *a * u(v, j, l);
      if (SingleEntry(d)) continue;
      for (int k = 0; k < m; ++k)
        for (int l = 0; l < n; ++l, ++a)
          for (int v = 0; v < b.size(); ++v)
            b(v, j, k) += *a * u(v, i, l);
    }
  }
}

void Matrix::matVecs(Vectors &u, const Vectors &v) const {
  u = 0;
  plusMatVecs(u, v);
}

void Matrix::plusMatTVecs(Vectors &b, const Vectors &u) const {
  for (size_t i = 0; i < u.size(); ++i) {
    if (!u[i].GetAccumulateFlag()) {
      Warning("At least one Vector in Vectors not accumulated!");
      break;
    }
  }
  const Scalar *a = data();
  for (int i = 0; i < nR(); ++i) {
    int d = Diag(i);
    int n = Dof(i);
    for (int k = 0; k < n; ++k)
      for (int l = 0; l < n; ++l, ++a)
        for (int v = 0; v < b.size(); ++v)
          b(v, i, l) += *a * u(v, i, k);
    for (++d; d < Diag(i + 1); ++d) {
      int j = Column(d);
      int m = Dof(j);
      for (int k = 0; k < n; ++k)
        for (int l = 0; l < m; ++l, ++a)
          for (int v = 0; v < b.size(); ++v)
            b(v, j, l) += *a * u(v, i, k);
      if (SingleEntry(d)) continue;
      for (int k = 0; k < m; ++k)
        for (int l = 0; l < n; ++l, ++a)
          for (int v = 0; v < b.size(); ++v)
            b(v, i, l) += *a * u(v, j, k);
    }
  }
}

void Matrix::copy(Scalar *a, int *ind, int *col) const {
  if (nR() == 0) return;
  const Scalar *A = data();
  for (int i = 0; i < Index(nR()); ++i)
    ind[i] = 0;
  for (int i = 0; i < nR(); ++i)
    for (int k = 0; k < Dof(i); ++k) {
      ind[Index(i) + k] += Dof(i);
      for (int d = Diag(i) + 1; d < Diag(i + 1); ++d) {
        int j = Column(d);
        int nc = Dof(j);
        for (int l = 0; l < nc; ++l) {
          ++ind[Index(i) + k];
          ++ind[Index(j) + l];
        }
      }
    }
  for (int i = 1; i < Index(nR()); ++i)
    ind[i] += ind[i - 1];
  for (int i = Index(nR()); i > 0; --i)
    ind[i] = ind[i - 1];
  ind[0] = 0;
  for (int i = 0; i < nR(); ++i) {
    int nc = Dof(i);
    for (int k = 0; k < nc; ++k) {
      int d = Diag(i);
      a[ind[Index(i) + k]] = A[Entry(d) + k * nc + k];
      col[ind[Index(i) + k]++] = Index(i) + k;
      for (int l = 0; l < nc; ++l) {
        if (l == k) continue;
        a[ind[Index(i) + k]] = A[Entry(d) + k * nc + l];
        col[ind[Index(i) + k]++] = Index(i) + l;
      }
    }
  }
  for (int i = 0; i < nR(); ++i) {
    int nr = Dof(i);
    for (int k = 0; k < nr; ++k) {
      int d = Diag(i);
      for (++d; d < Diag(i + 1); ++d) {
        int j = Column(d);
        int nc = Dof(j);
        for (int l = 0; l < nc; ++l) {
          a[ind[Index(i) + k]] = A[Entry(d) + k * nc + l];
          col[ind[Index(i) + k]++] = Index(j) + l;
          if (SingleEntry(d)) continue;
          a[ind[Index(j) + l]] = A[Entry(d) + nc * nr + l * nr + k];
          col[ind[Index(j) + l]++] = Index(i) + k;
        }
      }
    }
  }
  for (int i = Index(nR()); i > 0; --i)
    ind[i] = ind[i - 1];
  ind[0] = 0;
}

void Matrix::ClearDirichletValues() {
  const bool *Dirichlet = Vec.D();
  Scalar *A = data();
  int ind = 0;
  for (int i = 0; i < nR(); ++i) {
    int nr = Dof(i);
    for (int k = 0; k < nr; ++k) {
      int d = Diag(i);
      if (Dirichlet[ind++]) {
        int e = Entry(d);
        A[e + k * nr + k] = 1.0;
        for (int l = 0; l < nr; ++l)
          if (l != k) A[e + k * nr + l] = 0.0;
        for (++d; d < Diag(i + 1); ++d) {
          int j = Column(d);
          int nc = Dof(j);
          e = Entry(d);
          for (int l = 0; l < nc; ++l)
            A[e + k * nc + l] = 0.0;
          e += nr * nc;
          for (int l = 0; l < nc; ++l)
            if (Dirichlet[Index(j) + l]) A[e + l * nr + k] = 0.0;
        }
      } else {
        for (++d; d < Diag(i + 1); ++d) {
          int j = Column(d);
          int nc = Dof(j);
          int e = Entry(d) + nr * nc;
          for (int l = 0; l < nc; ++l)
            if (Dirichlet[Index(j) + l]) A[e + l * nr + k] = 0.0;
        }
      }
    }
  }
  // linear combinations of DoFs at bnd (for Argyris element)
  if (!Vec.BC()) return;
  if (!addBC) return;
  const BCEquations &bc = Vec.GetBCEquations();
  ind = 0;
  for (int i = 0; i < nR(); ++i) {
    int nr = Dof(i);
    for (int k = 0; k < nr; ++k, ind++) {
      int d = Diag(i);
      if (bc.Exists(ind)) {
        int e = Entry(d);
        A[e + k * nr + k] = 1.0;
        for (int l = 0; l < nr; ++l)
          if (l != k) {
            if (bc.Exists(ind, l)) {
              A[e + k * nr + l] = bc.Lambda(ind, l);
            } else {
              A[e + k * nr + l] = 0.0;
            }
          }
        for (++d; d < Diag(i + 1); ++d) {
          int j = Column(d);
          int nc = Dof(j);
          e = Entry(d);
          for (int l = 0; l < nc; ++l)
            A[e + k * nc + l] = 0.0;
          e += nr * nc;
          for (int l = 0; l < nc; ++l)
            if (bc.Exists(Index(j) + l)) A[e + l * nr + k] = 0.0;
        }
      } else {
        for (++d; d < Diag(i + 1); ++d) {
          int j = Column(d);
          int nc = Dof(j);
          int e = Entry(d) + nr * nc;
          for (int l = 0; l < nc; ++l)
            if (bc.Exists(Index(j) + l)) A[e + l * nr + k] = 0.0;
        }
      }
    }
  }
}

void Matrix::EliminateDirichlet() {
  const bool *Dirichlet = Vec.D();
  Scalar *A = data();
  if (identifysets() != identifysets_end())
    for (row r = rows(); r != rows_end(); ++r) {
      int n = r.n();
      for (entry e = r.entries(); e != r.entries_end(); ++e) {
        identifyset is = find_identifyset(e());
        if (is == identifysets_end()) continue;
        if (is.master()) continue;
        row r0 = find_row(is[0]);
        row r2 = find_row(e());
        Scalar *A00 = (*this)(r0, r0);
        Scalar *A22 = (*this)(r2, r2);
        Scalar *A01 = (*this)(r0, r);
        Scalar *A21 = (*this)(r2, r);
        int m = r2.n();
        for (int k = 0; k < n; ++k)
          for (int l = 0; l < m; ++l) {
            A01[k * m + l] += A21[k * m + l];
            A21[k * m + l] = 0;
            A00[k * m + l] += A22[k * m + l];
            A22[k * m + l] = 0;
          }
      }
    }
  for (identifyset is = identifysets(); is != identifysets_end(); ++is) {
    if (is.master()) continue;
    row r2 = find_row(is());
    row r0 = find_row(is[0]);
    Scalar *A00 = (*this)(r0, r0);
    Scalar *A22 = (*this)(r2, r2);
    Scalar *A20 = (*this)(r2, r0);
    int m = r2.n();
    for (entry e = r2.entries(); e != r2.entries_end(); ++e) {
      row r = find_row(e());
      int n = r.n();
      Scalar *A01 = (*this)(r0, r);
      Scalar *A21 = (*this)(r2, r);
      for (int k = 0; k < n; ++k)
        for (int l = 0; l < m; ++l) {
          A01[k * m + l] += A21[k * m + l];
          A21[k * m + l] = 0;
          A00[k * m + l] += A22[k * m + l];
          A22[k * m + l] = 0;
        }
      continue;

      int e00 = GetEntryX(r0, r0);
      int e22 = GetEntryX(r2, r2);
      int e10 = GetEntryX(r, r0);
      int e12 = GetEntryX(r, r2);

      mout << "Algebra " << OUT(e10) << OUT(e22) << OUT(e00) << OUT(e12) << endl;
    }
    for (int k = 0; k < m; ++k) {
      A22[k + m * k] = 1;
      A20[k + m * k] = -1;
    }
  }
  int ind = 0;
  for (int i = 0; i < nR(); ++i) {
    int nr = Dof(i);
    for (int k = 0; k < nr; ++k) {
      int d = Diag(i);
      if (Dirichlet[ind++]) {
        int e = Entry(d);
        A[e + k * nr + k] = 1.0;
        for (int l = 0; l < nr; ++l)
          if (l != k) {
            A[e + k * nr + l] = 0.0;
            A[e + l * nr + k] = 0.0;
          }
        for (++d; d < Diag(i + 1); ++d) {
          int j = Column(d);
          int nc = Dof(j);
          e = Entry(d);
          for (int l = 0; l < nc; ++l) {
            A[e + k * nc + l] = 0.0;
          }
          e += nr * nc;
          for (int l = 0; l < nc; ++l)
            if (Dirichlet[Index(j) + l]) {
              A[e + l * nr + k] = 0.0;
            } else A[e + l * nr + k] = 0.0;
        }
      } else {
        for (++d; d < Diag(i + 1); ++d) {
          int j = Column(d);
          int nc = Dof(j);
          int e = Entry(d);
          for (int l = 0; l < nc; ++l) {
            if (Dirichlet[Index(j) + l]) A[e + k * nc + l] = 0.0;
          }
          e += nr * nc;
          for (int l = 0; l < nc; ++l)
            if (Dirichlet[Index(j) + l]) A[e + l * nr + k] = 0.0;
        }
      }
    }
  }
}

void Matrix::Accumulate() {
  if (identify()) AccumulateIdentify();
  if (parallel()) AccumulateParallel();
}

void Matrix::CommunicateMatrix(ExchangeBuffer &exBuffer) {
  exBuffer.Communicate();
  Scalar *a = data();
  for (short q = 0; q < PPM->Size(CommSplit()); ++q)
    while (exBuffer.Receive(q).size() < exBuffer.ReceiveSize(q)) {
      Point x;
      exBuffer.Receive(q) >> x;
      row r = find_row(x);
      if (r == rows_end()) THROW("no row in Communicate Matrix")
      int id = r.Id();
      int d = Entry(Diag(id));
      int n = Dof(id);
      for (int j = 0; j < n * n; ++j) {
        Scalar b;
        exBuffer.Receive(q) >> b;
        a[d + j] += b;
      }
      int s;
      exBuffer.Receive(q) >> s;
      for (int k = 0; k < s; ++k) {
        Point y;
        int m;
        exBuffer.Receive(q) >> y >> m;
        int dd = r.GetEntryX(y);
        if (dd == -1) {
          Scalar tmp;
          for (int j = 0; j < m; ++j)
            exBuffer.Receive(q) >> tmp;
        } else {
          Scalar tmp;
          for (int j = 0; j < m; ++j) {
            exBuffer.Receive(q) >> tmp;
            a[dd + j] += tmp;
          }
        }
      }
    }
  exBuffer.ClearBuffers();
}

void Matrix::AccumulateIdentify() {
  ExchangeBuffer &exBuffer = Buffers().AccumulateMatrixIdentifyBuffer();
  Scalar *a = data();
  int q = PPM->Proc(CommSplit());
  for (identifyset is = identifysets(); is != identifysets_end(); ++is) {
    const Point &z = is();
    row r = find_row(z);
    int id = r.Id();
    int d = Entry(Diag(id));
    int n = Dof(id);
    for (int i = 0; i < is.size(); ++i) {
      Point y = is[i];
      exBuffer.Send(q) << y;
      for (int j = 0; j < n * n; ++j)
        exBuffer.Send(q) << a[d + j];
      exBuffer.Send(q) << int(r.size());
      for (entry e = r.entries(); e != r.entries_end(); ++e) {
        exBuffer.Send(q) << e();
        int m = 2 * n * Dof(e.Id());
        exBuffer.Send(q) << m;
        int dd = e.GetEntry();
        for (int j = 0; j < m; ++j)
          exBuffer.Send(q) << a[dd + j];
      }
    }
  }
  CommunicateMatrix(exBuffer);
}

void Matrix::AccumulateParallel() {
  ExchangeBuffer &exBuffer = Buffers().AccumulateMatrixBuffer();
  Scalar *a = data();
  for (procset p = procsets(); p != procsets_end(); ++p) {
    const Point &x = p();
    row r = find_row(x);
    int id = r.Id();
    int d = Entry(Diag(id));
    int n = Dof(id);
    for (int i = 0; i < p.size(); ++i) {
      int q = p[i];
      if (q == PPM->Proc(CommSplit())) continue;
      exBuffer.Send(q) << x;
      for (int j = 0; j < n * n; ++j)
        exBuffer.Send(q) << a[d + j];
      exBuffer.Send(q) << int(r.size());
      for (entry e = r.entries(); e != r.entries_end(); ++e) {
        exBuffer.Send(q) << e();
        int m = 2 * n * Dof(e.Id());
        exBuffer.Send(q) << m;
        int dd = e.GetEntry();
        for (int j = 0; j < m; ++j)
          exBuffer.Send(q) << a[dd + j];
      }
    }
  }
  CommunicateMatrix(exBuffer);
}

void Matrix::Transpose() {
  Scalar *a = data();
  for (int i = 0; i < nR(); ++i) {
    int d = Diag(i);
    int n = Dof(i);
    int e = Entry(d);

    vector<vector<Scalar>> mem_diag(n);

    for (int k = 0; k < n; ++k) {
      mem_diag[k].resize(n);
      for (int l = 0; l < n; ++l) {
        mem_diag[k][l] = a[e + k * n + l];
        a[e + k * n + l] = 0.0;
      }
    }
    for (int k = 0; k < n; ++k) {
      for (int l = 0; l < n; ++l) {
        a[e + k * n + l] += mem_diag[l][k];
      }
    }
    for (d = Diag(i) + 1; d < Diag(i + 1); ++d) {
      int e_col = Entry(d);
      int j = Column(d);
      int m = Dof(j);
      vector<vector<Scalar>> mem_0(n);
      vector<vector<Scalar>> mem_1(m);
      for (int k = 0; k < n; ++k) {
        mem_0[k].resize(m);
        for (int l = 0; l < m; ++l) {
          mem_0[k][l] = a[e_col + k * m + l];
          a[e_col + k * m + l] = 0.0;
        }
      }
      for (int l = 0; l < m; ++l) {
        mem_1[l].resize(n);
        for (int k = 0; k < n; ++k) {
          mem_1[l][k] = a[e_col + n * m + l * n + k];
          a[e_col + n * m + l * n + k] = 0.0;
        }
      }
      int cnt = 0;
      for (int l = 0; l < n; ++l) {
        for (int k = 0; k < m; ++k) {
          a[e_col + cnt] += mem_1[k][l];
          ++cnt;
        }
      }
      for (int l = 0; l < m; ++l) {
        for (int k = 0; k < n; ++k) {
          a[e_col + cnt] += mem_0[k][l];
          ++cnt;
        }
      }
    }
  }
}

void Matrix::Symmetric() {
  Scalar *A = data();
  for (int i = 0; i < nR(); ++i) {
    int n = Dof(i);
    int k = Diag(i);
    int e = Entry(k);
    for (int i0 = 0; i0 < n; ++i0)
      for (int i1 = 0; i1 < n; ++i1)
        A[e + i1 * n + i0] = A[e + i0 * n + i1];
    for (int k = Diag(i) + 1; k < Diag(i + 1); ++k) {
      int e = Entry(k);
      int j = Column(k);
      int m = Dof(j);
      for (int i0 = 0; i0 < n; ++i0) {
        for (int i1 = 0; i1 < m; ++i1) {
          A[e + n * m + i1 * m + i0] = A[e + i0 * n + i1];
        }
      }
    }
  }
}

bool Matrix::IsSymmetric() {
  const Scalar *A = data();
  for (int i = 0; i < nR(); ++i) {
    int n = Dof(i);
    int k = Diag(i);
    int e = Entry(k);
    for (int i0 = 0; i0 < n; ++i0)
      for (int i1 = 0; i1 < n; ++i1)
        if (!mpp_ba::isNear(A[e + i1 * n + i0], A[e + i0 * n + i1])) return false;
    for (int k = Diag(i) + 1; k < Diag(i + 1); ++k) {
      int e = Entry(k);
      int j = Column(k);
      int m = Dof(j);
      for (int i0 = 0; i0 < n; ++i0) {
        for (int i1 = 0; i1 < m; ++i1)
          if (!mpp_ba::isNear(A[e + n * m + i1 * m + i0], A[e + i0 * n + i1])) return false;
      }
    }
  }
  return true;
}

void Matrix::print() const {


  const Scalar *a = data();
  int zz = 0;
  for (int i = 0; i < nR(); ++i) {
    int d = Diag(i);
    int n = Dof(i);
    mout << n * n << "      ";
    for (int k = 0; k < n; ++k) {
      for (int l = 0; l < n; ++l, ++a) {
        mout << " " << *a;
        zz++;
      }
    }
    mout << endl;
    for (++d; d < Diag(i + 1); ++d) {
      int j = Column(d);
      int m = Dof(j);
      mout << n * m << "      ";
      for (int k = 0; k < n; ++k) {
        for (int l = 0; l < m; ++l, ++a) {
          mout << " " << *a;
          zz++;
        }
      }
      mout << endl;
      if (SingleEntry(d)) continue;
      mout << m * n << "      ";
      for (int k = 0; k < m; ++k) {
        for (int l = 0; l < n; ++l, ++a) {
          mout << " " << *a;
          zz++;
        }
      }
      mout << endl;
    }
  }
  //    for (row r = rows(); r != rows_end(); ++r) {
  //      mout << r() << " : " << (r.NumberOfEntries() + 1) << endl;
  //      mout << "      ";
  //      const Scalar *AA = this->operator()(r, r);
  //      for (int i = 0; i < r.n() * r.n(); ++i)
  //        mout << " " << AA[i];
  //      mout << endl;
  //      for (auto e = r.entries(); e != r.entries_end(); e++) {
  //        mout << "      ";
  //        row rr = find_row(e());
  //        const Scalar *AA = this->operator()(r, rr);
  //        for (int i = 0; i < r.n() * rr.n(); ++i) {
  //          mout << " " << AA[i];
  //        }
  //        mout << endl;
  //      }
  //      mout << endl;
  //      continue;
  //
  //
  //      mout << r() << " :";
  //      procset p = find_procset(r());
  //      if (p != procsets_end()) mout << " p " << p << endl;
  //      int id = r.Id();
  //      int d = Entry(Diag(id));
  //      int n = Dof(id);
  //      for (int i = 0; i < n; ++i)
  //        mout << " " << A[d + i];
  //      mout << endl;
  //    }
}

void Matrix::setCollectAfterProd(bool collect) { Matrix::collectAfterProd = collect; }

std::ostream &operator<<(std::ostream &s, const Matrix &A) {
  SparseMatrix A_sparse(A);
  return s << A_sparse;
}

constAB<Operator, Vector> operator*(const Matrix &A, const Vector &v) { return {A, v}; }

constAB<Operator, Vectors> operator*(const Matrix &A, const Vectors &v) { return {A, v}; }

constAB<Scalar, Matrix> operator*(const Scalar &a, const Matrix &A) { return {a, A}; }
