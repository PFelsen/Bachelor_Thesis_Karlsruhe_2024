#include "Jacobi.hpp"

#include "Lapdef.hpp"
#include "Matrix.hpp"
#include "Vector.hpp"
extern "C" void dgemv_(char *TRANS, int *M, int *N, const Scalar *alpha, const void *A,
                       const int *LDA, const void *X, int *INCX, const Scalar *BETA, void *Y,
                       int *INCY);

void invertSmallMatrixInPlaceLAPACK(int n, Scalar *a) {
  int IPIV[n];
  int info;
  int lwork = 2 * n;
  double work[2 * lwork];
  GETRF(&n, &n, a, &n, IPIV, &info);
  GETRI(&n, a, &n, IPIV, work, &lwork, &info);
}

void applySmallMatrixLAPACK(int n, Scalar *u, const Scalar *a, const Scalar *b,
                            const double &theta) {
  int INC = 1;
  char TRANS = 'T';
  double beta = 0;
  LGEMV(&TRANS, &n, &n, &theta, a, &n, b, &INC, &beta, u, &INC);
}

Jacobi::Jacobi() : D(0), theta(1), shift_special(false) {
  Config::Get("PreconditionerDamp", theta);
  int spec = 0;
  Config::Get("SHIFTspecial", spec);
  if (spec == 0) shift_special = false;
  else shift_special = true;
}

void Jacobi::Construct(const Matrix &A) {
  D = new Vector(A.GetSharedDisc(), A.Level());
  // D = new Vector(A.GetVector());
  const Scalar *a = A.GetData()();
  for (int i = 0; i < A.nR(); ++i) {
    int d = A.Diag(i);
    int n = A.Dof(i);
    for (int k = 0; k < n; ++k)
      for (int l = 0; l < n; ++l, ++a)
        if (k == l) (*D)(i, k) = *a;
    for (++d; d < A.Diag(i + 1); ++d) {
      int j = A.Column(d);
      int m = A.Dof(j);
      a += 2 * m * n;
    }
  }
  D->Accumulate();
  vout(10) << "diag " << *D;

  for (int i = 0; i < D->size(); ++i) {
    if ((*D)[i] == 0.0) (*D)[i] = 1.0;
    else (*D)[i] = 1.0 / (*D)[i];
  }
}

void Jacobi::Destruct() {
  if (D) delete D;
  D = 0;
}

void Jacobi::multiply(Vector &u, const Vector &b) const {
  for (int i = 0; i < D->size(); ++i)
    u[i] = (*D)[i] * b[i];
  u *= theta;
  u.SetAccumulateFlag(false);
  u.Accumulate();
}

void DampedJacobi::multiply(Vector &u, const Vector &b) const {
  for (int i = 0; i < D->size(); ++i)
    u[i] = damp * (*D)[i] * b[i];
  u.Accumulate();
}

BlockJacobi::BlockJacobi() : dd(0), blocksize(2) {}

void BlockJacobi::Construct(const Matrix &A) {
  Construct_Matrix_free(A);
  return;
}

void BlockJacobi::Construct_Matrix_free(const Matrix &A) {
  int dofs_on_diag = 0;
  const Scalar *a = A.GetData()();

  for (int i = 0; i < A.nR(); ++i, ++i)
    dofs_on_diag += pow(A.Dof(i) + A.Dof(i + 1), 2);

  dd = new Scalar[dofs_on_diag];
  int dd_entry_cnt = 0;

  for (int i = 0; i < A.nR(); ++i, ++i) {
    int e1 = A.Entry(A.Diag(i));
    int e2 = A.Entry(A.Diag(i + 1));
    int n1 = A.Dof(i);
    int n2 = A.Dof(i + 1);

    int e1e2 = -1;

    for (int d = A.Diag(i + 1) + 1; d < A.Diag(i + 2); ++d) {
      if (A.Column(d) == i) {
        e1e2 = A.Entry(d);
        break;
      }
    }

    int n = n1 + n2;
    for (int k = 0; k < n; ++k) {
      bool zero_only = true;
      for (int l = 0; l < n; ++l) {
        if (k < n1 && l < n1) dd[dd_entry_cnt + n * k + l] = a[e1 + n1 * k + l];
        if (k >= n1 && l >= n1) dd[dd_entry_cnt + n * k + l] = a[e2 + n2 * (k - n1) + (l - n1)];
        if (e1e2 != -1) {
          if (k >= n1 && l < n1) dd[dd_entry_cnt + n * k + l] = a[e1e2 + n2 * (k - n1) + l];
          if (k < n1 && l >= n1)
            dd[dd_entry_cnt + n * k + l] = a[e1e2 + n1 * n2 + n1 * k + (l - n1)];
        } else {
          if (k >= n1 && l < n1) dd[dd_entry_cnt + n * k + l] = 0.0;
          if (k < n1 && l >= n1) dd[dd_entry_cnt + n * k + l] = 0.0;
        }

        if (dd[dd_entry_cnt + n * k + l] != 0.0) zero_only = false;
      }
      if (zero_only) dd[dd_entry_cnt + n * k + k] = 1;
    }

    // 	    mout<<"invert "<<i<<" "<<dd_entry_cnt<<" "<<n<<endl;
    invertSmallMatrix(n, dd + dd_entry_cnt);
    dd_entry[i] = dd_entry_cnt;
    dd_entry_cnt += n * n;
  }
}

void BlockJacobi::Destruct() {
  if (dd) {
    delete[] dd;
    dd_entry.clear();
  }
  dd = 0;
}

void BlockJacobi::multiply(Vector &u, const Vector &b) const {
  double theta = 1.0;
  for (int i = 0; i < u.nR(); ++i, ++i) {
    int n = u.Dof(i) + u.Dof(i + 1);
    const int dd_e = dd_entry.find(i)->second;
    applySmallMatrix(n, u(i), dd + dd_e, b(i));
  }
  u.Accumulate();
}

void BlockJacobi::multiply_transpose(Vector &u, const Vector &b) const { multiply(u, b); }

PointBlockJacobi::PointBlockJacobi() : theta(1), dd(0) {
  Config::Get("PreconditionerDamp", theta);
  Config::Get("PointBlockJacobiLapack", LAPACKFlag);
}

void PointBlockJacobi::apply(int n, Scalar *u, const Scalar *a, const Scalar *b) const {
  if (LAPACKFlag) applySmallMatrixLAPACK(n, u, a, b, theta);
  else applySmallMatrix(n, u, a, b);
}

void PointBlockJacobi::Construct(const Matrix &A) {
  vout(2) << "Construct " << Name() << " " << A.pSize() << endl;

  int dofs_on_diag = 0;
  const Scalar *a = A.GetData()();

  bool *master = new bool[A.nR()];
  for (row r = A.rows(); r != A.rows_end(); ++r)
    master[r.Id()] = A.master(r());

  for (int i = 0; i < A.nR(); ++i)
    if (master[i]) dofs_on_diag += A.Dof(i) * A.Dof(i);

  dd = new Scalar[dofs_on_diag];
  dd_entry = new int[A.nR()];
  int dd_entry_cnt = 0;

  for (int i = 0; i < A.nR(); ++i) {
    if (!master[i]) {
      dd_entry[i] = -1;
      continue;
    }
    int e = A.Entry(A.Diag(i));
    int n = A.Dof(i);
    for (int k = 0; k < n; ++k) {
      bool zero_only = true;
      for (int l = 0; l < n; ++l) {
        dd[dd_entry_cnt + n * k + l] = a[e + n * k + l];
        if (dd[dd_entry_cnt + n * k + l] != 0.0) zero_only = false;
      }
      if (zero_only) dd[dd_entry_cnt + n * k + k] = 1;
    }
    if (LAPACKFlag) invertSmallMatrixInPlaceLAPACK(n, dd + dd_entry_cnt);
    else invertSmallMatrix(n, dd + dd_entry_cnt);
    dd_entry[i] = dd_entry_cnt;
    dd_entry_cnt += n * n;
  }
  delete[] master;
}

void PointBlockJacobi::Destruct() {
  if (dd) {
    delete[] dd;
    delete[] dd_entry;
  }
  dd = 0;
}

void PointBlockJacobi::multiply(Vector &u, const Vector &b) const {
  //  u = 0;
  u.Clear();
  for (int i = 0; i < u.nR(); ++i) {
    if (dd_entry[i] == -1) continue;
    int e = u.Entry(u.Diag(i));
    int n = u.Dof(i);
    apply(n, u(i), dd + dd_entry[i], b(i));
  }
  if (!LAPACKFlag) u *= theta;
  u.Accumulate();
}

void PointBlockJacobi::multiply_transpose(Vector &u, const Vector &b) const { multiply(u, b); }

PointBlockJacobi_dG::PointBlockJacobi_dG() : D(0), theta(1) {
  Config::Get("PreconditionerDamp", theta);
}

void PointBlockJacobi_dG::Construct(const Matrix &A) { // this constructor inverts diagonal blocks
                                                       // in place and stores it in the Matrix D...

  D = new Matrix(A);

  D->Accumulate(); // should not be needed in this DG method

  Scalar *d = (*D).GetData()();
  const Vector &v = A.GetVector();
  for (int i = 0; i < A.nR(); ++i) {
    int e = A.Entry(A.Diag(i));
    int n = A.Dof(i);
    for (int k = 0; k < n; ++k)
      if (d[e + n * k] == 0.0) d[e + n * k] = 1;
  }
  for (int i = 0; i < A.nR(); ++i) {
    int e = A.Entry(A.Diag(i));
    int n = A.Dof(i);
    for (int k = 0; k < n; ++k)
      if (d[e + n * k] == 0.0) d[e + n * k] = 1;
    invertSmallMatrix(n, d + e);
  }
}

void PointBlockJacobi_dG::Destruct() {
  if (D) delete D;
  D = 0;
}

void PointBlockJacobi_dG::multiply(Vector &u, const Vector &b) const {
  const Scalar *a = D->GetData()();
  for (int i = 0; i < u.nR(); ++i) {
    int e = u.Entry(u.Diag(i));
    int n = u.Dof(i);
    applySmallMatrix(n, u(i), a + e, b(i));
  }
  u *= theta;
  u.Accumulate();
}

CellWiseJacobi::CellWiseJacobi() : theta(1) { Config::Get("PreconditionerDamp", theta); }

RMatrix getRMatrixFromMatrixAndCell(const Matrix &A, const Cell &C) {
  std::vector<Point> z = A.GetDoF().GetStoragePoints(C);
  std::vector<row> R(z.size());
  int n = 0;
  for (int i = 0; i < R.size(); ++i) {
    R[i] = A.find_row(z[i]);
    n += R[i].n();
  }
  RMatrix a(0.0, n, n);
  std::vector<double> &data = a.Data();
  for (int i = 0, j = 0; i < R.size(); ++i)
    for (int k = 0; k < R[i].n(); ++k, ++j)
      for (int i1 = 0, j1 = 0; i1 < R.size(); ++i1)
        for (int k1 = 0; k1 < R[i1].n(); ++k1, ++j1)
          data[n * j + j1] = A(R[i], R[i1])[k * R[i1].n() + k1];

  return a;
}

void CellWiseJacobi::Construct(const Matrix &A) {
  InvMap.reserve(A.GetMesh().CellCount());
  for (cell c = A.cells(); c != A.cells_end(); ++c) {
    const Cell &C = *c;
    InvMap.emplace(C(), getRMatrixFromMatrixAndCell(A, C).Invert());
  }
}

void CellWiseJacobi::Destruct() {}

void CellWiseJacobi::multiply(Vector &u, const Vector &b) const {
  u.Clear();
  for (cell c = u.cells(); c != u.cells_end(); ++c) {
    const SparseRMatrix &SpSm = InvMap.at(c());
    int RowIdU = u.find_row(c()).Id();
    int RowIdB = u.find_row(c()).Id();
    for (int i = 0; i < SpSm.rows(); ++i) {
      for (int k = SpSm.rowind(i); k < SpSm.rowind(i + 1); ++k) {
        u(RowIdU, i) += SpSm.nzval(k) * b(RowIdB, SpSm.colptr(k));
      }
    }
  }
  u *= theta;
  u.Accumulate();
}