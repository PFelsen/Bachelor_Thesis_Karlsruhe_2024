#include "GaussSeidel.hpp"

#include "Sparse.hpp"
#include "Vector.hpp"

// TODO: Use BasicAlgbra
void applySmallMatrix2(int n, Scalar *u, const Scalar *a, const Scalar *b) {
  for (int i = 0; i < n; ++i) {
    Scalar s = 0;
    for (int k = 0; k < n; ++k)
      s += a[i * n + k] * b[k];
    u[i] = s;
  }
}

GaussSeidel::GaussSeidel() : Sp(0), A(0), shift_special(false) {
  int spec = 0;
  Config::Get("SHIFTspecial", spec);
  if (spec == 1) shift_special = true;
}

void GaussSeidel::Destruct() {
  if (Sp) delete Sp;
  Sp = 0;
  if (A) delete A;
  A = 0;
}

void GaussSeidel::Construct(const Matrix &_A) {
  Destruct();
  A = new Matrix(_A);
  A->Accumulate();
  Sp = new SparseMatrix(*A);
  Sp->CheckDiagonal();
}

void GaussSeidel::multiply(Vector &u, const Vector &b) const {
  u.Clear();
  Vector r(b);
  Sp->GaussSeidel(u, r, shift_special);
  u.Accumulate();
}

void GaussSeidel::multiply_transpose(Vector &u, const Vector &b) const {
  Vector r(b);
  r.Accumulate();
  Sp->BackwardGaussSeidel(u, r, shift_special);
  u.MakeAdditive();
  u.Accumulate();
}

std::ostream &operator<<(std::ostream &s, const GaussSeidel &GS) { return s << *(GS.Sp); }

void BlockGaussSeidel::Construct(const Matrix &_A) {
  A = &_A;
  BlockJacobi::Construct_Matrix_free(_A);
}

void BlockGaussSeidel::Destruct() {
  if (transposed_nj) delete[] transposed_nj;
  transposed_nj = 0;
  if (transposed_j) delete[] transposed_j;
  transposed_j = 0;
  if (transposed_entries) delete[] transposed_entries;
  transposed_entries = 0;
  BlockJacobi::Destruct();
}

void BlockGaussSeidel::multiply(Vector &u, const Vector &b) const {
  Vector r(b);
  const Scalar *diag = dd;
  const Scalar *a = A->GetData()();

  for (int i = 0; i < u.nR(); ++i, ++i) {
    {
      int d = u.Diag(i);
      int n1 = u.Dof(i);
      int n2 = u.Dof(i + 1);
      int n = n1 + n2;

      if (dd_entry.find(i) == dd_entry.end()) THROW("Entry not found")
      const int dd_e = dd_entry.find(i)->second;

      for (int k = 0; k < n; ++k) {

        Scalar R = b(i, k);
        if (k < n1)
          for (int d = u.Diag(i) + 1; d < u.Diag(i + 1); ++d) {
            int e_col = u.Entry(d);
            int j = u.Column(d);
            int m = u.Dof(j);
            for (int l = 0; l < m; ++l)
              R -= a[e_col + k * m + l] * u(j, l);
          }
        else
          for (int d = u.Diag(i + 1) + 1; d < u.Diag(i + 2); ++d) {
            int e_col = u.Entry(d);
            int j = u.Column(d);
            if (j == i) continue;
            int m = u.Dof(j);
            for (int l = 0; l < m; ++l)
              R -= a[e_col + (k - n1) * m + l] * u(j, l);
          }
        r(i, k) = R;
      }
      applySmallMatrix2(n, u(i), diag + dd_e, r(i));
    }
  }

  u.MakeAdditive();
  u.Accumulate();
}

PointBlockGaussSeidel::PointBlockGaussSeidel(bool forward) { Config::Get("PGSForward", forward); }

void PointBlockGaussSeidel::Construct(const Matrix &_A) {
  A = &_A;
  PointBlockJacobi::Construct(_A);

  return;
  // for debugging
  const Vector &u = A->GetVector();
  Point x[u.nR()];
  for (row r = u.rows(); r != u.rows_end(); ++r) {
    x[r.Id()] = r();
    mout << OUT(r()) << OUT(r.Id()) << OUT(u.master(r())) << endl;
  }
  for (int i = 0; i < u.nR(); ++i) {
    int n = u.Dof(i);
    mout << OUT(i) << "\t" << x[i] << endl;
    for (int d = u.Diag(i) + 1; d < u.Diag(i + 1); ++d) {
      int j = u.Column(d);
      mout << "  " << OUT(j) << "\t" << x[j] << " " << OUT(A->SingleEntry(d)) << endl;
    }
  }
}

void PointBlockGaussSeidel::Destruct() { PointBlockJacobi::Destruct(); }

void PointBlockGaussSeidel::multiply(Vector &u, const Vector &b) const {
  Vector r(b);
  const Scalar *a = A->GetData()();
  u.Clear();
  if (forward)
    for (int i = 0; i < u.nR(); ++i) {
      if (dd_entry[i] == -1) continue;
      int n = u.Dof(i);
      for (int d = u.Diag(i) + 1; d < u.Diag(i + 1); ++d) {
        const Scalar *aa = a + u.Entry(d);
        int j = u.Column(d);
        int m = u.Dof(j);
        for (int k = 0; k < n; ++k) {
          Scalar R = 0;
          for (int l = 0; l < m; ++l, ++aa)
            R += *aa * u(j, l);
          r(i, k) -= R;
        }
      }
      apply(n, u(i), dd + dd_entry[i], r(i));
    }
  else
    for (int i = u.nR() - 1; i >= 0; --i) {
      if (dd_entry[i] == -1) continue;
      int n = u.Dof(i);
      apply(n, u(i), dd + dd_entry[i], r(i));
      for (int d = u.Diag(i) + 1; d < u.Diag(i + 1); ++d) {
        if (u.SingleEntry(d)) continue;
        const Scalar *aa = a + u.Entry(d) + n * n;
        int j = u.Column(d);
        int m = u.Dof(j);
        for (int k = 0; k < m; ++k) {
          Scalar R = 0;
          for (int l = 0; l < n; ++l, ++aa)
            R += *aa * u(i, l);
          r(j, k) -= R;
        }
      }
    }
  u.Accumulate();
}

void PointBlockGaussSeidel::multiply_transpose(Vector &u, const Vector &b) const {
  Vector r(b);
  r.MakeAdditive();
  const Scalar *a = A->GetData()();
  u = 0;
  if (forward)
    for (int i = u.nR() - 1; i >= 0; --i) {
      if (dd_entry[i] == -1) continue;
      int n = u.Dof(i);
      apply(n, u(i), dd + dd_entry[i], r(i));
      for (int d = u.Diag(i) + 1; d < u.Diag(i + 1); ++d) {
        if (u.SingleEntry(d)) continue;
        const Scalar *aa = a + u.Entry(d) + n * n;
        int j = u.Column(d);
        int m = u.Dof(j);
        for (int k = 0; k < m; ++k) {
          Scalar R = 0;
          for (int l = 0; l < n; ++l, ++aa)
            R += *aa * u(i, l);
          r(j, k) -= R;
        }
      }
    }
  else
    for (int i = 0; i < u.nR(); ++i) {
      if (dd_entry[i] == -1) continue;
      int n = u.Dof(i);
      for (int d = u.Diag(i) + 1; d < u.Diag(i + 1); ++d) {
        const Scalar *aa = a + u.Entry(d);
        int j = u.Column(d);
        int m = u.Dof(j);
        for (int k = 0; k < n; ++k) {
          Scalar R = 0;
          for (int l = 0; l < m; ++l, ++aa)
            R += *aa * u(j, l);
          r(i, k) -= R;
        }
      }
      apply(n, u(i), dd + dd_entry[i], r(i));
    }
  u.MakeAdditive();
  u.Accumulate();
}