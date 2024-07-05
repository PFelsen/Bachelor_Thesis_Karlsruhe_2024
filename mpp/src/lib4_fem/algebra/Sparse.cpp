#include "Sparse.hpp"
#include "Parallel.hpp"

#include "Matrix.hpp"
#include "Vector.hpp"

#include <algorithm>
#include <iomanip>
#include <map>
#include <vector>

int NonZeros(const SparseMatrix &M) {
  int nz = 0;
  for (int i = 0; i < M.Size(); ++i)
    if (isnotZero(M.nzval(i))) ++nz;
  return nz;
}

int SparseMatrix::size_(const SparseMatrix &M, const vector<bool> &mask) const {
  int n = 0;
  for (int i = 0; i < M.size(); ++i)
    if (mask[i]) ++n;
  return n;
}

int SparseMatrix::Size_(const SparseMatrix &M, const vector<bool> &mask) const {
  int n = 0;
  for (int i = 0; i < M.size(); ++i)
    if (mask[i])
      for (int k = M.rowind(i); k < M.rowind(i + 1); ++k)
        if (mask[M.colptr(k)]) ++n;
  return n;
}

int SparseMatrix::Size_(const SparseMatrix &M, const vector<int> &rindex) const {
  int n = 0;
  for (int i = 0; i < M.size(); ++i)
    if (rindex[i] != -1)
      for (int k = M.rowind(i); k < M.rowind(i + 1); ++k)
        if (rindex[M.colptr(k)] != -1) ++n;
  return n;
}

int SparseMatrix::Size_(const SparseMatrix &M, const vector<int> &indexrow,
                        const vector<int> &indexcol) const {
  int nonzero = 0;
  vector<int> searchindex;
  searchindex.resize(M.size());
  for (int i = 0; i < searchindex.size(); ++i)
    searchindex[i] = -1;
  for (int i = 0; i < indexcol.size(); ++i)
    searchindex[indexcol[i]] = i;
  for (int i = 0; i < indexrow.size(); ++i)
    for (int k = M.rowind(indexrow[i]); k < M.rowind(indexrow[i] + 1); ++k)
      if (searchindex[M.colptr(k)] != -1) nonzero++;
  return nonzero;
}

int SparseMatrix::int_in_vector(int x, const vector<int> &y) const {
  if (x == -1) return -1;
  for (int i = 0; i < y.size(); ++i)
    if (y[i] == x) return i;
  return -1;
}

int SparseMatrix::_Size(const SparseMatrix &M, const vector<int> &Indices) const {
  int n = 0;
  for (int i = 0; i < Indices.size(); ++i)
    for (int k = M.rowind(Indices[i]); k < M.rowind(Indices[i] + 1); ++k)
      if (int_in_vector(M.colptr(k), Indices) != -1) n++;
  return n;
}

SparseMatrix::SparseMatrix(const Matrix &M) : SparseRMatrix(M.size(), M.Size()) {
  M.copy(nzval(), rowind(), colptr());

  SetComputedNcols();
}

SparseMatrix::SparseMatrix(const SparseMatrix &M) : SparseRMatrix(M.size(), NonZeros(M)) {
  const int n = size();

  Scalar *val = this->nzval();
  int *rowind = this->rowind();
  int *ind = this->colptr();

  int nz = 0;
  rowind[0] = 0;
  for (int r = 0; r < n; ++r) {
    for (int i = M.rowind(r); i < M.rowind(r + 1); ++i) {
      if (isZero(M.nzval(i))) continue;
      //             if (M.nzval(i) == 0.0) continue;
      ind[nz] = M.colptr(i);
      val[nz++] = M.nzval(i);
    }
    rowind[r + 1] = nz;
  }
}

SparseMatrix::SparseMatrix(const SparseMatrix &M, const vector<bool> &mask) :
    SparseRMatrix(size_(M, mask), Size_(M, mask)) {
  vector<int> index(M.size());
  int *d = rowind();
  int n = 0;
  d[0] = 0;
  for (int i = 0; i < M.size(); ++i)
    if (mask[i]) {
      index[i] = n;
      ++n;
      d[n] = d[n - 1];
      for (int k = M.rowind(i); k < M.rowind(i + 1); ++k) {
        int j = M.colptr(k);
        if (mask[j]) ++(d[n]);
      }
    }
  int m = 0;
  for (int i = 0; i < M.size(); ++i)
    if (mask[i])
      for (int k = M.rowind(i); k < M.rowind(i + 1); ++k) {
        int j = M.colptr(k);
        if (!mask[j]) continue;
        nzval(m) = M.nzval(k);
        colptr(m) = index[j];
        ++m;
      }
}

SparseMatrix::SparseMatrix(const SparseMatrix &M, const vector<int> &Indices) :
    SparseRMatrix(int(Indices.size()), _Size(M, Indices)) {
  int *d = rowind();
  int m = 0;
  d[0] = 0;
  for (int i = 0; i < Indices.size(); ++i) {
    d[i + 1] = d[i];
    for (int k = M.rowind(Indices[i]); k < M.rowind(Indices[i] + 1); ++k) {
      int akt = int_in_vector(M.colptr(k), Indices);
      if (akt != -1) {
        colptr(m) = akt;
        nzval(m) = M.nzval(k);
        ++d[i + 1];
        m++;
      }
    }
  }
}

SparseMatrix::SparseMatrix(const SparseMatrix &M, const vector<int> &index,
                           const vector<int> &rindex) :
    SparseRMatrix(int(index.size()), Size_(M, rindex)) {
  int *d = rowind();
  int n = 0;
  d[0] = 0;
  for (int i = 0; i < M.size(); ++i)
    if (rindex[i] != -1) {
      ++n;
      d[n] = d[n - 1];
      for (int k = M.rowind(i); k < M.rowind(i + 1); ++k) {
        int j = M.colptr(k);
        if (rindex[j] != -1) ++(d[n]);
      }
    }
  int m = 0;
  for (int i = 0; i < M.size(); ++i)
    if (rindex[i] != -1)
      for (int k = M.rowind(i); k < M.rowind(i + 1); ++k) {
        int j = M.colptr(k);
        if (rindex[j] == -1) continue;
        nzval(m) = M.nzval(k);
        colptr(m) = rindex[j];
        ++m;
      }
}

SparseMatrix::SparseMatrix(const SparseMatrix &M, const vector<int> &indexrow,
                           const vector<int> &indexcol, bool dummy) :
    SparseRMatrix(int(indexrow.size()), Size_(M, indexrow, indexcol)) {

  vector<int> searchindex;
  searchindex.resize(M.size());
  for (int i = 0; i < searchindex.size(); ++i)
    searchindex[i] = -1;
  for (int i = 0; i < indexcol.size(); ++i)
    searchindex[indexcol[i]] = i;

  int *d = rowind();
  Scalar *nzval_ = nzval();
  int *colptr_ = colptr();
  int mi = 0;
  d[0] = 0;
  for (int i = 0; i < size(); ++i) {
    d[i + 1] = d[i];
    for (int k = M.rowind(indexrow[i]); k < M.rowind(indexrow[i] + 1); ++k) {
      int akt = searchindex[M.colptr(k)];
      if (akt != -1) {
        colptr_[mi] = akt;
        nzval_[mi] = M.nzval(k);
        ++d[i + 1];
        mi++;
      }
    }
  }
}

void SparseMatrix::convert_sm(SparseMatT &M) const {
  const int n = size();
  const Scalar *val = this->nzval();
  const int *rowind = this->rowind();
  const int *ind = this->colptr();

  M.clear();
  M.resize(n);
  for (int r = 0; r < n; ++r) {
    SparseVecT &row = M[r];
    for (int i = rowind[r]; i < rowind[r + 1]; ++i)
      row[*(ind++)] = *(val++);
  }
}

void SparseMatrix::print_sm(const SparseMatT &M) const {
  for (int i = 0; i < M.size(); ++i) {
    const SparseVecT &r = M[i];
    mout << "row " << i << ": ";
    for (conSparVecIt it = r.begin(); it != r.end(); ++it)
      mout << '(' << it->first << ';' << it->second << "),";
    mout << endl;
  }
  mout << endl;
}

void SparseMatrix::print_sm_dense(const SparseMatT &M) const {
  printf("%3s ", "");
  for (int i = 0; i < M.size(); ++i)
    printf("%5i ", i);
  mout << endl;
  for (int i = 0; i < M.size(); ++i) {
    const SparseVecT &r = M[i];
    conSparVecIt rend = r.end();
    printf("%3i ", i);
    for (int j = 0; j < M.size(); ++j) {
      conSparVecIt it = r.find(j);
      printf("%5.2f ", (it == rend) ? 0.0 : it->second);
    }
    mout << endl;
  }
  mout << endl;
}

void SparseMatrix::convert_sm_back(const SparseMatT &M) {
  const int n = int(M.size());
  int rowpos = 0;
  Scalar *val = this->nzval();
  int *rowind = this->rowind();
  int *ind = this->colptr();

  for (int r = 0; r < n; ++r) {
    const SparseVecT &row = M[r];
    rowind[r] = rowpos;
    for (conSparVecIt it = row.begin(); it != row.end(); ++it) {
      *(val++) = it->second;
      *(ind++) = it->first;
      rowpos++;
    }
  }
  rowind[n] = rowpos;
}

void SparseMatrix::build_ident(const Matrix &u) {
  int nDof = u.size();

  ParDofId.resize(nDof);
  for (int i = 0; i < nDof; ++i)
    ParDofId[i] = i;
  for (identifyset is = u.identifysets(); is != u.identifysets_end(); ++is) {
    const Point &pmid = is();
    const int mid = u.Id(pmid);
    int i;
    for (i = 0; i < is.size(); ++i)
      if (pmid > is[i]) break;
    if (i == is.size()) {
      for (int i = 0; i < is.size(); ++i) {
        ParDofId[u.Id(is[i])] = mid;
      }
      nDof -= int(is.size());
    }
  }
  ParDofRes.resize(nDof);
  //   mout<<"rest DoFs="<<nDof<<" of "<<ParDofId.size()<<endl;
  //   mout<<"first run\n"<<ParDofId<<endl;
  int pos = 0;
  for (int i = 0; i < ParDofId.size(); ++i)
    if (ParDofId[i] == i) {
      ParDofId[i] = pos;
      ParDofRes[pos++] = i;
    } else if (ParDofId[i] < pos) ParDofId[i] = ParDofId[ParDofId[i]];
    else ParDofId[i] = -ParDofId[i];
  for (int i = 0; i < ParDofId.size(); ++i)
    if (ParDofId[i] < 0) ParDofId[i] = ParDofId[-ParDofId[i]];
}

void SparseMatrix::apply_PC_mat(SparseMatT &D, const SparseMatT &S) const {
  D.clear();
  D.resize(ParDofRes.size());
  for (int r = 0; r < S.size(); ++r) {
    const SparseVecT &srow = S[r];
    SparseVecT &drow = D[ParDofId[r]];
    for (conSparVecIt it = srow.begin(); it != srow.end(); ++it)
      drow[ParDofId[it->first]] += it->second;
  }
}

void SparseMatrix::multiply_plus(Vector &b, const Vector &u) const {
  plusMatVec(b, u);
  b.Collect();
}

void SparseMatrix::multiply_minus(Vector &b, const Vector &u) const {
  minusMatVec(b, u);
  b.Collect();
}

void SparseMatrix::plusMatVec(Vector &b, const Vector &u) const {
  SparseRMatrix::plusMatVec(b(), u());
}

void SparseMatrix::minusMatVec(Vector &b, const Vector &u) const {
  SparseRMatrix::minusMatVec(b(), u());
}

void SparseMatrix::pc_mat_convert(const Matrix &A) {
  build_ident(A);
  SparseMatT M, Mpc;
  convert_sm(M);
  apply_PC_mat(Mpc, M);
  resize(int(ParDofRes.size()), Size());
  convert_sm_back(Mpc);
}

void SparseMatrix::ExpandIdentify(Vector &u) const {
  Vector v(u);
  for (int i = 0; i < u.size(); ++i)
    u[i] = v[ParDofId[i]];
}

void SparseMatrix::ShrinkIdentify(Vector &u, const Vector &b) const {
  u = 0;
  for (int i = 0; i < u.size(); ++i)
    u[ParDofId[i]] += b[i];
}

Scalar *SparseMatrix::ref(int j) const {
  Scalar *a;
  a = new Scalar[rows()];
  for (int i = 0; i < rows(); ++i) {
    int e = find(i, j);
    if (e == -1) a[i] = 0.0;
    else a[i] = nzval(e);
  }
  return a;
}

constAB<Operator, Vector> operator*(const SparseMatrix &S, const Vector &v) {
  return constAB<Operator, Vector>(S, v);
}

constAB<Operator, Vectors> operator*(const SparseMatrix &S, const Vectors &v) {
  return constAB<Operator, Vectors>(S, v);
}
