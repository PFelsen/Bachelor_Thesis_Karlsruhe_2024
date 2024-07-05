#include "SparseRMatrix.hpp"

#include <cmath>

template<typename REAL>
void SparseRMatrixT<REAL>::SetZero() {
  int m = d[Nrows];
  for (int i = 0; i < m; ++i) {
    a[i] = REAL{};
    col[i] = 0;
  }
  for (int i = 0; i <= Nrows; ++i)
    d[i] = 0;
  d[Nrows] = m;
}

template<typename REAL>
void SparseRMatrixT<REAL>::resize(int rows, int numEntries) {
  Nrows = rows;
  Ncols = -1;
  a = std::vector<REAL>(numEntries);
  d = std::vector<int>(rows + 1);
  col = std::vector<int>(numEntries);
  d[rows] = numEntries;
}

template<typename REAL>
RVectorT<REAL> SparseRMatrixT<REAL>::diag() const {
  RVectorT<REAL> r(std::min(Nrows, Ncols));
  for (int i = 0; i < r.size(); ++i)
    r[i] = a[d[i]];
  return r;
}

template<typename REAL>
RMatrixT<REAL> SparseRMatrixT<REAL>::ToRMatrix() const {
  RMatrixT<REAL> A(Nrows, Ncols);
  int cnt = 0;
  for (int i = 0; i < Nrows; ++i)
    for (; cnt < d[i + 1]; ++cnt)
      A[i][col[cnt]] = a[cnt];
  return A;
}

template<typename REAL>
REAL SparseRMatrixT<REAL>::operator()(int i, int j) const {
  for (int cnt = d[i]; cnt < d[i + 1]; ++cnt) {
    if (col[cnt] == j) return a[cnt];
  }
  return REAL{};
}

template<typename REAL>
int SparseRMatrixT<REAL>::computeNcols() const {
  int cols = -1;
  for (int j = 0; j < d[Nrows]; ++j)
    cols = std::max(cols, col[j]);
  return ++cols;
}

template<typename REAL>
void SparseRMatrixT<REAL>::SetComputedNcols() {
  Ncols = computeNcols();
}

template<typename REAL>
void SparseRMatrixT<REAL>::CheckDiagonal() {
  for (int i = 0; i < Nrows; ++i)
    if (mpp_ba::isNearZero(a[d[i]])) a[d[i]] = 1.0;
}

template<typename REAL>
void SparseRMatrixT<REAL>::plusMatVec(REAL *b, const REAL *u) const {
  REAL val{};
  for (int i = 0; i < Nrows; ++i) {
    val = 0.0;
    for (int k = d[i]; k < d[i + 1]; ++k) {
      val += a[k] * u[col[k]];
    }
    b[i] += val;
  }
}

template<typename REAL>
void SparseRMatrixT<REAL>::minusMatVec(REAL *b, const REAL *u) const {
  REAL val{};
  for (int i = 0; i < Nrows; ++i) {
    val = 0.0;
    for (int k = d[i]; k < d[i + 1]; ++k) {
      val += a[k] * u[col[k]];
    }
    b[i] -= val;
  }
}

template<typename REAL>
void SparseRMatrixT<REAL>::compress() {
  int numEntries = std::min(Nrows, Ncols);
  for (int i = 0; i < Nrows; ++i)
    for (int k = d[i] + 1; k < d[i + 1]; ++k)
      if (!mpp_ba::isNearZero(a[k])) ++numEntries;
  std::vector<REAL> a0(numEntries);
  std::vector<int> d0(Nrows + 1);
  std::vector<int> col0(numEntries);
  d0[Nrows] = numEntries;

  int nn = 0;
  for (int i = 0; i < Nrows; ++i) {
    int k = d[i];
    d0[i] = nn;
    a0[nn] = a[k];
    col0[nn] = col[k];
    ++nn;
    for (int k = d[i] + 1; k < d[i + 1]; ++k)
      if (!mpp_ba::isNearZero(a[k])) {
        a0[nn] = a[k];
        col0[nn] = col[k];
        ++nn;
      }
  }
  //  mout << "Compressed matrix. Before: " << d[Nrows] << ", after: " << nn << "."
  //       << " Ratio: " << double(nn) / d[Nrows] << endl;

  a = std::move(a0);
  d = std::move(d0);
  col = std::move(col0);
}

template<typename REAL>
void SparseRMatrixT<REAL>::changeToDiagonal(int j) {
  std::vector<REAL> aa;
  aa.resize(a.size());
  std::vector<int> dd;
  dd.resize(d.size());
  std::vector<int> colcol;
  colcol.resize((col.size()));
  for (int i = 0; i < Nrows; ++i) {
    int k = (*this).colptr((*this).rowind(i));
    REAL val = (*this).nzval((*this).rowind(i));
    bool sorted = false;
    for (int l = (*this).rowind(i) + 1; l < (*this).rowind(i + 1); ++l) {
      if (colptr(l) == i + j) {
        colcol[this->rowind(i)] = i + j;
        aa[this->rowind(i)] = a[l];
        colcol[l] = colptr(l - 1);
        aa[l] = a[l - 1];
      } else {
        if (k > colptr(l)) {
          colcol[l] = col[l];
          aa[l] = a[l];
        } else if (k < colptr(l) && colptr(l) < (i + j)) {
          if (!sorted) {
            colcol[l] = k;
            aa[l] = a[(*this).rowind(i)];
            sorted = true;
          } else {
            colcol[l] = colptr(l - 1);
            aa[l] = a[l - 1];
          }
        } else {
          colcol[l] = colptr(l);
          aa[l] = a[l];
        }
      }
    }
  }
  a = aa;
  col = colcol;
}


template class SparseRMatrixT<double>;

#ifdef BUILD_IA
template class SparseRMatrixT<IAInterval>;

#endif
