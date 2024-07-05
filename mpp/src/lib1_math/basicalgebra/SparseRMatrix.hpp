#ifndef SPARSERMATRIX_HPP
#define SPARSERMATRIX_HPP

#include "GlobalDefinitions.hpp"
#include "RMatrix.hpp"


#ifdef BUILD_IA

#include "IAInterval.hpp"

#endif

template<typename REAL = double>
class SparseRMatrixT {
  /**
   * Uses Compressed sparse row (CSR) format to store the data
   * Note that all diagonal entries are stored (independent of their value) at the
   * the beginning of each row
   */
  /// contains entries (all diagonal entries and non zero off diagonal entries)
  std::vector<REAL> a;
  /// contains indices (wrt a) of the diagonal entries, i.e. a new row starts at each index
  std::vector<int> d;
  /// contains the column index of the corresponding entry (in a)
  std::vector<int> col;
  /// number of rows
  int Nrows;
  /// number of colums
  int Ncols;

  int computeNcols() const;

  template<typename REAL1>
  int computeNumberOfEntries(const RMatrixT<REAL1> &A) {
    int numEntries = std::min(A.rows(), A.cols());
    for (int i = 0; i < A.rows(); ++i)
      for (int j = 0; j < A.cols(); ++j)
        if (i != j)
          if (!mpp_ba::isNearZero(A[i][j])) ++numEntries;
    return numEntries;
  }
protected:
  void SetComputedNcols();
public:
  template<typename REAL1>
  SparseRMatrixT(const SparseRMatrixT<REAL1> &A) :
      a(A.Values()), d(A.Indices()), col(A.Columns()), Nrows(A.rows()), Ncols(A.cols()) {}

  template<typename REAL1>
  SparseRMatrixT(const RMatrixT<REAL1> &A) :
      a(computeNumberOfEntries(A)), d(A.rows() + 1), col(a.size()), Nrows(A.rows()),
      Ncols(A.cols()) {
    int cnt = 0;
    for (int i = 0; i < Nrows; ++i) {
      d[i] = cnt;
      if (i < Ncols) {
        a[cnt] = A[i][i];
        col[cnt] = i;
        ++cnt;
      }
      for (int j = 0; j < Ncols; ++j)
        if (i != j)
          if (!mpp_ba::isNearZero(A[i][j])) {
            a[cnt] = A[i][j];
            col[cnt] = j;
            ++cnt;
          }
    }
    d[Nrows] = cnt;
  }

  // TODO: should be deleted or at least protected
  SparseRMatrixT(int rows, int numEntries, int cols = -1) :
      a(numEntries, REAL{}), d(rows + 1), col(numEntries), Nrows(rows), Ncols(cols) {}

  RMatrixT<REAL> ToRMatrix() const;

  REAL operator()(int i, int j) const;

  void SetZero();

  void resize(int rows, int numEntries);

  const std::vector<REAL> &Values() const { return a; }

  const std::vector<int> &Indices() const { return d; }

  const std::vector<int> &Columns() const { return col; }

  int rows() const { return Nrows; }

  int cols() const { return Ncols; }

  int NumberOfEntries() const { return d[Nrows]; }

  RVectorT<REAL> diag() const;

  template<typename REAL1>
  void diag(const RVectorT<REAL1> &r) {
    if (r.size() != std::min(Ncols, Nrows)) THROW("Size does not fit")
    for (int i = 0; i < r.size(); ++i)
      d[d[i]] = r[i];
  }

  void CheckDiagonal();

  void plusMatVec(REAL *b, const REAL *u) const;

  void minusMatVec(REAL *b, const REAL *u) const;

  void compress();

  int find(int i, int j) const {
    for (int e = d[i]; e < d[i + 1]; ++e)
      if (j == col[e]) return e;
    return -1;
  }

  template<typename VECTOR>
  void GaussSeidel(VECTOR &u, const VECTOR &b, bool shift = false) const {
    for (int i = 0; i < Nrows; ++i) {
      double r = b[i], diag = 0;
      for (int k = d[i]; k < d[i + 1]; ++k) {
        int j = col[k];
        if (j < i) r -= a[k] * u[j];
        else if (j == i) diag = a[k];
      }
      if ((shift) && (diag == 0.0)) diag = 1;
      u[i] = r / diag;
    }
  }

  template<typename VECTOR>
  void BackwardGaussSeidel(VECTOR &u, const VECTOR &b, bool shift = false) const {
    for (int i = Nrows - 1; i >= 0; --i) {
      double r = b[i], diag = 0;
      for (int k = d[i + 1] - 1; k >= d[i]; --k) {
        int j = col[k];
        if (j > i) r -= a[k] * u[j];
        else if (j == i) diag = a[k];
      }
      if ((shift) && (diag == 0.0)) diag = 1;
      u[i] = r / diag;
    }
  }

  //================================================================================================
  // Deprecated code: Should be deleted or at least protected
  //================================================================================================
  const REAL *nzval() const { return a.data(); }

  REAL *nzval() { return a.data(); }

  REAL nzval(int i) const { return a[i]; }

  REAL &nzval(int i) { return a[i]; }

  const int *rowind() const { return d.data(); }

  int *rowind() { return d.data(); }

  int rowind(int i) const { return d[i]; }

  int &rowind(int i) { return d[i]; }

  const int *colptr() const { return col.data(); }

  int *colptr() { return col.data(); }

  int colptr(int i) const { return col[i]; }

  int &colptr(int i) { return col[i]; }

  int size() const { return rows(); } // TODO: remove

  int Size() const { return NumberOfEntries(); } // TODO: remove

  void changeToDiagonal(int j);
};

template<typename REAL>
std::ostream &operator<<(std::ostream &os, const SparseRMatrixT<REAL> &A) {
  if constexpr (std::is_same_v<REAL, double>) os << beginD;
  for (int i = 0; i < A.rows(); ++i) {
    if (i != 0) {
      os << "\n";
      if constexpr (std::is_same_v<REAL, double>) os << beginD;
    }
    for (int j = 0; j < A.cols(); ++j)
      os << A(i, j) << " ";
  }
  if constexpr (std::is_same_v<REAL, double>) os << endD;
  return os;
}

using SparseRMatrix = SparseRMatrixT<double>;

#endif // SPARSERMATRIX_HPP
