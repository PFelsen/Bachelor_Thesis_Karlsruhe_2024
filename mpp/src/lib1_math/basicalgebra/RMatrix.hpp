#ifndef RMATRIX_H
#define RMATRIX_H

#include "AntisymRMatrix.hpp"
#include "SymRMatrix.hpp"

template<typename REAL = double>
class RMatrixT {
protected:
  std::vector<REAL> z{}; // rowwise entries
  int Nrows;
  int Ncols;
public:
  RMatrixT() : z(0, REAL{}), Nrows(0), Ncols(0) {}

  explicit constexpr RMatrixT(int dim) : z(dim * dim, REAL{}), Nrows(dim), Ncols(dim) {}

  explicit constexpr RMatrixT(const REAL &b, int dim) : z(dim * dim, b), Nrows(dim), Ncols(dim) {}

  constexpr RMatrixT(int rows, int cols) : z(rows * cols, REAL{}), Nrows(rows), Ncols(cols) {}

  constexpr RMatrixT(const REAL &b, int rows, int cols) :
      z(rows * cols, b), Nrows(rows), Ncols(cols) {}

  template<typename REAL1>
  RMatrixT(const std::initializer_list<std::initializer_list<REAL1>> &v) {
    Ncols = (int)(v.begin())->size();
    Nrows = (int)v.size();
    z.resize(Nrows * Ncols);
    for (int i = 0; i < Nrows; i++) {
      for (int j = 0; j < Ncols; j++) {
        z[i * Ncols + j] = ((v.begin() + i)->begin())[j];
      }
    }
  }

  template<typename REAL1>
  RMatrixT(const std::vector<std::vector<REAL1>> &v) {
    Nrows = (int)(v.begin())->size();
    Ncols = (int)v.size();
    z.resize(Nrows * Ncols);
    for (int i = 0; i < Nrows; i++) {
      for (int j = 0; j < Ncols; j++) {
        z[i * Ncols + j] = ((v.begin() + i)->begin())[j];
      }
    }
  }

  template<typename REAL1>
  explicit RMatrixT(const RMatrixT<REAL1> &A) :
      z(A.begin(), A.end()), Nrows(A.rows()), Ncols(A.cols()) {}

  template<typename REAL1>
  explicit RMatrixT(const SymRMatrixT<REAL1> &A) :
      z(A.Dim() * A.Dim()), Nrows(A.Dim()), Ncols(A.Dim()) {
    for (int i = 0; i < A.Dim(); ++i) {
      z[i * Ncols + i] = A(i, i);
      for (int j = 0; j < i; ++j) {
        z[i * Ncols + j] = A(i, j);
        z[j * Ncols + i] = A(i, j);
      }
    }
  }

  template<typename REAL1>
  explicit RMatrixT(const AntisymRMatrixT<REAL1> &A) :
      z(A.Dim() * A.Dim()), Nrows(A.Dim()), Ncols(A.Dim()) {
    for (int i = 0; i < A.Dim(); ++i)
      for (int j = 0; j < i; ++j) {
        z[i * Ncols + j] = A(i, j);
        z[j * Ncols + i] = A(j, i);
      }
  }

  template<typename MATRIX_A, typename MATRIX_B>
  RMatrixT(const MATRIX_A &A, const MATRIX_B &B) :
      z(A.rows() * B.cols(), REAL{}), Nrows(A.rows()), Ncols(B.cols()) {
    if (A.cols() != B.rows()) THROW("Size in matrix multiplication does not fit!")
    for (int i = 0; i < Nrows; ++i)
      for (int j = 0; j < Ncols; ++j)
        for (int k = 0; k < A.cols(); ++k)
          z[i * Ncols + j] += A(i, k) * B(k, j);
  }

  template<typename REAL1>
  explicit constexpr RMatrixT(const RVectorT<REAL1> &d) :
      z(d.size() * d.size(), REAL{}), Nrows(d.size()), Ncols(d.size()) {
    for (int i = 0; i < Nrows; ++i)
      z[i * Ncols + i] = d[i];
  }

  RMatrixT &operator=(const REAL &b) {
    for (int i = 0; i < z.size(); ++i)
      z[i] = b;
    return *this;
  }

  template<typename REAL1>
  RMatrixT &operator=(const RMatrixT<REAL1> &A) {
    Nrows = A.rows();
    Ncols = A.cols();
    z = std::vector<REAL>(A.begin(), A.end());
    return *this;
  }

  //    RMatrixT &operator=(RMatrixT &&A) {
  //        (*this) = std::move(A);
  //        return *this;
  //    }

  template<typename REAL1>
  RMatrixT &operator=(const SymRMatrixT<REAL1> &A) {
    resize(A.Dim());
    for (int i = 0; i < A.Dim(); ++i) {
      z[i * Ncols + i] = A(i, i);
      for (int j = 0; j < i; ++j) {
        z[i * Ncols + j] = A(i, j);
        z[j * Ncols + i] = A(i, j);
      }
    }
    return *this;
  }

  template<typename REAL1>
  RMatrixT &operator=(const AntisymRMatrixT<REAL1> &A) {
    resize(A.Dim());
    for (int i = 0; i < A.Dim(); ++i)
      for (int j = 0; j < i; ++j) {
        z[i * Ncols + j] = A(i, j);
        z[j * Ncols + i] = A(j, i);
      }
    return *this;
  }

  template<typename REAL1>
  RMatrixT &operator+=(const RMatrixT<REAL1> &A) {
    if (Nrows != A.rows() || Ncols != A.cols()) THROW("Size does not fit!")
    for (int i = 0; i < z.size(); ++i)
      z[i] += A.Data()[i];
    return *this;
  }

  template<typename REAL1>
  RMatrixT &operator+=(const SymRMatrixT<REAL1> &A) {
    if (Nrows != A.Dim() || Ncols != A.Dim()) THROW("Size does not fit!")
    for (int i = 0; i < Nrows; ++i)
      for (int j = 0; j < Ncols; ++j)
        z[i * Ncols + j] += A(i, j);
    return *this;
  }

  template<typename REAL1>
  RMatrixT &operator+=(const AntisymRMatrixT<REAL1> &A) {
    if (Nrows != A.Dim() || Ncols != A.Dim()) THROW("Size does not fit!")
    for (int i = 0; i < Nrows; ++i)
      for (int j = 0; j < Ncols; ++j)
        z[i * Ncols + j] += A(i, j);
    return *this;
  }

  template<typename REAL1>
  RMatrixT &operator-=(const RMatrixT<REAL1> &A) {
    if (Nrows != A.rows() || Ncols != A.cols()) THROW("Size does not fit!")
    for (int i = 0; i < z.size(); ++i)
      z[i] -= A.Data()[i];
    return *this;
  }

  template<typename REAL1>
  RMatrixT &operator-=(const SymRMatrixT<REAL1> &A) {
    if (Nrows != A.Dim() || Ncols != A.Dim()) THROW("Size does not fit!")
    for (int i = 0; i < Nrows; ++i)
      for (int j = 0; j < Ncols; ++j)
        z[i * Ncols + j] -= A(i, j);
    return *this;
  }

  template<typename REAL1>
  RMatrixT &operator-=(const AntisymRMatrixT<REAL1> &A) {
    if (Nrows != A.Dim() || Ncols != A.Dim()) THROW("Size does not fit!")
    for (int i = 0; i < Nrows; ++i)
      for (int j = 0; j < Ncols; ++j)
        z[i * Ncols + j] -= A(i, j);
    return *this;
  }

  RMatrixT &operator*=(const REAL &b) {
    for (int i = 0; i < z.size(); ++i)
      z[i] *= b;
    return *this;
  }

  RMatrixT &operator/=(const REAL &b) {
    for (int i = 0; i < z.size(); ++i)
      z[i] /= b;
    return *this;
  }

  auto operator[](int i) { return z.begin() + i * Ncols; }

  auto operator[](int i) const { return z.begin() + i * Ncols; }

  REAL &operator()(int row, int col) { return z[row * Ncols + col]; }

  REAL &operator()(int index) { return z[index]; }

  const REAL &operator()(int row, int col) const { return z[row * Ncols + col]; }

  void operator()(const REAL &b, int row, int col) { z[row * Ncols + col] = b; }

  RVectorT<REAL> row(int i) const {
    RVectorT<REAL> r(Ncols);
    for (int j = 0; j < Ncols; ++j)
      r[j] = z[i * Ncols + j];
    return r;
  }

  RVectorT<REAL> col(int j) const {
    RVectorT<REAL> c(Nrows);
    for (int i = 0; i < Nrows; ++i)
      c[i] = z[i * Ncols + j];
    return c;
  }

  int size() const { return z.size(); }

  int rows() const { return Nrows; }

  int cols() const { return Ncols; }

  auto begin() const { return z.begin(); }

  auto begin() { return z.begin(); }

  auto end() const { return z.end(); }

  auto end() { return z.end(); }

  const std::vector<REAL> &Data() const { return z; }

  std::vector<REAL> &Data() { return z; }

  REAL normSqr() const;

  REAL norm() const;

  RVectorT<REAL> diag();

  template<typename REAL1>
  void diag(const RVectorT<REAL1> &d) {
    resize(d.size(), d.size());
    for (int i = 0; i < d.size(); ++i)
      z[i * Ncols + i] = d[i];
  }

  template<typename REAL1>
  RVectorT<REAL> multiplyWith(const RVectorT<REAL1> &v) const {
    if (v.Dim() != Ncols) {
      THROW("Dimension does not fit: "
            "Matcols: "
            + std::to_string(Ncols)
            + " "
              "Vectordim: "
            + std::to_string(v.Dim()))
    }
    RVectorT<REAL> w(Nrows);
    for (int i = 0; i < Nrows; ++i) {
      for (int j = 0; j < Ncols; ++j)
        w[i] += z[i * Ncols + j] * v[j];
    }
    return w;
  }

  template<typename REAL1>
  RVectorT<REAL> multipliedWith(const RVectorT<REAL1> &v) const {
    if (v.Dim() != Nrows) THROW("Dimension does not fit!")
    RVectorT<REAL> w(Ncols);
    for (int i = 0; i < Ncols; ++i) {
      w[i] = v * col(i);
    }
    return w;
  }

  void resize(int rows, int cols);

  void resize(int dim) { resize(dim, dim); }

  void Accumulate(int commSplit = 0);

  template<typename REAL1>
  void Insert(const RMatrixT<REAL1> &B, int row, int col) {
    if (row + B.rows() >= Nrows || col + B.cols() >= Ncols) {
      THROW("Cannot insert matrix: size of B too large")
    }
    for (int i = 0; i < B.rows(); ++i)
      for (int j = 0; j < B.cols(); ++j)
        z[(row + i) * Ncols + col + j] = B[i][j];
  }

  template<typename REAL1>
  void InsertRow(const RVectorT<REAL1> &v, int row, int startCol = 0) {
    if (row >= Nrows || startCol + v.size() > Ncols) {
      THROW("Cannot insert row: size of v too large")
    }
    for (int j = 0; j < v.size(); ++j)
      z[row * Ncols + startCol + j] = v[j];
  }

  template<typename REAL1>
  void InsertCol(const RVectorT<REAL1> &v, int col, int startRow = 0) {
    if (col >= Ncols || startRow + v.size() > Nrows) {
      THROW("Cannot insert column: size of v too large")
    }
    for (int i = 0; i < v.size(); ++i)
      z[(startRow + i) * Ncols + col] = v[i];
  }

  REAL Mean() const;

  REAL Variance() const;

  REAL Trace() const;

  REAL NormOne() const;

  REAL NormInfty() const;

  RMatrixT<REAL> &transpose();

  RMatrixT<REAL> &Identity();

  RMatrixT<REAL> &Invert();

  RMatrixT<REAL> &Exp(bool deg13 = false);

  RMatrixT<REAL> &Phi1(RVector &Evaluated, bool deg13 = false);

  RMatrixT<REAL> &Phi1();

  RMatrixT<REAL> &Phi2();

  RMatrixT<REAL> &Phi3();


  void SaddlePoint(const RMatrixT &A, const RMatrixT &B);

  Saver &save(Saver &saver) const;

  Loader &load(Loader &loader);

  RMatrixT<double> &Exp(RVector &Evaluated, bool deg13 = false);
};

template<typename REAL>
bool operator==(const RMatrixT<REAL> &A, const RMatrixT<REAL> &B) {
  if (A.rows() != B.rows() || A.cols() != B.cols()) return false;
  for (int i = 0; i < A.Data().size(); ++i)
    if (!mpp_ba::isNear(A.Data()[i], B.Data()[i])) return false;
  return true;
}

template<typename REAL>
bool operator==(const RMatrixT<REAL> &A, const SymRMatrixT<REAL> &B) {
  return A == RMatrixT<REAL>(B);
}

template<typename REAL>
bool operator==(const SymRMatrixT<REAL> &A, const RMatrixT<REAL> &B) {
  return RMatrixT<REAL>(A) == B;
}

template<typename REAL>
bool operator==(const RMatrixT<REAL> &A, const AntisymRMatrixT<REAL> &B) {
  return A == RMatrixT<REAL>(B);
}

template<typename REAL>
bool operator==(const AntisymRMatrixT<REAL> &A, const RMatrixT<REAL> &B) {
  return RMatrixT<REAL>(A) == B;
}

template<typename REAL>
bool operator!=(const RMatrixT<REAL> &A, const RMatrixT<REAL> &B) {
  return !(A == B);
}

template<typename REAL>
bool operator!=(const RMatrixT<REAL> &A, const SymRMatrixT<REAL> &B) {
  return !(A == B);
}

template<typename REAL>
bool operator!=(const SymRMatrixT<REAL> &A, const RMatrixT<REAL> &B) {
  return !(A == B);
}

template<typename REAL>
bool operator!=(const RMatrixT<REAL> &A, const AntisymRMatrixT<REAL> &B) {
  return !(A == B);
}

template<typename REAL>
bool operator!=(const AntisymRMatrixT<REAL> &A, const RMatrixT<REAL> &B) {
  return !(A == B);
}

template<typename REAL>
inline RMatrixT<REAL> operator-(const RMatrixT<REAL> &a) {
  RMatrixT<REAL> c(a);
  return c *= -1;
}

template<typename REAL>
inline RMatrixT<REAL> operator+(const RMatrixT<REAL> &a, const RMatrixT<REAL> &b) {
  RMatrixT<REAL> c(a);
  return c += b;
}

template<typename REAL>
inline RMatrixT<REAL> operator+(const RMatrixT<REAL> &a, const SymRMatrixT<REAL> &b) {
  RMatrixT<REAL> c(a);
  return c += b;
}

template<typename REAL>
inline RMatrixT<REAL> operator+(const SymRMatrixT<REAL> &a, const RMatrixT<REAL> &b) {
  RMatrixT<REAL> c(a);
  return c += b;
}

template<typename REAL>
inline RMatrixT<REAL> operator+(const RMatrixT<REAL> &a, const AntisymRMatrixT<REAL> &b) {
  RMatrixT<REAL> c(a);
  return c += b;
}

template<typename REAL>
inline RMatrixT<REAL> operator+(const AntisymRMatrixT<REAL> &a, const RMatrixT<REAL> &b) {
  RMatrixT<REAL> c(a);
  return c += b;
}

template<typename REAL>
inline RMatrixT<REAL> operator-(const RMatrixT<REAL> &a, const RMatrixT<REAL> &b) {
  RMatrixT<REAL> c(a);
  return c -= b;
}

template<typename REAL>
inline RMatrixT<REAL> operator-(const RMatrixT<REAL> &a, const SymRMatrixT<REAL> &b) {
  RMatrixT<REAL> c(a);
  return c -= b;
}

template<typename REAL>
inline RMatrixT<REAL> operator-(const SymRMatrixT<REAL> &a, const RMatrixT<REAL> &b) {
  RMatrixT<REAL> c(a);
  return c -= b;
}

template<typename REAL>
inline RMatrixT<REAL> operator-(const RMatrixT<REAL> &a, const AntisymRMatrixT<REAL> &b) {
  RMatrixT<REAL> c(a);
  return c -= b;
}

template<typename REAL>
inline RMatrixT<REAL> operator-(const AntisymRMatrixT<REAL> &a, const RMatrixT<REAL> &b) {
  RMatrixT<REAL> c(a);
  return c -= b;
}

template<typename REAL>
inline RMatrixT<REAL> operator*(const RMatrixT<REAL> &a, const RMatrixT<REAL> &b) {
  return RMatrixT<REAL>(a, b);
}

template<typename REAL>
inline RMatrixT<REAL> operator*(const SymRMatrixT<REAL> &a, const SymRMatrixT<REAL> &b) {
  return RMatrixT<REAL>(a, b);
}

template<typename REAL, typename REAL1>
inline RMatrixT<REAL> operator*(const RMatrixT<REAL> &a, const SymRMatrixT<REAL1> &b) {
  return RMatrixT<REAL>(a, b);
}

template<typename REAL, typename REAL1>
inline RMatrixT<REAL> operator*(const SymRMatrixT<REAL1> &a, const RMatrixT<REAL> &b) {
  return RMatrixT<REAL>(a, b);
}

template<typename REAL>
inline RMatrixT<REAL> operator*(const AntisymRMatrixT<REAL> &a, const AntisymRMatrixT<REAL> &b) {
  return RMatrixT<REAL>(a, b);
}

template<typename REAL, typename REAL1>
inline RMatrixT<REAL> operator*(const RMatrixT<REAL> &a, const AntisymRMatrixT<REAL1> &b) {
  return RMatrixT<REAL>(a, b);
}

template<typename REAL, typename REAL1>
inline RMatrixT<REAL> operator*(const AntisymRMatrixT<REAL1> &a, const RMatrixT<REAL> &b) {
  return RMatrixT<REAL>(a, b);
}

template<typename REAL>
RVectorT<REAL> operator*(const RMatrixT<REAL> &a, const RVectorT<REAL> &v) {
  return a.multiplyWith(v);
}

template<typename REAL>
RVectorT<REAL> operator*(const RVectorT<REAL> &v, const RMatrixT<REAL> &a) {
  return a.multipliedWith(v);
}

template<typename REAL, typename REAL1>
inline RMatrixT<REAL> operator*(const RMatrixT<REAL> &a, const REAL1 &b) {
  RMatrixT<REAL> c(a);
  return c *= b;
}

template<typename REAL, typename REAL1>
inline RMatrixT<REAL> operator*(const REAL1 &b, const RMatrixT<REAL> &a) {
  RMatrixT<REAL> c(a);
  return c *= b;
}

template<typename REAL, typename REAL1>
inline RMatrixT<REAL> operator/(const RMatrixT<REAL> &a, const REAL1 &b) {
  RMatrixT<REAL> c(a);
  return c /= b;
}

template<typename REAL>
inline RMatrixT<REAL> transpose(const RMatrixT<REAL> &a) {
  RMatrixT<REAL> b(a);
  return b.transpose();
}

template<typename REAL>
inline RMatrixT<REAL> invert(const RMatrixT<REAL> &a) {
  RMatrixT<REAL> b(a);
  return b.Invert();
}

template<typename REAL>
inline RMatrixT<REAL> product(const RVectorT<REAL> &v, const RVectorT<REAL> &w) {
  RMatrixT<REAL> A(v.size(), w.size());
  for (int r = 0; r < v.size(); ++r) {
    for (int c = 0; c < w.size(); ++c) {
      A(r, c) = v[r] * w[c];
    }
  }
  return A;
}

template<typename REAL>
inline Saver &operator<<(Saver &saver, const RMatrixT<REAL> &A) {
  return A.save(saver);
}

template<typename REAL>
inline Loader &operator>>(Loader &loader, RMatrixT<REAL> &A) {
  return A.load(loader);
}

template<typename REAL>
std::ostream &operator<<(std::ostream &os, const RMatrixT<REAL> &A) {
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

// LAPACK-like Methods because the real LAPACK seems to slow down the calculations...
void LUDecomp(int n, double *Mat, int *ipv);
void invertSmallMatrix(int n, double *a);
void FwdBwdSubstitution(int n, double *Mat, int *piv, double *rhs);
void applySmallMatrix(int n, double *u, const double *a, const double *b);

using RMatrix = RMatrixT<double>;

#ifdef BUILD_IA

using IARMatrix = RMatrixT<IAInterval>;

RMatrix mid(const IARMatrix &);

RMatrix sup(const IARMatrix &);

RMatrix inf(const IARMatrix &);

inline IARMatrix operator+(const IARMatrix &A, const RMatrix &B) {
  IARMatrix C(A);
  return C += B;
}

inline IARMatrix operator+(const RMatrix &A, const IARMatrix &B) {
  IARMatrix C(A);
  return C += B;
}

inline IARMatrix operator+(const IARMatrix &A, const SymRMatrix &B) {
  IARMatrix C(A);
  return C += B;
}

inline IARMatrix operator+(const SymRMatrix &A, const IARMatrix &B) {
  IARMatrix C(A);
  return C += B;
}

inline IARMatrix operator+(const IARMatrix &A, const AntisymRMatrix &B) {
  IARMatrix C(A);
  return C += B;
}

inline IARMatrix operator+(const AntisymRMatrix &A, const IARMatrix &B) {
  IARMatrix C(A);
  return C += B;
}

inline IARMatrix operator-(const IARMatrix &A, const RMatrix &B) {
  IARMatrix C(A);
  return C -= B;
}

inline IARMatrix operator-(const RMatrix &A, const IARMatrix &B) {
  IARMatrix C(A);
  return C -= B;
}

inline IARMatrix operator-(const IARMatrix &A, const SymRMatrix &B) {
  IARMatrix C(A);
  return C -= B;
}

inline IARMatrix operator-(const SymRMatrix &A, const IARMatrix &B) {
  IARMatrix C(A);
  return C -= B;
}

inline IARMatrix operator-(const IARMatrix &A, const AntisymRMatrix &B) {
  IARMatrix C(A);
  return C -= B;
}

inline IARMatrix operator-(const AntisymRMatrix &A, const IARMatrix &B) {
  IARMatrix C(A);
  return C -= B;
}

inline IARVector operator*(const IARMatrix &A, const RVector &v) { return A.multiplyWith(v); }

inline IARVector operator*(const RMatrix &A, const IARVector &v) {
  IARMatrix B(A);
  return B.multiplyWith(v);
}

inline IARMatrix operator*(const IARMatrix &a, const RMatrix &b) { return IARMatrix(a, b); }

inline IARMatrix operator*(const RMatrix &a, const IARMatrix &b) { return IARMatrix(a, b); }

inline IARMatrix operator*(const IARMatrix &a, const SymRMatrix &b) { return IARMatrix(a, b); }

inline IARMatrix operator*(const SymRMatrix &a, const IARMatrix &b) { return IARMatrix(a, b); }

inline IARMatrix operator*(const IASymRMatrix &a, const RMatrix &b) { return IARMatrix(a, b); }

inline IARMatrix operator*(const RMatrix &a, const IASymRMatrix &b) { return IARMatrix(a, b); }

inline IARMatrix operator*(const IARMatrix &a, const AntisymRMatrix &b) { return IARMatrix(a, b); }

inline IARMatrix operator*(const AntisymRMatrix &a, const IARMatrix &b) { return IARMatrix(a, b); }

inline IARMatrix operator*(const IAAntisymRMatrix &a, const RMatrix &b) { return IARMatrix(a, b); }

inline IARMatrix operator*(const RMatrix &a, const IAAntisymRMatrix &b) { return IARMatrix(a, b); }

#endif // BUILD_IA

#endif
