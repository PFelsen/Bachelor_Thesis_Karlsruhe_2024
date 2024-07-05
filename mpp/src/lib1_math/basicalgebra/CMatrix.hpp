#ifndef CMATRIX_H
#define CMATRIX_H

#include "HermCMatrix.hpp"
#include "RMatrix.hpp"

template<typename REAL = double>
class CMatrixT {
  using COMPLEX = COMPLEX_TYPE<REAL>;
protected:
  std::vector<COMPLEX> z; // rowwise entries
  int Nrows;
  int Ncols;
public:
  CMatrixT() : z(0, COMPLEX{}), Nrows(0), Ncols(0) {}

  explicit constexpr CMatrixT(int dim) : z(dim * dim, COMPLEX{}), Nrows(dim), Ncols(dim) {}

  constexpr CMatrixT(int rows, int cols) : z(rows * cols, COMPLEX{}), Nrows(rows), Ncols(cols) {}

  constexpr CMatrixT(const COMPLEX &b, int rows, int cols) :
      z(rows * cols, b), Nrows(rows), Ncols(cols) {}

  template<typename COMPLEX1>
  CMatrixT(const std::initializer_list<std::initializer_list<COMPLEX1>> &v) {
    Nrows = (int)(v.begin())->size();
    Ncols = (int)v.size();
    z.resize(Nrows * Ncols);
    for (int i = 0; i < Nrows; i++) {
      for (int j = 0; j < Ncols; j++) {
        z[i * Ncols + j] = ((v.begin() + i)->begin())[j];
      }
    }
  }

  template<typename COMPLEX1>
  CMatrixT(const std::vector<std::vector<COMPLEX1>> &v) {
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
  explicit CMatrixT(const CMatrixT<REAL1> &A) :
      z(A.begin(), A.end()), Nrows(A.rows()), Ncols(A.cols()) {}

  template<typename REAL1>
  explicit CMatrixT(const RMatrixT<REAL1> &A) :
      z(A.begin(), A.end()), Nrows(A.rows()), Ncols(A.cols()) {}

  template<typename REAL1>
  explicit CMatrixT(const HermCMatrixT<REAL1> &A) :
      z(A.Dim() * A.Dim(), COMPLEX{}), Nrows(A.Dim()), Ncols(A.Dim()) {
    for (int i = 0; i < A.Dim(); ++i) {
      z[i * Ncols + i] = A(i, i);
      for (int j = 0; j < i; ++j) {
        z[i * Ncols + j] = A(i, j);
        z[j * Ncols + i] = A(j, i);
      }
    }
  }

  template<typename REAL1>
  explicit CMatrixT(const SymRMatrixT<REAL1> &A) :
      z(A.Dim() * A.Dim(), COMPLEX{}), Nrows(A.Dim()), Ncols(A.Dim()) {
    for (int i = 0; i < A.Dim(); ++i) {
      z[i * Ncols + i] = A(i, i);
      for (int j = 0; j < i; ++j) {
        z[i * Ncols + j] = A(i, j);
        z[j * Ncols + i] = A(j, i);
      }
    }
  }

  template<typename REAL1>
  explicit CMatrixT(const AntisymRMatrixT<REAL1> &A) :
      z(A.Dim() * A.Dim(), COMPLEX{}), Nrows(A.Dim()), Ncols(A.Dim()) {
    for (int i = 0; i < A.Dim(); ++i)
      for (int j = 0; j < i; ++j) {
        z[i * Ncols + j] = A(i, j);
        z[j * Ncols + i] = A(j, i);
      }
  }

  template<typename MATRIX_A, typename MATRIX_B>
  CMatrixT(const MATRIX_A &A, const MATRIX_B &B) :
      z(A.rows() * B.cols(), COMPLEX{}), Nrows(A.rows()), Ncols(B.cols()) {
    if (A.cols() != B.rows()) THROW("Size in matrix multiplication does not fit!")
    for (int i = 0; i < Nrows; ++i)
      for (int j = 0; j < Ncols; ++j)
        for (int k = 0; k < A.cols(); ++k)
          z[i * Ncols + j] += A(i, k) * B(k, j);
  }

  template<typename REAL1>
  explicit CMatrixT(const CVectorT<REAL1> &d) :
      z(d.size() * d.size(), COMPLEX{}), Nrows(d.size()), Ncols(d.size()) {
    for (int i = 0; i < Nrows; ++i)
      z[i * Ncols + i] = d[i];
  }

  template<typename REAL1>
  explicit CMatrixT(const RVectorT<REAL1> &d) :
      z(d.size() * d.size(), COMPLEX{}), Nrows(d.size()), Ncols(d.size()) {
    for (int i = 0; i < Nrows; ++i)
      z[i * Ncols + i] = d[i];
  }

  CMatrixT &operator=(const COMPLEX &b) {
    for (int i = 0; i < z.size(); ++i)
      z[i] = b;
    return *this;
  }

  template<typename REAL1>
  CMatrixT &operator=(const CMatrixT<REAL1> &A) {
    Nrows = A.rows();
    Ncols = A.cols();
    z = std::vector<COMPLEX>(A.begin(), A.end());
    return *this;
  }

  template<typename REAL1>
  CMatrixT &operator=(const RMatrixT<REAL1> &A) {
    Nrows = A.rows();
    Ncols = A.cols();
    z = std::vector<COMPLEX>(A.begin(), A.end());
    return *this;
  }

  template<typename REAL1>
  CMatrixT &operator=(const HermCMatrixT<REAL1> &A) {
    resize(A.Dim());
    for (int i = 0; i < A.Dim(); ++i) {
      z[i * Ncols + i] = A(i, i);
      for (int j = 0; j < i; ++j) {
        z[i * Ncols + j] = A(i, j);
        z[j * Ncols + i] = A(j, i);
      }
    }
    return *this;
  }

  template<typename REAL1>
  CMatrixT &operator=(const SymRMatrixT<REAL1> &A) {
    resize(A.Dim());
    for (int i = 0; i < A.Dim(); ++i) {
      z[i * Ncols + i] = A(i, i);
      for (int j = 0; j < i; ++j) {
        z[i * Ncols + j] = A(i, j);
        z[j * Ncols + i] = A(j, i);
      }
    }
    return *this;
  }

  template<typename REAL1>
  CMatrixT &operator=(const AntisymRMatrixT<REAL1> &A) {
    resize(A.Dim());
    for (int i = 0; i < A.Dim(); ++i)
      for (int j = 0; j < i; ++j) {
        z[i * Ncols + j] = A(i, j);
        z[j * Ncols + i] = A(j, i);
      }
    return *this;
  }

  template<typename REAL1>
  CMatrixT &operator+=(const CMatrixT<REAL1> &b) {
    if (Nrows != b.rows() || Ncols != b.cols()) THROW("Size does not fit!")
    for (int i = 0; i < z.size(); ++i)
      z[i] += b.Data()[i];
    return *this;
  }

  template<typename REAL1>
  CMatrixT &operator+=(const HermCMatrixT<REAL1> &b) {
    if (Nrows != b.Dim() || Ncols != b.Dim()) THROW("Size does not fit!")
    for (int i = 0; i < Nrows; ++i)
      for (int j = 0; j < Ncols; ++j)
        z[i * Ncols + j] += b(i, j);
    return *this;
  }

  template<typename REAL1>
  CMatrixT &operator+=(const RMatrixT<REAL1> &b) {
    if (Nrows != b.rows() || Ncols != b.cols()) THROW("Size does not fit!")
    for (int i = 0; i < z.size(); ++i)
      z[i] += b.Data()[i];
    return *this;
  }

  template<typename REAL1>
  CMatrixT &operator+=(const SymRMatrixT<REAL1> &b) {
    if (Nrows != b.Dim() || Ncols != b.Dim()) THROW("Size does not fit!")
    for (int i = 0; i < Nrows; ++i)
      for (int j = 0; j < Ncols; ++j)
        z[i * Ncols + j] += b(i, j);
    return *this;
  }

  template<typename REAL1>
  CMatrixT &operator+=(const AntisymRMatrixT<REAL1> &b) {
    if (Nrows != b.Dim() || Ncols != b.Dim()) THROW("Size does not fit!")
    for (int i = 0; i < Nrows; ++i)
      for (int j = 0; j < Ncols; ++j)
        z[i * Ncols + j] += b(i, j);
    return *this;
  }

  template<typename REAL1>
  CMatrixT &operator-=(const CMatrixT<REAL1> &b) {
    if (Nrows != b.rows() || Ncols != b.cols()) THROW("Size does not fit!")
    for (int i = 0; i < z.size(); ++i)
      z[i] -= b.Data()[i];
    return *this;
  }

  template<typename REAL1>
  CMatrixT &operator-=(const HermCMatrixT<REAL1> &b) {
    if (Nrows != b.Dim() || Ncols != b.Dim()) THROW("Size does not fit!")
    for (int i = 0; i < Nrows; ++i)
      for (int j = 0; j < Ncols; ++j)
        z[i * Ncols + j] -= b(i, j);
    return *this;
  }

  template<typename REAL1>
  CMatrixT &operator-=(const RMatrixT<REAL1> &b) {
    if (Nrows != b.rows() || Ncols != b.cols()) THROW("Size does not fit!")
    for (int i = 0; i < z.size(); ++i)
      z[i] -= b.Data()[i];
    return *this;
  }

  template<typename REAL1>
  CMatrixT &operator-=(const SymRMatrixT<REAL1> &b) {
    if (Nrows != b.Dim() || Ncols != b.Dim()) THROW("Size does not fit!")
    for (int i = 0; i < Nrows; ++i)
      for (int j = 0; j < Ncols; ++j)
        z[i * Ncols + j] -= b(i, j);
    return *this;
  }

  template<typename REAL1>
  CMatrixT &operator-=(const AntisymRMatrixT<REAL1> &b) {
    if (Nrows != b.Dim() || Ncols != b.Dim()) THROW("Size does not fit!")
    for (int i = 0; i < Nrows; ++i)
      for (int j = 0; j < Ncols; ++j)
        z[i * Ncols + j] -= b(i, j);
    return *this;
  }

  CMatrixT &operator*=(const COMPLEX &b) {
    for (int i = 0; i < z.size(); ++i)
      z[i] *= b;
    return *this;
  }

  // This is a component wise operation - not a matrix multiplication!
  template<typename REAL1>
  CMatrixT &operator*=(const CMatrixT<REAL1> &b) {
    if (Nrows != b.rows() || Ncols != b.cols()) THROW("Size does not fit!")
    for (int i = 0; i < z.size(); ++i)
      z[i] *= b.Data()[i];
    return *this;
  }

  CMatrixT &operator/=(const COMPLEX &b) {
    for (int i = 0; i < z.size(); ++i)
      z[i] /= b;
    return *this;
  }

  auto operator[](int row) { return begin() + row * Ncols; }

  auto operator[](int row) const { return begin() + row * Ncols; }

  COMPLEX &operator()(int row, int col) { return z[row * Ncols + col]; }

  COMPLEX &operator()(int index) { return z[index]; }

  const COMPLEX &operator()(int row, int col) const { return z[row * Ncols + col]; }

  template<typename COMPLEX1>
  void operator()(const COMPLEX1 &b, int row, int col) {
    z[row * Ncols + col] = b;
  }

  CVectorT<REAL> row(int i) const {
    CVectorT<REAL> r(Ncols);
    for (int j = 0; j < Ncols; ++j)
      r[j] = z[i * Ncols + j];
    return r;
  }

  CVectorT<REAL> col(int j) const {
    CVectorT<REAL> c(Nrows);
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

  const std::vector<COMPLEX> &Data() const { return z; }

  std::vector<COMPLEX> &Data() { return z; }

  void resize(int rows, int cols);

  void resize(int dim) { resize(dim, dim); }

  CVectorT<REAL> diag();

  template<typename REAL1>
  void diag(const CVectorT<REAL1> &d) {
    resize(d.size());
    for (int i = 0; i < d.size(); ++i)
      z[i * Ncols + i] = d[i];
  }

  template<typename REAL1>
  void diag(const RVectorT<REAL1> &d) {
    resize(d.size());
    for (int i = 0; i < d.size(); ++i)
      z[i * Ncols + i] = d[i];
  }

  void Accumulate(int commSpit = 0);

  CMatrixT &conj();

  RMatrixT<REAL> real() const;

  RMatrixT<REAL> imag() const;

  CMatrixT &transpose();

  COMPLEX Mean() const;

  REAL Variance() const;

  template<typename REAL1>
  void Insert(const CMatrixT<REAL1> &B, int row, int col) {
    if (row + B.rows() >= Nrows || col + B.cols() >= Ncols) {
      THROW("Cannot insert matrix: size of B too large")
    }
    for (int i = 0; i < B.rows(); ++i)
      for (int j = 0; j < B.cols(); ++j)
        z[(row + i) * Ncols + col + j] = B[i][j];
  }

  template<typename REAL1>
  void InsertRow(const CVectorT<REAL1> &v, int row, int startCol = 0) {
    if (row >= Nrows || startCol + v.size() > Ncols) {
      THROW("Cannot insert row: size of v too large")
    }
    for (int j = 0; j < v.size(); ++j)
      z[row * Ncols + startCol + j] = v[j];
  }

  template<typename REAL1>
  void InsertCol(const CVectorT<REAL1> &v, int col, int startRow = 0) {
    if (col >= Ncols || startRow + v.size() > Nrows) {
      THROW("Cannot insert column: size of v too large")
    }
    for (int i = 0; i < v.size(); ++i)
      z[(startRow + i) * Ncols + col] = v[i];
  }

  CMatrixT &adjoint();

  CMatrixT &Identity();

  CMatrixT &Invert();

  Saver &save(Saver &saver) const;

  Loader &load(Loader &loader);
};

template<typename REAL>
bool operator==(const CMatrixT<REAL> &A, const CMatrixT<REAL> &B) {
  if (A.rows() != B.rows() || A.cols() != B.cols()) return false;
  return (A.real() == B.real()) && (A.imag() == B.imag());
}

template<typename REAL>
bool operator==(const CMatrixT<REAL> &A, const HermCMatrixT<REAL> &B) {
  return A == CMatrixT<REAL>(B);
}

template<typename REAL>
bool operator==(const HermCMatrixT<REAL> &A, const CMatrixT<REAL> &B) {
  return CMatrixT<REAL>(A) == B;
}

template<typename REAL>
bool operator==(const CMatrixT<REAL> &A, const RMatrixT<REAL> &B) {
  return A == CMatrixT<REAL>(B);
}

template<typename REAL>
bool operator==(const RMatrixT<REAL> &A, const CMatrixT<REAL> &B) {
  return CMatrixT<REAL>(A) == B;
}

template<typename REAL>
bool operator==(const CMatrixT<REAL> &A, const SymRMatrixT<REAL> &B) {
  return A == CMatrixT<REAL>(B);
}

template<typename REAL>
bool operator==(const SymRMatrixT<REAL> &A, const CMatrixT<REAL> &B) {
  return CMatrixT<REAL>(A) == B;
}

template<typename REAL>
bool operator==(const CMatrixT<REAL> &A, const AntisymRMatrixT<REAL> &B) {
  return A == CMatrixT<REAL>(B);
}

template<typename REAL>
bool operator==(const AntisymRMatrixT<REAL> &A, const CMatrixT<REAL> &B) {
  return CMatrixT<REAL>(A) == B;
}

template<typename REAL>
bool operator!=(const CMatrixT<REAL> &A, const CMatrixT<REAL> &B) {
  return !(A == B);
}

template<typename REAL>
bool operator!=(const CMatrixT<REAL> &A, const HermCMatrixT<REAL> &B) {
  return !(A == B);
}

template<typename REAL>
bool operator!=(const HermCMatrixT<REAL> &A, const CMatrixT<REAL> &B) {
  return !(A == B);
}

template<typename REAL>
bool operator!=(const CMatrixT<REAL> &A, const RMatrixT<REAL> &B) {
  return !(A == B);
}

template<typename REAL>
bool operator!=(const RMatrixT<REAL> &A, const CMatrixT<REAL> &B) {
  return !(A == B);
}

template<typename REAL>
bool operator!=(const CMatrixT<REAL> &A, const SymRMatrixT<REAL> &B) {
  return !(A == B);
}

template<typename REAL>
bool operator!=(const SymRMatrixT<REAL> &A, const CMatrixT<REAL> &B) {
  return !(A == B);
}

template<typename REAL>
bool operator!=(const CMatrixT<REAL> &A, const AntisymRMatrixT<REAL> &B) {
  return !(A == B);
}

template<typename REAL>
bool operator!=(const AntisymRMatrixT<REAL> &A, const CMatrixT<REAL> &B) {
  return !(A == B);
}

template<typename REAL>
inline CMatrixT<REAL> operator-(const CMatrixT<REAL> &a) {
  CMatrixT<REAL> c(a);
  return c *= -1;
}

template<typename REAL>
inline CMatrixT<REAL> operator+(const CMatrixT<REAL> &a, const CMatrixT<REAL> &b) {
  CMatrixT<REAL> c(a);
  return c += b;
}

template<typename REAL>
inline CMatrixT<REAL> operator+(const CMatrixT<REAL> &a, const RMatrixT<REAL> &b) {
  CMatrixT<REAL> c(a);
  return c += b;
}

template<typename REAL>
inline CMatrixT<REAL> operator+(const RMatrixT<REAL> &a, const CMatrixT<REAL> &b) {
  CMatrixT<REAL> c(a);
  return c += b;
}

template<typename REAL>
inline CMatrixT<REAL> operator+(const CMatrixT<REAL> &a, const HermCMatrixT<REAL> &b) {
  CMatrixT<REAL> c(a);
  return c += b;
}

template<typename REAL>
inline CMatrixT<REAL> operator+(const HermCMatrixT<REAL> &a, const CMatrixT<REAL> &b) {
  CMatrixT<REAL> c(a);
  return c += b;
}

template<typename REAL>
inline CMatrixT<REAL> operator+(const CMatrixT<REAL> &a, const SymRMatrixT<REAL> &b) {
  CMatrixT<REAL> c(a);
  return c += b;
}

template<typename REAL>
inline CMatrixT<REAL> operator+(const SymRMatrixT<REAL> &a, const CMatrixT<REAL> &b) {
  CMatrixT<REAL> c(a);
  return c += b;
}

template<typename REAL>
inline CMatrixT<REAL> operator+(const CMatrixT<REAL> &a, const AntisymRMatrixT<REAL> &b) {
  CMatrixT<REAL> c(a);
  return c += b;
}

template<typename REAL>
inline CMatrixT<REAL> operator+(const AntisymRMatrixT<REAL> &a, const CMatrixT<REAL> &b) {
  CMatrixT<REAL> c(a);
  return c += b;
}

template<typename REAL>
inline CMatrixT<REAL> operator-(const CMatrixT<REAL> &a, const CMatrixT<REAL> &b) {
  CMatrixT<REAL> c(a);
  return c -= b;
}

template<typename REAL>
inline CMatrixT<REAL> operator-(const CMatrixT<REAL> &a, const RMatrixT<REAL> &b) {
  CMatrixT<REAL> c(a);
  return c -= b;
}

template<typename REAL>
inline CMatrixT<REAL> operator-(const RMatrixT<REAL> &a, const CMatrixT<REAL> &b) {
  CMatrixT<REAL> c(a);
  return c -= b;
}

template<typename REAL>
inline CMatrixT<REAL> operator-(const CMatrixT<REAL> &a, const HermCMatrixT<REAL> &b) {
  CMatrixT<REAL> c(a);
  return c -= b;
}

template<typename REAL>
inline CMatrixT<REAL> operator-(const HermCMatrixT<REAL> &a, const CMatrixT<REAL> &b) {
  CMatrixT<REAL> c(a);
  return c -= b;
}

template<typename REAL>
inline CMatrixT<REAL> operator-(const CMatrixT<REAL> &a, const SymRMatrixT<REAL> &b) {
  CMatrixT<REAL> c(a);
  return c -= b;
}

template<typename REAL>
inline CMatrixT<REAL> operator-(const SymRMatrixT<REAL> &a, const CMatrixT<REAL> &b) {
  CMatrixT<REAL> c(a);
  return c -= b;
}

template<typename REAL>
inline CMatrixT<REAL> operator-(const CMatrixT<REAL> &a, const AntisymRMatrixT<REAL> &b) {
  CMatrixT<REAL> c(a);
  return c -= b;
}

template<typename REAL>
inline CMatrixT<REAL> operator-(const AntisymRMatrixT<REAL> &a, const CMatrixT<REAL> &b) {
  CMatrixT<REAL> c(a);
  return c -= b;
}

template<typename REAL>
inline CMatrixT<REAL> operator*(const CMatrixT<REAL> &a, const CMatrixT<REAL> &b) {
  return CMatrixT<REAL>(a, b);
}

template<typename REAL>
inline CMatrixT<REAL> operator*(const CMatrixT<REAL> &a, const RMatrixT<REAL> &b) {
  return CMatrixT<REAL>(a, b);
}

template<typename REAL>
inline CMatrixT<REAL> operator*(const RMatrixT<REAL> &a, const CMatrixT<REAL> &b) {
  return CMatrixT<REAL>(a, b);
}

template<typename REAL>
inline CMatrixT<REAL> operator*(const HermCMatrixT<REAL> &a, const HermCMatrixT<REAL> &b) {
  return CMatrixT<REAL>(a, b);
}

template<typename REAL>
inline CMatrixT<REAL> operator*(const CMatrixT<REAL> &a, const HermCMatrixT<REAL> &b) {
  return CMatrixT<REAL>(a, b);
}

template<typename REAL>
inline CMatrixT<REAL> operator*(const HermCMatrixT<REAL> &a, const CMatrixT<REAL> &b) {
  return CMatrixT<REAL>(a, b);
}

template<typename REAL>
inline CMatrixT<REAL> operator*(const CMatrixT<REAL> &a, const SymRMatrixT<REAL> &b) {
  return CMatrixT<REAL>(a, b);
}

template<typename REAL>
inline CMatrixT<REAL> operator*(const SymRMatrixT<REAL> &a, const CMatrixT<REAL> &b) {
  return CMatrixT<REAL>(a, b);
}

template<typename REAL>
inline CMatrixT<REAL> operator*(const CMatrixT<REAL> &a, const AntisymRMatrixT<REAL> &b) {
  return CMatrixT<REAL>(a, b);
}

template<typename REAL>
inline CMatrixT<REAL> operator*(const AntisymRMatrixT<REAL> &a, const CMatrixT<REAL> &b) {
  return CMatrixT<REAL>(a, b);
}

template<typename REAL>
CVectorT<REAL> operator*(const CMatrixT<REAL> &a, const CVectorT<REAL> &v) {
  CVectorT<REAL> w(a.rows());
  for (int i = 0; i < w.size(); ++i)
    for (int k = 0; k < a.cols(); ++k)
      w[i] += a[i][k] * v[k];
  return w;
}

template<typename REAL>
inline CVectorT<REAL> operator*(const CMatrixT<REAL> &a, const RVectorT<REAL> &v) {
  CVectorT<REAL> w(a.rows());
  for (int i = 0; i < w.size(); ++i)
    for (int k = 0; k < a.cols(); ++k)
      w[i] += a[i][k] * v[k];
  return w;
}

template<typename REAL, typename T>
inline CMatrixT<REAL> operator*(const CMatrixT<REAL> &a, const T &b) {
  CMatrixT<REAL> c(a);
  return c *= b;
}

template<typename REAL, typename T>
inline CMatrixT<REAL> operator*(const T &b, const CMatrixT<REAL> &a) {
  CMatrixT<REAL> c(a);
  return c *= b;
}

template<typename REAL, typename T>
inline CMatrixT<REAL> operator/(const CMatrixT<REAL> &a, const T &b) {
  CMatrixT<REAL> c(a);
  return c /= b;
}

template<typename REAL>
inline CMatrixT<REAL> conj(const CMatrixT<REAL> &a) {
  CMatrixT<REAL> b(a);
  return b.conj();
}

template<typename REAL>
inline RMatrixT<REAL> real(const CMatrixT<REAL> &A) {
  return A.real();
}

template<typename REAL>
inline RMatrixT<REAL> imag(const CMatrixT<REAL> &v) {
  return v.imag();
}

template<typename REAL>
inline CMatrixT<REAL> transpose(const CMatrixT<REAL> &a) {
  CMatrixT<REAL> b(a);
  return b.transpose();
}

template<typename REAL>
inline CMatrixT<REAL> invert(const CMatrixT<REAL> &a) {
  CMatrixT<REAL> b(a);
  return b.Invert();
}

template<typename REAL>
inline CMatrixT<REAL> adjoint(const CMatrixT<REAL> &a) {
  CMatrixT<REAL> b(a);
  return b.adjoint();
}

template<typename REAL>
inline Saver &operator<<(Saver &saver, const CMatrixT<REAL> &A) {
  return A.save(saver);
}

template<typename REAL>
inline Loader &operator>>(Loader &loader, CMatrixT<REAL> &A) {
  return A.load(loader);
}

template<typename REAL>
std::ostream &operator<<(std::ostream &os, const CMatrixT<REAL> &v) {
  for (int i = 0; i < v.rows(); ++i) {
    if (i != 0) os << "\n";
    for (int j = 0; j < v.cols(); ++j)
      os << v(i, j) << " ";
  }
  return os;
}

using CMatrix = CMatrixT<double>;

#ifdef BUILD_IA

using IACMatrix = CMatrixT<IAInterval>;

CMatrix mid(const IACMatrix &);

inline IACMatrix operator+(const IACMatrix &A, const CMatrix &B) {
  IACMatrix C(A);
  return C += B;
}

inline IACMatrix operator+(const CMatrix &A, const IACMatrix &B) {
  IACMatrix C(A);
  return C += B;
}

inline IACMatrix operator+(const IACMatrix &A, const HermCMatrix &B) {
  IACMatrix C(A);
  return C += B;
}

inline IACMatrix operator+(const HermCMatrix &A, const IACMatrix &B) {
  IACMatrix C(A);
  return C += B;
}

inline IACMatrix operator+(const IACMatrix &A, const RMatrix &B) {
  IACMatrix C(A);
  return C += B;
}

inline IACMatrix operator+(const RMatrix &A, const IACMatrix &B) {
  IACMatrix C(A);
  return C += B;
}

inline IACMatrix operator+(const IACMatrix &A, const SymRMatrix &B) {
  IACMatrix C(A);
  return C += B;
}

inline IACMatrix operator+(const SymRMatrix &A, const IACMatrix &B) {
  IACMatrix C(A);
  return C += B;
}

inline IACMatrix operator+(const IACMatrix &A, const AntisymRMatrix &B) {
  IACMatrix C(A);
  return C += B;
}

inline IACMatrix operator+(const AntisymRMatrix &A, const IACMatrix &B) {
  IACMatrix C(A);
  return C += B;
}

inline IACMatrix operator-(const IACMatrix &A, const CMatrix &B) {
  IACMatrix C(A);
  return C -= B;
}

inline IACMatrix operator-(const CMatrix &A, const IACMatrix &B) {
  IACMatrix C(A);
  return C -= B;
}

inline IACMatrix operator-(const IACMatrix &A, const HermCMatrix &B) {
  IACMatrix C(A);
  return C -= B;
}

inline IACMatrix operator-(const HermCMatrix &A, const IACMatrix &B) {
  IACMatrix C(A);
  return C -= B;
}

inline IACMatrix operator-(const IACMatrix &A, const RMatrix &B) {
  IACMatrix C(A);
  return C -= B;
}

inline IACMatrix operator-(const RMatrix &A, const IACMatrix &B) {
  IACMatrix C(A);
  return C -= B;
}

inline IACMatrix operator-(const IACMatrix &A, const SymRMatrix &B) {
  IACMatrix C(A);
  return C -= B;
}

inline IACMatrix operator-(const SymRMatrix &A, const IACMatrix &B) {
  IACMatrix C(A);
  return C -= B;
}

inline IACMatrix operator-(const IACMatrix &A, const AntisymRMatrix &B) {
  IACMatrix C(A);
  return C -= B;
}

inline IACMatrix operator-(const AntisymRMatrix &A, const IACMatrix &B) {
  IACMatrix C(A);
  return C -= B;
}

inline IACMatrix operator*(const IACMatrix &A, const CMatrix &B) { return IACMatrix(A, B); }

inline IACMatrix operator*(const CMatrix &A, const IACMatrix &B) { return IACMatrix(A, B); }

inline IACMatrix operator*(const IACMatrix &A, const RMatrix &B) { return IACMatrix(A, B); }

inline IACMatrix operator*(const RMatrix &A, const IACMatrix &B) { return IACMatrix(A, B); }

inline IACMatrix operator*(const IACMatrix &A, const HermCMatrix &B) { return IACMatrix(A, B); }

inline IACMatrix operator*(const HermCMatrix &A, const IACMatrix &B) { return IACMatrix(A, B); }

inline IACMatrix operator*(const IACMatrix &A, const SymRMatrix &B) { return IACMatrix(A, B); }

inline IACMatrix operator*(const SymRMatrix &A, const IACMatrix &B) { return IACMatrix(A, B); }

inline IACMatrix operator*(const IACMatrix &A, const AntisymRMatrix &B) { return IACMatrix(A, B); }

inline IACMatrix operator*(const AntisymRMatrix &A, const IACMatrix &B) { return IACMatrix(A, B); }

inline IACVector operator*(const IACMatrix &a, const CVector &v) {
  IACVector w(a.rows());
  for (int i = 0; i < w.size(); ++i)
    for (int k = 0; k < a.cols(); ++k)
      w[i] += a[i][k] * v[k];
  return w;
}

inline IACVector operator*(const CMatrix &a, const IACVector &v) {
  IACVector w(a.rows());
  for (int i = 0; i < w.size(); ++i)
    for (int k = 0; k < a.cols(); ++k)
      w[i] += a[i][k] * v[k];
  return w;
}

inline IACVector operator*(const IACMatrix &a, const RVector &v) {
  IACVector w(a.rows());
  for (int i = 0; i < w.size(); ++i)
    for (int k = 0; k < a.cols(); ++k)
      w[i] += a[i][k] * v[k];
  return w;
}

inline IACVector operator*(const CMatrix &a, const IARVector &v) {
  IACVector w(a.rows());
  for (int i = 0; i < w.size(); ++i)
    for (int k = 0; k < a.cols(); ++k)
      w[i] += IACInterval(a[i][k]) * v[k];
  return w;
}

#endif // BUILD_IA

#endif
