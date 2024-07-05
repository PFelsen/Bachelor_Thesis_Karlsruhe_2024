#ifndef SYMRMATRIX_H
#define SYMRMATRIX_H

#include "RVector.hpp"

template<typename REAL = double>
class SymRMatrixT {
  std::vector<REAL> z;
  int dim;
public:
  SymRMatrixT() : dim(0){};

  explicit SymRMatrixT(int dim);

  template<typename REAL1>
  constexpr SymRMatrixT(const REAL1 &a, int dim) : dim(dim), z((dim * (dim + 1)) / 2) {
    for (int i = 0; i < z.size(); ++i)
      z[i] = a;
  }

  template<typename REAL1>
  explicit SymRMatrixT(const SymRMatrixT<REAL1> &A) : dim(A.Dim()), z(A.size()) {
    for (int i = 0; i < z.size(); ++i)
      z[i] = A[i];
  }

  template<typename REAL1>
  SymRMatrixT &operator=(const REAL1 &a) {
    for (int i = 0; i < z.size(); ++i)
      z[i] = a;
    return *this;
  }

  template<typename REAL1>
  SymRMatrixT &operator=(const SymRMatrixT<REAL1> &A) {
    dim = A.Dim();
    z.resize(A.size());
    for (int i = 0; i < z.size(); ++i)
      z[i] = A[i];
    return *this;
  }

  template<typename REAL1>
  SymRMatrixT &operator+=(const SymRMatrixT<REAL1> &A) {
    if (dim != A.Dim()) THROW("Size does not fit!")
    for (int i = 0; i < z.size(); ++i)
      z[i] += A[i];
    return *this;
  }

  template<typename REAL1>
  SymRMatrixT &operator-=(const SymRMatrixT<REAL1> &A) {
    if (dim != A.Dim()) THROW("Size does not fit!")
    for (int i = 0; i < z.size(); ++i)
      z[i] -= A[i];
    return *this;
  }

  template<typename REAL1>
  SymRMatrixT &operator*=(const REAL1 &a) {
    for (int i = 0; i < z.size(); ++i)
      z[i] *= a;
    return *this;
  }

  template<typename REAL1>
  SymRMatrixT &operator/=(const REAL1 &a) {
    for (int i = 0; i < z.size(); ++i)
      z[i] /= a;
    return *this;
  }

  RVectorT<REAL> row(int i) const {
    RVectorT<REAL> r(dim);
    for (int j = 0; j < dim; ++j)
      r[j] = (*this)(i, j);
    return r;
  }

  RVectorT<REAL> col(int j) const {
    RVectorT<REAL> c(dim);
    for (int i = 0; i < dim; ++i)
      c[i] = (*this)(i, j);
    return c;
  }

  int size() const { return z.size(); }

  int Dim() const { return dim; }

  int rows() const { return dim; }

  int cols() const { return dim; }

  REAL &operator()(int i, int j) {
    if (i >= j) return z[(i * (i + 1)) / 2 + j];
    else return z[(j * (j + 1)) / 2 + i];
  }

  const REAL &operator()(int i, int j) const {
    if (i >= j) return z[(i * (i + 1)) / 2 + j];
    else return z[(j * (j + 1)) / 2 + i];
  }

  REAL &operator[](int i) { return z[i]; }

  const REAL &operator[](int i) const { return z[i]; }

  template<typename REAL1>
  void operator()(const REAL1 &a, int i, int j) {
    if (i >= j) z[(i * (i + 1)) / 2 + j] = a;
    else z[(j * (j + 1)) / 2 + i] = a;
  }

  auto begin() const { return z.begin(); }

  auto begin() { return z.begin(); }

  auto end() const { return z.end(); }

  auto end() { return z.end(); }

  template<typename REAL1>
  RVectorT<REAL> multiplyWith(const RVectorT<REAL1> &v) const {
    if (v.Dim() != dim) THROW("Dimension does not fit!")
    RVectorT<REAL> y(dim);
    for (int n = 0; n < dim; ++n)
      for (int k = 0; k < dim; ++k)
        y[n] += (*this)(n, k) * v[k];
    return y;
  }

  void resize(int dim);

  void Accumulate(int commSplit = 0);

  SymRMatrixT<REAL> &Identity();

  SymRMatrixT<REAL> &Invert();

  Saver &save(Saver &saver) const;

  Loader &load(Loader &loader);
};

template<typename REAL>
bool operator==(const SymRMatrixT<REAL> &A, const SymRMatrixT<REAL> &B) {
  if (A.Dim() != B.Dim()) return false;
  for (int i = 0; i < A.size(); ++i)
    if (!mpp_ba::isNear(A[i], B[i])) return false;
  return true;
}

template<typename REAL>
bool operator!=(const SymRMatrixT<REAL> &A, const SymRMatrixT<REAL> &B) {
  return !(A == B);
}

template<typename REAL>
inline SymRMatrixT<REAL> operator+(const SymRMatrixT<REAL> &A, const SymRMatrixT<REAL> &B) {
  SymRMatrixT<REAL> AB(A);
  return AB += B;
}

template<typename REAL>
inline SymRMatrixT<REAL> operator-(const SymRMatrixT<REAL> &A, const SymRMatrixT<REAL> &B) {
  SymRMatrixT<REAL> AB(A);
  return AB -= B;
}

template<typename REAL>
inline SymRMatrixT<REAL> operator-(const SymRMatrixT<REAL> &A) {
  SymRMatrixT<REAL> B(A);
  return B *= -1.0;
}

template<typename REAL, typename REAL1>
inline SymRMatrixT<REAL> operator*(const REAL1 &b, const SymRMatrixT<REAL> &A) {
  SymRMatrixT<REAL> B(A);
  return B *= b;
}

template<typename REAL, typename REAL1>
inline SymRMatrixT<REAL> operator*(const SymRMatrixT<REAL> &A, const REAL1 &b) {
  SymRMatrixT<REAL> B(A);
  return B *= b;
}

template<typename REAL, typename REAL1>
inline SymRMatrixT<REAL> operator/(const SymRMatrixT<REAL> &A, const REAL1 &b) {
  SymRMatrixT<REAL> B(A);
  return B /= b;
}

template<typename REAL>
inline RVectorT<REAL> operator*(const SymRMatrixT<REAL> &A, const RVectorT<REAL> &v) {
  return A.multiplyWith(v);
}

template<typename REAL>
inline SymRMatrixT<REAL> invert(const SymRMatrixT<REAL> &A) {
  SymRMatrixT<REAL> B(A);
  return B.Invert();
}

template<typename REAL>
inline Saver &operator<<(Saver &saver, const SymRMatrixT<REAL> &A) {
  return A.save(saver);
}

template<typename REAL>
inline Loader &operator>>(Loader &loader, SymRMatrixT<REAL> &A) {
  return A.load(loader);
}

template<typename REAL>
std::ostream &operator<<(std::ostream &os, const SymRMatrixT<REAL> &A) {
  if constexpr (std::is_same_v<REAL, double>) os << beginD;
  for (int i = 0; i < A.Dim(); ++i) {
    if (i != 0) {
      os << "\n";
      if constexpr (std::is_same_v<REAL, double>) os << beginD;
    }
    for (int j = 0; j < A.Dim(); ++j)
      os << A(i, j) << " ";
  }
  if constexpr (std::is_same_v<REAL, double>) os << endD;
  return os;
}

typedef SymRMatrixT<> SymRMatrix;

#ifdef BUILD_IA

using IASymRMatrix = SymRMatrixT<IAInterval>;

SymRMatrix mid(const IASymRMatrix &);

SymRMatrix sup(const IASymRMatrix &);

SymRMatrix inf(const IASymRMatrix &);

inline IASymRMatrix operator+(const IASymRMatrix &A, const SymRMatrix &B) {
  IASymRMatrix C(A);
  return C += B;
}

inline IASymRMatrix operator+(const SymRMatrix &A, const IASymRMatrix &B) {
  IASymRMatrix C(A);
  return C += B;
}

inline IASymRMatrix operator-(const IASymRMatrix &A, const SymRMatrix &B) {
  IASymRMatrix C(A);
  return C -= B;
}

inline IASymRMatrix operator-(const SymRMatrix &A, const IASymRMatrix &B) {
  IASymRMatrix C(A);
  return C -= B;
}

inline IARVector operator*(const IASymRMatrix &A, const RVector &v) { return A.multiplyWith(v); }

inline IARVector operator*(const SymRMatrix &A, const IARVector &v) {
  IASymRMatrix B(A);
  return B.multiplyWith(v);
}

#endif // BUILD_IA

#endif
