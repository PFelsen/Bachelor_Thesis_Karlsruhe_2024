#ifndef HERMCMATRIX_H
#define HERMCMATRIX_H

#include "AntisymRMatrix.hpp"
#include "CVector.hpp"
#include "SymRMatrix.hpp"

template<typename REAL = double>
class HermCMatrixT {
  using COMPLEX = COMPLEX_TYPE<REAL>;
  std::vector<COMPLEX> z;
  int dim;
public:
  HermCMatrixT() : dim(0) {}

  explicit HermCMatrixT(int dim);

  template<typename REAL1>
  explicit HermCMatrixT(const SymRMatrixT<REAL1> &A) : dim(A.Dim()), z(A.size()) {
    for (int i = 0; i < z.size(); ++i)
      z[i] = A[i];
  }

  template<typename REAL1>
  explicit HermCMatrixT(const HermCMatrixT<REAL1> &A) : dim(A.Dim()), z(A.size()) {
    for (int i = 0; i < z.size(); ++i)
      z[i] = A[i];
  }

  template<typename REAL1>
  HermCMatrixT &operator=(const SymRMatrixT<REAL1> &A) {
    dim = A.Dim();
    z.resize(A.size());
    for (int i = 0; i < z.size(); ++i)
      z[i] = A[i];
    return *this;
  }

  template<typename REAL1>
  HermCMatrixT &operator=(const HermCMatrixT<REAL1> &A) {
    dim = A.Dim();
    z.resize(A.size());
    for (int i = 0; i < z.size(); ++i)
      z[i] = A[i];
    return *this;
  }

  template<typename REAL1>
  HermCMatrixT &operator+=(const SymRMatrixT<REAL1> &A) {
    if (dim != A.Dim()) THROW("Size does not fit!")
    for (int i = 0; i < z.size(); ++i)
      z[i] += A[i];
    return *this;
  }

  template<typename REAL1>
  HermCMatrixT &operator+=(const HermCMatrixT<REAL1> &A) {
    if (dim != A.Dim()) THROW("Size does not fit!")
    for (int i = 0; i < z.size(); ++i)
      z[i] += A[i];
    return *this;
  }

  template<typename REAL1>
  HermCMatrixT &operator-=(const SymRMatrixT<REAL1> &A) {
    if (dim != A.Dim()) THROW("Size does not fit!")
    for (int i = 0; i < z.size(); ++i)
      z[i] -= A[i];
    return *this;
  }

  template<typename REAL1>
  HermCMatrixT &operator-=(const HermCMatrixT<REAL1> &A) {
    if (dim != A.Dim()) THROW("Size does not fit!")
    for (int i = 0; i < z.size(); ++i)
      z[i] -= A[i];
    return *this;
  }

  template<typename REAL1>
  HermCMatrixT &operator*=(const REAL1 &a) {
    for (int i = 0; i < z.size(); ++i)
      z[i] *= REAL(a);
    return *this;
  }

  template<typename REAL1>
  HermCMatrixT &operator/=(const REAL1 &a) {
    for (int i = 0; i < z.size(); ++i)
      z[i] /= REAL(a);
    return *this;
  }

  template<typename REAL1>
  CVectorT<REAL> multiplyWith(const RVectorT<REAL1> &v) const {
    if (v.Dim() != dim) THROW("Dimension does not fit!")
    CVectorT<REAL> y(dim);
    for (int n = 0; n < dim; ++n)
      for (int k = 0; k < dim; ++k)
        y[n] += (*this)(n, k) * v[k];
    return y;
  }

  template<typename REAL1>
  CVectorT<REAL> multiplyWith(const CVectorT<REAL1> &v) const {
    if (v.Dim() != dim) THROW("Dimension does not fit!")
    CVectorT<REAL> y(dim);
    for (int n = 0; n < dim; ++n)
      for (int k = 0; k < dim; ++k)
        y[n] += (*this)(n, k) * v[k];
    return y;
  }

  CVectorT<REAL> row(int i) const {
    CVectorT<REAL> r(dim);
    for (int j = 0; j < dim; ++j)
      r[j] = (*this)(i, j);
    return r;
  }

  CVectorT<REAL> col(int j) const {
    CVectorT<REAL> c(dim);
    for (int i = 0; i < dim; ++i)
      c[i] = (*this)(i, j);
    return c;
  }

  int size() const { return z.size(); }

  int Dim() const { return dim; }

  int rows() const { return dim; }

  int cols() const { return dim; }

  COMPLEX operator()(int i, int j);

  COMPLEX operator()(int i, int j) const;

  COMPLEX &operator[](int i) { return z[i]; }

  const COMPLEX &operator[](int i) const { return z[i]; }

  void operator()(const COMPLEX &a, int i, int j);

  auto begin() const { return z.begin(); }

  auto begin() { return z.begin(); }

  auto end() const { return z.end(); }

  auto end() { return z.end(); }

  void resize(int dim);

  HermCMatrixT &conj();

  SymRMatrixT<REAL> real() const;

  AntisymRMatrixT<REAL> imag() const;

  void Accumulate(int commSplit = 0);

  HermCMatrixT &transpose();

  const HermCMatrixT &adjoint() const;

  HermCMatrixT &Identity();

  HermCMatrixT &Invert();

  Saver &save(Saver &saver) const;

  Loader &load(Loader &loader);
};

template<typename REAL>
bool operator==(const HermCMatrixT<REAL> &A, const HermCMatrixT<REAL> &B) {
  if (A.Dim() != B.Dim()) return false;
  return (A.real() == B.real()) && (A.imag() == B.imag());
}

template<typename REAL, typename REAL1>
bool operator==(const HermCMatrixT<REAL> &A, const SymRMatrixT<REAL1> &B) {
  return A == HermCMatrixT<REAL>(B);
}

template<typename REAL, typename REAL1>
bool operator==(const SymRMatrixT<REAL1> &A, const HermCMatrixT<REAL> &B) {
  return HermCMatrixT<REAL>(A) == B;
}

template<typename REAL>
bool operator!=(const HermCMatrixT<REAL> &A, const HermCMatrixT<REAL> &B) {
  return !(A == B);
}

template<typename REAL, typename REAL1>
bool operator!=(const HermCMatrixT<REAL> &A, const SymRMatrixT<REAL1> &B) {
  return !(A == B);
}

template<typename REAL, typename REAL1>
bool operator!=(const SymRMatrixT<REAL1> &A, const HermCMatrixT<REAL> &B) {
  return !(A == B);
}

template<typename REAL>
inline HermCMatrixT<REAL> operator+(const HermCMatrixT<REAL> &A, const HermCMatrixT<REAL> &B) {
  HermCMatrixT<REAL> AB(A);
  return AB += B;
}

template<typename REAL>
inline HermCMatrixT<REAL> operator+(const SymRMatrixT<REAL> &A, const HermCMatrixT<REAL> &B) {
  HermCMatrixT<REAL> AB(A);
  return AB += B;
}

template<typename REAL>
inline HermCMatrixT<REAL> operator+(const HermCMatrixT<REAL> &A, const SymRMatrixT<REAL> &B) {
  HermCMatrixT<REAL> AB(A);
  return AB += B;
}

template<typename REAL>
inline HermCMatrixT<REAL> operator-(const HermCMatrixT<REAL> &A, const HermCMatrixT<REAL> &B) {
  HermCMatrixT<REAL> AB(A);
  return AB -= B;
}

template<typename REAL>
inline HermCMatrixT<REAL> operator-(const SymRMatrixT<REAL> &A, const HermCMatrixT<REAL> &B) {
  HermCMatrixT<REAL> AB(A);
  return AB -= B;
}

template<typename REAL>
inline HermCMatrixT<REAL> operator-(const HermCMatrixT<REAL> &A, const SymRMatrixT<REAL> &B) {
  HermCMatrixT<REAL> AB(A);
  return AB -= B;
}

template<typename REAL>
inline HermCMatrixT<REAL> operator-(const HermCMatrixT<REAL> &A) {
  HermCMatrixT<REAL> B(A);
  return B *= -1.0;
}

template<typename REAL, typename REAL1>
inline HermCMatrixT<REAL> operator*(const REAL1 &b, const HermCMatrixT<REAL> &A) {
  HermCMatrixT<REAL> B(A);
  return B *= b;
}

template<typename REAL, typename REAL1>
inline HermCMatrixT<REAL> operator*(const HermCMatrixT<REAL> &A, const REAL1 &b) {
  HermCMatrixT<REAL> B(A);
  return B *= b;
}

template<typename REAL, typename REAL1>
inline HermCMatrixT<REAL> operator/(const HermCMatrixT<REAL> &A, const REAL1 &b) {
  HermCMatrixT<REAL> B(A);
  return B /= b;
}

template<typename REAL>
inline CVectorT<REAL> operator*(const HermCMatrixT<REAL> &A, const RVectorT<REAL> &v) {
  return A.multiplyWith(v);
}

template<typename REAL>
inline CVectorT<REAL> operator*(const HermCMatrixT<REAL> &A, const CVectorT<REAL> &v) {
  return A.multiplyWith(v);
}

template<typename REAL>
inline SymRMatrixT<REAL> real(const HermCMatrixT<REAL> &A) {
  return A.real();
}

template<typename REAL>
inline AntisymRMatrixT<REAL> imag(const HermCMatrixT<REAL> &v) {
  return v.imag();
}

template<typename REAL>
inline HermCMatrixT<REAL> conj(const HermCMatrixT<REAL> &A) {
  HermCMatrixT<REAL> B(A);
  return B.conj();
}

template<typename REAL>
inline HermCMatrixT<REAL> transpose(const HermCMatrixT<REAL> &A) {
  HermCMatrixT<REAL> B(A);
  return B.transpose();
}

template<typename REAL>
inline HermCMatrixT<REAL> invert(const HermCMatrixT<REAL> &A) {
  HermCMatrixT<REAL> B(A);
  return B.Invert();
}

template<typename REAL>
inline const HermCMatrixT<REAL> &adjoint(const HermCMatrixT<REAL> &A) {
  return A.adjoint();
}

template<typename REAL>
inline Saver &operator<<(Saver &saver, const HermCMatrixT<REAL> &A) {
  return A.save(saver);
}

template<typename REAL>
inline Loader &operator>>(Loader &loader, HermCMatrixT<REAL> &A) {
  return A.load(loader);
}

template<typename REAL>
std::ostream &operator<<(std::ostream &os, const HermCMatrixT<REAL> &v) {
  for (int i = 0; i < v.Dim(); ++i) {
    if (i != 0) os << "\n";
    for (int j = 0; j < v.Dim(); ++j)
      os << v(i, j) << " ";
  }
  return os;
}

using HermCMatrix = HermCMatrixT<double>;

#ifdef BUILD_IA

using IAHermCMatrix = HermCMatrixT<IAInterval>;

HermCMatrix mid(const IAHermCMatrix &);

inline IAHermCMatrix operator+(const IAHermCMatrix &A, const HermCMatrix &B) {
  IAHermCMatrix C(A);
  return C += B;
}

inline IAHermCMatrix operator+(const HermCMatrix &A, const IAHermCMatrix &B) {
  IAHermCMatrix C(A);
  return C += B;
}

inline IAHermCMatrix operator+(const IAHermCMatrix &A, const SymRMatrix &B) {
  IAHermCMatrix C(A);
  return C += B;
}

inline IAHermCMatrix operator+(const SymRMatrix &A, const IAHermCMatrix &B) {
  IAHermCMatrix C(A);
  return C += B;
}

inline IAHermCMatrix operator-(const IAHermCMatrix &A, const HermCMatrix &B) {
  IAHermCMatrix C(A);
  return C -= B;
}

inline IAHermCMatrix operator-(const HermCMatrix &A, const IAHermCMatrix &B) {
  IAHermCMatrix C(A);
  return C -= B;
}

inline IAHermCMatrix operator-(const IAHermCMatrix &A, const SymRMatrix &B) {
  IAHermCMatrix C(A);
  return C -= B;
}

inline IAHermCMatrix operator-(const SymRMatrix &A, const IAHermCMatrix &B) {
  IAHermCMatrix C(A);
  return C -= B;
}

inline IACVector operator*(const IAHermCMatrix &A, const RVector &v) { return A.multiplyWith(v); }

inline IACVector operator*(const HermCMatrix &A, const IARVector &v) {
  IAHermCMatrix B(A);
  return B.multiplyWith(v);
}

inline IACVector operator*(const IAHermCMatrix &A, const CVector &v) { return A.multiplyWith(v); }

inline IACVector operator*(const HermCMatrix &A, const IACVector &v) {
  IAHermCMatrix B(A);
  return B.multiplyWith(v);
}

#endif // BUILD_IA

#endif
