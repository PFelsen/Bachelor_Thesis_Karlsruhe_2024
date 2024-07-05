#ifndef ANTISYMRMATRIX_H
#define ANTISYMRMATRIX_H

#include "RVector.hpp"

//==================================================================================================
// BLAS/LAPACK routines
//==================================================================================================
// LU factorization
extern "C" void dgetrf_(int *M, int *N, void *A, int *LDA, int *IPIV, int *INFO);
// Invert matrix using LU factorization
extern "C" void dgetri_(int *N, void *A, int *LDA, int *IPIV, void *WORK, int *LWORK, int *INFO);

//==================================================================================================

template<typename REAL = double>
class AntisymRMatrixT {
  std::vector<REAL> z;
  int dim;
public:
  AntisymRMatrixT() : dim(0) {}

  explicit AntisymRMatrixT(int dim);

  template<typename REAL1>
  explicit AntisymRMatrixT(const AntisymRMatrixT<REAL1> &A) : dim(A.Dim()), z(A.size()) {
    for (int i = 0; i < z.size(); ++i)
      z[i] = A[i];
  }

  template<typename REAL1>
  AntisymRMatrixT &operator=(const AntisymRMatrixT<REAL1> &A) {
    dim = A.Dim();
    z.resize(A.size());
    for (int i = 0; i < z.size(); ++i)
      z[i] = A[i];
    return *this;
  }

  template<typename REAL1>
  AntisymRMatrixT &operator+=(const AntisymRMatrixT<REAL1> &A) {
    if (dim != A.Dim()) THROW("Size does not fit!")
    for (int i = 0; i < z.size(); ++i)
      z[i] += A[i];
    return *this;
  }

  template<typename REAL1>
  AntisymRMatrixT &operator-=(const AntisymRMatrixT<REAL1> &A) {
    if (dim != A.Dim()) THROW("Size does not fit!")
    for (int i = 0; i < z.size(); ++i)
      z[i] -= A[i];
    return *this;
  }

  template<typename REAL1>
  AntisymRMatrixT &operator*=(const REAL1 &a) {
    for (int i = 0; i < z.size(); ++i)
      z[i] *= a;
    return *this;
  }

  template<typename REAL1>
  AntisymRMatrixT &operator/=(const REAL1 &a) {
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

  REAL operator()(int i, int j);

  REAL operator()(int i, int j) const;

  REAL &operator[](int i) { return z[i]; }

  const REAL &operator[](int i) const { return z[i]; }

  template<typename REAL1>
  void operator()(const REAL1 &b, int i, int j) {
    if (i == j) return;
    if (i > j) z[(i * (i - 1)) / 2 + j] = b;
    else z[(j * (j - 1)) / 2 + i] = -b;
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
      for (int k = 0; k < dim; ++k) {
        if (n == k) continue;
        y[n] += (*this)(n, k) * v[k];
      }
    return y;
  }

  void resize(int dim);

  AntisymRMatrixT &transpose();

  AntisymRMatrixT &Invert();

  void Accumulate(int commSplit = 0);

  Saver &save(Saver &saver) const;

  Loader &load(Loader &loader);
};

template<typename REAL>
bool operator==(const AntisymRMatrixT<REAL> &A, const AntisymRMatrixT<REAL> &B) {
  if (A.Dim() != B.Dim()) return false;
  for (int i = 0; i < A.size(); ++i)
    if (!mpp_ba::isNear(A[i], B[i])) return false;
  return true;
}

template<typename REAL>
bool operator!=(const AntisymRMatrixT<REAL> &A, const AntisymRMatrixT<REAL> &B) {
  return !(A == B);
}

template<typename REAL>
inline AntisymRMatrixT<REAL> operator+(const AntisymRMatrixT<REAL> &A,
                                       const AntisymRMatrixT<REAL> &B) {
  AntisymRMatrixT<REAL> AB(A);
  return AB += B;
}

template<typename REAL>
inline AntisymRMatrixT<REAL> operator-(const AntisymRMatrixT<REAL> &A,
                                       const AntisymRMatrixT<REAL> &B) {
  AntisymRMatrixT<REAL> AB(A);
  return AB -= B;
}

template<typename REAL>
inline AntisymRMatrixT<REAL> operator-(const AntisymRMatrixT<REAL> &A) {
  AntisymRMatrixT<REAL> B(A);
  return B *= -1.0;
}

template<typename REAL, typename REAL1>
inline AntisymRMatrixT<REAL> operator*(const REAL1 &b, const AntisymRMatrixT<REAL> &A) {
  AntisymRMatrixT<REAL> B(A);
  return B *= b;
}

template<typename REAL, typename REAL1>
inline AntisymRMatrixT<REAL> operator*(const AntisymRMatrixT<REAL> &A, const REAL1 &b) {
  AntisymRMatrixT<REAL> B(A);
  return B *= b;
}

template<typename REAL, typename REAL1>
inline AntisymRMatrixT<REAL> operator/(const AntisymRMatrixT<REAL> &A, const REAL1 &b) {
  AntisymRMatrixT<REAL> B(A);
  return B /= b;
}

template<typename REAL>
inline RVectorT<REAL> operator*(const AntisymRMatrixT<REAL> &A, const RVectorT<REAL> &v) {
  return A.multiplyWith(v);
}

template<typename REAL>
inline AntisymRMatrixT<REAL> transpose(const AntisymRMatrixT<REAL> &A) {
  AntisymRMatrixT<REAL> B(A);
  return B.transpose();
}

template<typename REAL>
inline AntisymRMatrixT<REAL> invert(const AntisymRMatrixT<REAL> &a) {
  AntisymRMatrixT<REAL> b(a);
  return b.Invert();
}

template<typename REAL>
inline Saver &operator<<(Saver &saver, const AntisymRMatrixT<REAL> &A) {
  return A.save(saver);
}

template<typename REAL>
inline Loader &operator>>(Loader &loader, AntisymRMatrixT<REAL> &A) {
  return A.load(loader);
}

template<typename REAL>
std::ostream &operator<<(std::ostream &os, const AntisymRMatrixT<REAL> &A) {
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

typedef AntisymRMatrixT<> AntisymRMatrix;

#ifdef BUILD_IA

using IAAntisymRMatrix = AntisymRMatrixT<IAInterval>;

template<typename REAL>
class RMatrixT;

AntisymRMatrix mid(const IAAntisymRMatrix &);

// RMatrix is needed since upper and lower bound do not have to coincide
RMatrixT<double> sup(const IAAntisymRMatrix &);

// RMatrix is needed since upper and lower bound do not have to coincide
RMatrixT<double> inf(const IAAntisymRMatrix &);

inline IAAntisymRMatrix operator+(const IAAntisymRMatrix &A, const AntisymRMatrix &B) {
  IAAntisymRMatrix C(A);
  return C += B;
}

inline IAAntisymRMatrix operator+(const AntisymRMatrix &A, const IAAntisymRMatrix &B) {
  IAAntisymRMatrix C(A);
  return C += B;
}

inline IAAntisymRMatrix operator-(const IAAntisymRMatrix &A, const AntisymRMatrix &B) {
  IAAntisymRMatrix C(A);
  return C -= B;
}

inline IAAntisymRMatrix operator-(const AntisymRMatrix &A, const IAAntisymRMatrix &B) {
  IAAntisymRMatrix C(A);
  return C -= B;
}

inline IARVector operator*(const IAAntisymRMatrix &A, const RVector &v) {
  return A.multiplyWith(v);
}

inline IARVector operator*(const AntisymRMatrix &A, const IARVector &v) {
  IAAntisymRMatrix B(A);
  return B.multiplyWith(v);
}

#endif // BUILD_IA

#endif
