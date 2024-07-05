#ifndef _SYMTENSOR_H_
#define _SYMTENSOR_H_

#include "VectorField.hpp"

template<typename T = double, int dim = SpaceDimension>
class SymTensorT {
  T z[(dim * (dim + 1)) / 2]{};
public:
  constexpr SymTensorT() {
    for (int n = 0; n < (dim * (dim + 1)) / 2; ++n)
      z[n] = T{};
  }

  constexpr SymTensorT(const T &a) {
    for (int n = 0; n < (dim * (dim + 1)) / 2; ++n)
      z[n] = a;
  }

  template<typename TT>
  constexpr SymTensorT(const SymTensorT<TT, dim> &S) {
    for (int n = 0; n < (dim * (dim + 1)) / 2; ++n)
      z[n] = S[n];
  }

  constexpr SymTensorT(T a00, T a10, T a11, T a20 = T{}, T a21 = T{}, T a22 = T(1.0)) {
    z[0] = a00;
    if constexpr (dim > 1) {
      z[1] = a10;
      z[2] = a11;
    }
    if constexpr (dim > 2) {
      z[3] = a20;
      z[4] = a21;
      z[5] = a22;
    }
  }

  constexpr int Dim() const { return dim; }

  constexpr int size() const { return (dim * (dim + 1)) / 2; }

  constexpr T &operator[](int i) { return z[i]; }

  constexpr const T &operator[](int i) const { return z[i]; }

  constexpr const T &operator()(int i, int j) const {
    if (i >= j) return z[(i * (i + 1)) / 2 + j];
    return z[(j * (j + 1)) / 2 + i];
  }

  constexpr T &operator()(int i, int j) {
    if (i >= j) return z[(i * (i + 1)) / 2 + j];
    return z[(j * (j + 1)) / 2 + i];
  }

  constexpr SymTensorT operator=(const T &a) {
    for (int n = 0; n < (dim * (dim + 1)) / 2; ++n)
      z[n] = a;
    return *this;
  }

  template<typename TT>
  constexpr SymTensorT operator=(const SymTensorT<TT, dim> &S) {
    for (int n = 0; n < (dim * (dim + 1)) / 2; ++n)
      z[n] = S[n];
    return *this;
  }

  template<typename TT>
  constexpr SymTensorT operator+=(const SymTensorT<TT, dim> &S) {
    for (int n = 0; n < (dim * (dim + 1)) / 2; ++n)
      z[n] += S[n];
    return *this;
  }

  template<typename TT>
  constexpr SymTensorT operator-=(const SymTensorT<TT, dim> &S) {
    for (int n = 0; n < (dim * (dim + 1)) / 2; ++n)
      z[n] -= S[n];
    return *this;
  }

  template<typename TT>
  constexpr SymTensorT operator*=(const TT &a) {
    for (int n = 0; n < (dim * (dim + 1)) / 2; ++n)
      z[n] *= a;
    return *this;
  }

  template<typename TT>
  constexpr SymTensorT operator/=(const TT &a) {
    for (int n = 0; n < (dim * (dim + 1)) / 2; ++n)
      z[n] /= a;
    return *this;
  }

  template<typename TT>
  constexpr VectorFieldT<T, dim> multiplyWith(const VectorFieldT<TT, dim> &v) const {
    VectorFieldT<T, dim> y;
    for (int n = 0; n < dim; ++n)
      for (int k = 0; k < dim; ++k)
        y[n] += (*this)(n, k) * v[k];
    return y;
  }

  constexpr T det() const {
    if constexpr (dim == 1) {
      return z[0];
    } else if constexpr (dim == 2) {
      return z[0] * z[2] - z[1] * z[1];
    } else if constexpr (dim == 3) {
      return z[0] * (z[2] * z[5] - z[4] * z[4]) + z[1] * (z[4] * z[3] - z[1] * z[5])
             + z[3] * (z[1] * z[4] - z[2] * z[3]);
    } else THROW("Determinant not defined for dimension greater than 3!")
  }

  constexpr T trace() const {
    T sum{};
    int cnt = 0;
    for (int i = 0; i < dim; ++i) {
      sum += z[cnt];
      cnt += i + 2;
    }
    return sum;
  }

  friend constexpr bool operator==(const SymTensorT<T, dim> &R, const SymTensorT<T, dim> &S) {
    for (int n = 0; n < (dim * (dim + 1)) / 2; ++n)
      if (!mpp_geometry::isNear(R[n], S[n])) return false;
    return true;
  }

  friend constexpr bool operator!=(const SymTensorT<T, dim> &R, const SymTensorT<T, dim> &S) {
    return !(R == S);
  }

  template<typename OSTREAM>
  constexpr OSTREAM &print(OSTREAM &os) const {
    if constexpr (std::is_same<OSTREAM, TextLogger>::value) os << beginD;
    os << "(" << (*this)(0, 0);
    for (int i = 1; i < dim; ++i)
      os << ", " << (*this)(0, i);
    os << ")";
    for (int i = 1; i < dim; ++i) {
      os << "\n(" << (*this)(i, 0);
      for (int j = 1; j < dim; ++j)
        os << ", " << (*this)(i, j);
      os << ")";
    }
    if constexpr (std::is_same<OSTREAM, TextLogger>::value) os << endD;
    return os;
  }

  Saver &save(Saver &saver) const {
    saver << dim;
    for (int n = 0; n < (dim * (dim + 1)) / 2; ++n)
      saver << z[n];
    return saver;
  }

  Loader &load(Loader &loader) {
    int loadedDim;
    loader >> loadedDim;
    if (loadedDim != dim) THROW("SymTensor Dimension does not fit")
    for (int n = 0; n < (dim * (dim + 1)) / 2; ++n)
      loader >> z[n];
    return loader;
  }
};

template<typename T, int dim>
constexpr inline SymTensorT<T, dim> operator+(const SymTensorT<T, dim> &R,
                                              const SymTensorT<T, dim> &S) {
  SymTensorT<T, dim> RS(R);
  return RS += S;
}

template<typename T, int dim>
constexpr inline SymTensorT<T, dim> operator-(const SymTensorT<T, dim> &R,
                                              const SymTensorT<T, dim> &S) {
  SymTensorT<T, dim> RS(R);
  return RS -= S;
}

template<typename T, int dim>
constexpr inline SymTensorT<T, dim> operator-(const SymTensorT<T, dim> &S) {
  SymTensorT<T, dim> R(S);
  return R *= T(-1.0);
}

template<typename T, typename TT, int dim>
constexpr inline SymTensorT<T, dim> operator*(const TT &b, const SymTensorT<T, dim> &S) {
  SymTensorT<T, dim> R(S);
  return R *= b;
}

template<typename T, typename TT, int dim>
constexpr inline SymTensorT<T, dim> operator*(const SymTensorT<T, dim> &S, const TT &b) {
  SymTensorT<T, dim> R(S);
  return R *= b;
}

template<typename T, typename TT, int dim>
constexpr inline SymTensorT<T, dim> operator/(const SymTensorT<T, dim> &S, const TT &b) {
  SymTensorT<T, dim> R(S);
  return R /= b;
}

template<typename T, int dim>
constexpr VectorFieldT<T, dim> operator*(const SymTensorT<T, dim> &S,
                                         const VectorFieldT<T, dim> &v) {
  return S.multiplyWith(v);
}

template<typename T, int dim>
constexpr inline T det(const SymTensorT<T, dim> &S) {
  return S.det();
}

template<typename T, int dim>
constexpr inline T trace(const SymTensorT<T, dim> &S) {
  return S.trace();
}

template<typename T, int dim>
constexpr T norm(const SymTensorT<T, dim> &S) {
  T sum{};
  for (int i = 0; i < dim; ++i)
    for (int j = 0; j < dim; ++j)
      sum += S(i, j) * S(i, j);
  return sqrt(sum);
}

/// templates cannot be used for the stream since Saver has the same stucture
template<typename T, int dim>
constexpr inline std::ostream &operator<<(std::ostream &os, const SymTensorT<T, dim> &S) {
  return S.print(os);
}

template<typename T, int dim>
inline Saver &operator<<(Saver &saver, const SymTensorT<T, dim> &S) {
  return S.save(saver);
}

template<typename T, int dim>
inline Loader &operator>>(Loader &loader, SymTensorT<T, dim> &S) {
  return S.load(loader);
}

using SymTensor = SymTensorT<double, SpaceDimension>;

#ifdef BUILD_IA

using IASymTensor = SymTensorT<IAInterval, SpaceDimension>;

template<int sDim>
SymTensorT<double, sDim> mid(const SymTensorT<IAInterval, sDim> &IA_T) {
  SymTensorT<double, sDim> T;
  for (int n = 0; n < (sDim * (sDim + 1) / 2); ++n)
    T[n] = mid(IA_T[n]);
  return T;
}

template<int sDim>
SymTensorT<double, sDim> sup(const SymTensorT<IAInterval, sDim> &IA_T) {
  SymTensorT<double, sDim> T;
  for (int n = 0; n < (sDim * (sDim + 1) / 2); ++n)
    T[n] = sup(IA_T[n]);
  return T;
}

template<int sDim>
SymTensorT<double, sDim> inf(const SymTensorT<IAInterval, sDim> &IA_T) {
  SymTensorT<double, sDim> T;
  for (int n = 0; n < (sDim * (sDim + 1) / 2); ++n)
    T[n] = inf(IA_T[n]);
  return T;
}

template<int sDim>
inline SymTensorT<IAInterval, sDim> operator+(const SymTensorT<double, sDim> &x,
                                              const SymTensorT<IAInterval, sDim> &y) {
  SymTensorT<IAInterval, sDim> IAx(x);
  return IAx += y;
}

template<int sDim>
inline SymTensorT<IAInterval, sDim> operator+(const SymTensorT<IAInterval, sDim> &x,
                                              const SymTensorT<double, sDim> &y) {
  SymTensorT<IAInterval, sDim> IAx(x);
  return IAx += y;
}

template<int sDim>
inline SymTensorT<IAInterval, sDim> operator-(const SymTensorT<double, sDim> &x,
                                              const SymTensorT<IAInterval, sDim> &y) {
  SymTensorT<IAInterval, sDim> IAx(x);
  return IAx -= y;
}

template<int sDim>
inline SymTensorT<IAInterval, sDim> operator-(const SymTensorT<IAInterval, sDim> &x,
                                              const SymTensorT<double, sDim> &y) {
  SymTensorT<IAInterval, sDim> IAx(x);
  return IAx -= y;
}

template<int sDim>
inline SymTensorT<IAInterval, sDim> operator*(const SymTensorT<double, sDim> &x,
                                              const IAInterval &b) {
  SymTensorT<IAInterval, sDim> IAx(x);
  return IAx *= b;
}

template<int sDim>
inline SymTensorT<IAInterval, sDim> operator*(const IAInterval &b,
                                              const SymTensorT<double, sDim> &x) {
  SymTensorT<IAInterval, sDim> IAx(x);
  return IAx *= b;
}

template<int sDim>
inline SymTensorT<IAInterval, sDim> operator/(const SymTensorT<double, sDim> &x,
                                              const IAInterval &b) {
  SymTensorT<IAInterval, sDim> IAx(x);
  return IAx /= b;
}

template<int sDim>
inline VectorFieldT<IAInterval, sDim> operator*(const SymTensorT<double, sDim> &A,
                                                const VectorFieldT<IAInterval, sDim> &v) {
  SymTensorT<IAInterval, sDim> B(A);
  return B.multiplyWith(v);
}

template<int sDim>
inline VectorFieldT<IAInterval, sDim> operator*(const SymTensorT<IAInterval, sDim> &A,
                                                const VectorFieldT<double, sDim> &v) {
  return A.multiplyWith(v);
}

#endif // BUILD_IA

#endif // of #ifndef _SYMTENSOR_H_
