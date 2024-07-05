#ifndef _VECTORFIELD_H_
#define _VECTORFIELD_H_

#include "VectorFieldComponent.hpp"

template<typename T = double, int dim = SpaceDimension, int tDim = TimeDimension>
class VectorFieldT {
  T z[dim + tDim]{};
public:
  constexpr VectorFieldT() {
    for (int n = 0; n < dim + tDim; ++n)
      z[n] = T{};
  }

  constexpr VectorFieldT(const T &a) {
    for (int n = 0; n < dim + tDim; ++n)
      z[n] = a;
  }

  constexpr VectorFieldT(T a, T b, T c = T{}, T d = T{}) {
    if constexpr (dim + tDim > 0) z[0] = a;
    if constexpr (dim + tDim > 1) z[1] = b;
    if constexpr (dim + tDim > 2) z[2] = c;
    if constexpr (dim + tDim > 3) z[3] = d;
  }

  constexpr VectorFieldT(int k, T a) {
    for (int n = 0; n < dim + tDim; ++n)
      z[n] = T{};
    z[k] = a;
  }

  template<typename TT>
  constexpr VectorFieldT(const VectorFieldT<TT, dim> &y) {
    for (int n = 0; n < dim; ++n)
      z[n] = T(y[n]);
  }

  template<int tDim2>
  constexpr VectorFieldT(const PointT<T, dim, tDim2> &y) {
    for (int n = 0; n < dim; ++n)
      z[n] = y[n];
  }

  constexpr explicit VectorFieldT(const VectorFieldComponentT<T> &y) {
    for (int n = 0; n < dim + tDim; ++n)
      z[n] = T{};
    z[y.Component()] = y.Value();
  }

  constexpr VectorFieldT(T *a) {
    for (int n = 0; n < dim + tDim; ++n)
      z[n] = a[n];
  }

  constexpr int Dim() const { return dim; }

  const T &t() const {
    if constexpr (tDim > 0) return z[dim];
    THROW("Do not call t() with tDim = 0.")
  }

  T &t() {
    if constexpr (tDim > 0) return z[dim];
    THROW("Do not call t() with tDim = 0.")
  }

  constexpr T operator[](const int i) const { return z[i]; }

  constexpr T &operator[](const int i) { return z[i]; }

  constexpr T *operator()() const { return z; }

  constexpr VectorFieldT &operator=(const T &a) {
    for (int n = 0; n < dim + tDim; ++n)
      z[n] = a;
    return *this;
  }

  template<typename TT>
  constexpr VectorFieldT &operator=(const VectorFieldT<TT, dim> &y) {
    for (int n = 0; n < dim; ++n)
      z[n] = y[n];
    return *this;
  }

  template<typename TT>
  constexpr VectorFieldT &operator=(const PointT<TT, dim, tDim> &y) {
    for (int n = 0; n < dim + tDim; ++n)
      z[n] = y[n];
    return *this;
  }

  template<typename TT>
  constexpr VectorFieldT &operator+=(const VectorFieldT<TT, dim> &y) {
    for (int n = 0; n < dim; ++n)
      z[n] += y[n];
    return *this;
  }

  template<typename TT>
  constexpr VectorFieldT &operator+=(const VectorFieldComponentT<TT> &y) {
    z[y.Component()] += y.Value();
    return *this;
  }

  template<typename TT, int tDim2>
  constexpr VectorFieldT &operator+=(const PointT<TT, dim, tDim2> &y) {
    for (int n = 0; n < dim; ++n)
      z[n] += y[n];
    return *this;
  }

  template<typename TT>
  constexpr VectorFieldT &operator+=(const TT &a) {
    for (int n = 0; n < dim + tDim; ++n)
      z[n] += a;
    return *this;
  }

  template<typename TT>
  constexpr VectorFieldT &operator-=(const VectorFieldT<TT, dim> &y) {
    for (int n = 0; n < dim; ++n)
      z[n] -= y[n];
    return *this;
  }

  template<typename TT>
  constexpr VectorFieldT &operator-=(const VectorFieldComponentT<TT> &y) {
    z[y.Component()] -= y.Value();
    return *this;
  }

  template<typename TT, int tDim2>
  constexpr VectorFieldT &operator-=(const PointT<TT, dim, tDim2> &y) {
    for (int n = 0; n < dim; ++n)
      z[n] -= y[n];
    return *this;
  }

  template<typename TT>
  constexpr VectorFieldT &operator-=(const TT &a) {
    for (int n = 0; n < dim + tDim; ++n)
      z[n] -= a;
    return *this;
  }

  template<typename TT>
  constexpr VectorFieldT &operator*=(const TT &a) {
    for (int n = 0; n < dim + tDim; ++n)
      z[n] *= a;
    return *this;
  }

  template<typename TT>
  constexpr VectorFieldT &operator/=(const TT &a) {
    for (int n = 0; n < dim + tDim; ++n)
      z[n] /= a;
    return *this;
  }

  constexpr bool operator==(const VectorFieldT<T, dim> &other) const {
    for (int i = 0; i < dim + tDim; ++i)
      if (!mpp_geometry::isNear(z[i], other[i])) return false;
    return true;
  }

  constexpr bool operator!=(const VectorFieldT<T, dim> &other) const { return !((*this) == other); }

  template<typename OSTREAM>
  constexpr OSTREAM &print(OSTREAM &os) const {
    if constexpr (std::is_same<OSTREAM, TextLogger>::value) os << beginD;
    os << "(" << z[0];
    for (int i = 1; i < dim + tDim; ++i)
      os << ", " << z[i];
    os << ")";
    if constexpr (std::is_same<OSTREAM, TextLogger>::value) os << endD;
    return os;
  }

  Saver &save(Saver &saver) const {
    saver << dim + tDim;
    for (int n = 0; n < dim + tDim; ++n)
      saver << z[n];
    return saver;
  }

  Loader &load(Loader &loader) {
    int loadedDim;
    loader >> loadedDim;
    if (loadedDim != dim + tDim) THROW("VectorField Dimension does not fit")
    for (int n = 0; n < dim + tDim; ++n)
      loader >> z[n];
    return loader;
  }

  constexpr VectorFieldT<T, dim> sliceTime() {
    if constexpr (dim == 1) {
      return VectorFieldT<T, dim>(z[0]);
    } else if constexpr (dim == 2) {
      return VectorFieldT<T, dim>(z[0], z[1]);
    } else if constexpr (dim == 3) {
      return VectorFieldT<T, dim>(z[0], z[1], z[2]);
    }
    Exit("Not implemented for dim = " + std::to_string(dim));
  }
};

template<typename T, int dim>
inline constexpr VectorFieldT<T, dim> operator+(const VectorFieldT<T, dim> &x,
                                                const VectorFieldT<T, dim> &y) {
  VectorFieldT<T, dim> z(x);
  return z += y;
}

template<typename T, int dim>
inline constexpr VectorFieldT<T, dim> operator+(const VectorFieldT<T, dim> &x,
                                                const VectorFieldComponentT<T> &y) {
  VectorFieldT<T, dim> z(x);
  return z += y;
}

template<typename T, int dim>
inline constexpr VectorFieldT<T, dim> operator+(const VectorFieldComponentT<T> &x,
                                                const VectorFieldT<T, dim> &y) {
  VectorFieldT<T, dim> z(y);
  return z += x;
}

template<typename T, int dim, int tDim>
inline constexpr VectorFieldT<T, dim> operator+(const VectorFieldT<T, dim> &x,
                                                const PointT<T, dim, tDim> &y) {
  VectorFieldT<T, dim> z(x);
  return z += y;
}

template<typename T, int dim, int tDim>
inline constexpr VectorFieldT<T, dim> operator+(const PointT<T, dim, tDim> &x,
                                                const VectorFieldT<T, dim> &y) {
  VectorFieldT<T, dim> z(x);
  return z += y;
}

template<typename T, int dim>
inline constexpr VectorFieldT<T, dim> operator-(const VectorFieldT<T, dim> &x,
                                                const VectorFieldT<T, dim> &y) {
  VectorFieldT<T, dim> z(x);
  return z -= y;
}

template<typename T, int dim>
inline constexpr VectorFieldT<T, dim> operator-(const VectorFieldT<T, dim> &x,
                                                const VectorFieldComponentT<T> &y) {
  VectorFieldT<T, dim> z(x);
  return z -= y;
}

template<typename T, int dim>
inline constexpr VectorFieldT<T, dim> operator-(const VectorFieldComponentT<T> &x,
                                                const VectorFieldT<T, dim> &y) {
  VectorFieldT<T, dim> z(x);
  return z -= y;
}

template<typename T, int dim, int tDim>
inline constexpr VectorFieldT<T, dim> operator-(const PointT<T, dim, tDim> &x,
                                                const VectorFieldT<T, dim> &y) {
  VectorFieldT<T, dim> z(x);
  return z -= y;
}

template<typename T, int dim, int tDim>
inline constexpr VectorFieldT<T, dim> operator-(const VectorFieldT<T, dim> &x,
                                                const PointT<T, dim, tDim> &y) {
  VectorFieldT<T, dim> z(x);
  return z -= y;
}

template<typename T, int dim>
inline constexpr VectorFieldT<T, dim> operator-(const VectorFieldT<T, dim> &x) {
  VectorFieldT<T, dim> z(x);
  return z *= T(-1.0);
}

template<typename T, int sDim>
inline constexpr VectorFieldT<T, sDim> operator^(const VectorFieldT<T, sDim> &x,
                                                 const VectorFieldT<T, sDim> &y) {
  if constexpr (sDim == 3)
    return VectorFieldT<T, 3>(x[1] * y[2] - x[2] * y[1], x[2] * y[0] - x[0] * y[2],
                              x[0] * y[1] - x[1] * y[0]);
  else THROW("Cross product not implemented for dimension other than 3!")
}

template<typename T, int sDim, int tDim>
inline constexpr VectorFieldT<T, sDim> operator^(const VectorFieldT<T, sDim> &x,
                                                 const PointT<T, sDim, tDim> &y) {
  if constexpr (sDim == 3)
    return VectorFieldT<T, 3>(x[1] * y[2] - x[2] * y[1], x[2] * y[0] - x[0] * y[2],
                              x[0] * y[1] - x[1] * y[0]);
  else THROW("Cross product not implemented for dimension other than 3!")
}

template<typename T, int sDim, int tDim>
inline constexpr VectorFieldT<T, sDim> operator^(const PointT<T, sDim, tDim> &x,
                                                 const VectorFieldT<T, sDim> &y) {
  if constexpr (sDim == 3)
    return VectorFieldT<T, 3>(x[1] * y[2] - x[2] * y[1], x[2] * y[0] - x[0] * y[2],
                              x[0] * y[1] - x[1] * y[0]);
  else THROW("Cross product not implemented for dimension other than 3!")
}

template<typename T, int sDim>
inline constexpr VectorFieldT<T, sDim> curl(const VectorFieldT<T, sDim> &x,
                                            const VectorFieldT<T, sDim> &y) {
  return x ^ y;
}

template<typename T, int sDim, int tDim>
inline constexpr VectorFieldT<T, sDim> curl(const VectorFieldT<T, sDim> &x,
                                            const PointT<T, sDim, tDim> &y) {
  return x ^ y;
}

template<typename T, int sDim, int tDim>
inline constexpr VectorFieldT<T, sDim> curl(const PointT<T, sDim, tDim> &x,
                                            const VectorFieldT<T, sDim> &y) {
  return x ^ y;
}

template<typename T, int dim>
constexpr T operator*(const VectorFieldT<T, dim> &x, const VectorFieldT<T, dim> &y) {
  T sum{};
  for (int n = 0; n < dim; ++n)
    sum += x[n] * y[n];
  return sum;
}

template<typename T, int dim, int tDim>
constexpr T operator*(const VectorFieldT<T, dim> &x, const PointT<T, dim, tDim> &y) {
  T sum{};
  for (int n = 0; n < dim; ++n)
    sum += x[n] * y[n];
  return sum;
}

template<typename T, int dim, int tDim>
constexpr T operator*(const PointT<T, dim, tDim> &x, const VectorFieldT<T, dim> &y) {
  T sum{};
  for (int n = 0; n < dim; ++n)
    sum += x[n] * y[n];
  return sum;
}

template<typename T, int dim>
constexpr T operator*(const VectorFieldT<T, dim> &x, const VectorFieldComponentT<T> &y) {
  return x[y.Component()] * y.Value();
}

template<typename T, int dim>
constexpr T operator*(const VectorFieldComponentT<T> &x, const VectorFieldT<T, dim> &y) {
  return x.Value() * y[x.Component()];
}

template<typename T, typename TT, int dim, int tDim>
inline constexpr std::enable_if_t<std::is_scalar<TT>::value, VectorFieldT<T, dim, tDim>>
operator*(const TT &b, const VectorFieldT<T, dim, tDim> &x) {
  VectorFieldT<T, dim, tDim> z(x);
  return z *= b;
}

template<typename T, typename TT, int dim>
inline constexpr std::enable_if_t<std::is_scalar<TT>::value, VectorFieldT<T, dim>>
operator*(const VectorFieldT<T, dim> &x, const TT &b) {
  VectorFieldT<T, dim> z(x);
  return z *= b;
}

template<typename T, typename TT, int dim>
inline constexpr std::enable_if_t<std::is_scalar<TT>::value, VectorFieldT<T, dim>>
operator/(const VectorFieldT<T, dim> &x, const TT &b) {
  VectorFieldT<T, dim> z(x);
  return z /= b;
}

template<typename T, int dim>
inline constexpr T norm(const VectorFieldT<T, dim> &x) {
  return sqrt(abs(x * x));
}

template<typename T, int dim>
inline constexpr T absNorm(const VectorFieldT<T, dim> &x) {
  double norm = 0;
  for (int d = 0; d < x.Dim(); d++) {
    norm += std::abs(x[d]);
  }
  return norm;
}

template<typename T, int dim>
inline constexpr T maxNorm(const VectorFieldT<T, dim> &x) {
  double norm = std::abs(x[0]);
  for (int d = 1; d < x.Dim(); d++) {
    norm = std::max(norm, std::abs(x[d]));
  }
  return norm;
}

/// templates cannot be used for the stream since Saver has the same stucture
template<typename T, int dim>
constexpr inline std::ostream &operator<<(std::ostream &os, const VectorFieldT<T, dim> &x) {
  return x.print(os);
}

template<typename T, int dim>
inline Saver &operator<<(Saver &saver, const VectorFieldT<T, dim> &x) {
  return x.save(saver);
}

template<typename T, int dim>
inline Loader &operator>>(Loader &loader, VectorFieldT<T, dim> &x) {
  return x.load(loader);
}

template<typename T, int sDim, int tDim>
void VFtoPoint(const VectorFieldT<T, sDim> &x, PointT<T, sDim, tDim> &z) {
  for (int n = 0; n < sDim; ++n)
    z[n] = x[n];
}

template<typename T = double, int dim = SpaceDimension>
using VelocityT = VectorFieldT<T, dim>;

using VectorField = VectorFieldT<>;

using Displacement = VectorField;

using Deformation = VectorField;

using Velocity = VectorField;

template<typename T = double, int sDim = SpaceDimension, int tDim = TimeDimension>
using SpaceTimeVectorFieldT = VectorFieldT<T, sDim, tDim>;

using SpaceTimeVectorField = SpaceTimeVectorFieldT<>;

constexpr VectorField zero{};


#ifdef BUILD_IA

using IAVectorField = VectorFieldT<IAInterval, SpaceDimension>;

using IADisplacement = IAVectorField;

using IADeformation = IAVectorField;

using IAVelocity = IAVectorField;

template<int sDim>
VectorFieldT<double, sDim> mid(const VectorFieldT<IAInterval, sDim> &IAx) {
  VectorFieldT<double, sDim> x;
  for (int n = 0; n < sDim; ++n)
    x[n] = mid(IAx[n]);
  return x;
}

template<int sDim>
VectorFieldT<double, sDim> sup(const VectorFieldT<IAInterval, sDim> &IAx) {
  VectorFieldT<double, sDim> x;
  for (int n = 0; n < sDim; ++n)
    x[n] = sup(IAx[n]);
  return x;
}

template<int sDim>
VectorFieldT<double, sDim> inf(const VectorFieldT<IAInterval, sDim> &IAx) {
  VectorFieldT<double, sDim> x;
  for (int n = 0; n < sDim; ++n)
    x[n] = inf(IAx[n]);
  return x;
}

template<int sDim>
inline VectorFieldT<IAInterval, sDim> operator+(const VectorFieldT<double, sDim> &x,
                                                const VectorFieldT<IAInterval, sDim> &y) {
  VectorFieldT<IAInterval, sDim> IAx(x);
  return IAx += y;
}

template<int sDim>
inline VectorFieldT<IAInterval, sDim> operator+(const VectorFieldT<IAInterval, sDim> &x,
                                                const VectorFieldT<double, sDim> &y) {
  VectorFieldT<IAInterval, sDim> IAx(x);
  return IAx += y;
}

template<int sDim>
inline VectorFieldT<IAInterval, sDim> operator-(const VectorFieldT<double, sDim> &x,
                                                const VectorFieldT<IAInterval, sDim> &y) {
  VectorFieldT<IAInterval, sDim> IAx(x);
  return IAx -= y;
}

template<int sDim>
inline VectorFieldT<IAInterval, sDim> operator-(const VectorFieldT<IAInterval, sDim> &x,
                                                const VectorFieldT<double, sDim> &y) {
  VectorFieldT<IAInterval, sDim> IAx(x);
  return IAx -= y;
}

template<int sDim>
inline VectorFieldT<IAInterval, sDim> operator*(const VectorFieldT<double, sDim> &x,
                                                const IAInterval &b) {
  VectorFieldT<IAInterval, sDim> IAx(x);
  return IAx *= b;
}

template<int sDim>
inline VectorFieldT<IAInterval, sDim> operator*(const IAInterval &b,
                                                const VectorFieldT<double, sDim> &x) {
  VectorFieldT<IAInterval, sDim> IAx(x);
  return IAx *= b;
}

template<int sDim>
inline VectorFieldT<IAInterval, sDim> operator/(const VectorFieldT<double, sDim> &x,
                                                const IAInterval &b) {
  VectorFieldT<IAInterval, sDim> IAx(x);
  return IAx /= b;
}

template<int sDim>
inline IAInterval operator*(const VectorFieldT<IAInterval, sDim> &x,
                            const VectorFieldT<double, sDim> &y) {
  VectorFieldT<IAInterval, sDim> IAy(y);
  return x * IAy;
}

template<int sDim>
inline IAInterval operator*(const VectorFieldT<double, sDim> &x,
                            const VectorFieldT<IAInterval, sDim> &y) {
  VectorFieldT<IAInterval, sDim> IAx(x);
  return IAx * y;
}

#endif // BUILD_IA

#endif // of #ifndef _VECTORFIELD_H_
