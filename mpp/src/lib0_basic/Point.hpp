#ifndef _POINT_H_
#define _POINT_H_

#include <cmath>
#include "Assertion.hpp"
#include "GlobalDefinitions.hpp"
#include "Logging.hpp"
#include "SaveLoad.hpp"

#include <iosfwd>
#include <numeric>


#ifdef BUILD_IA

#include "IAInterval.hpp"

#endif

namespace mpp_geometry {
extern double tolerance;

inline void SetTolerance(double tol = GeometricTolerance) { tolerance = tol; }

inline double GetTolerance() { return tolerance; }

inline bool isNear(double a, double b, double tol = GeometricTolerance) { return abs(a - b) < tol; }

#ifdef BUILD_IA
inline bool isNear(const IAInterval &a, const IAInterval &b) {
  return (abs(inf(a) - inf(b)) < tolerance) && (abs(sup(a) - sup(b)) < tolerance);
}
#endif
} // namespace mpp_geometry


template<typename T, int sDim, int tDim>
class PointT;

template<typename T, int sDim, int tDim>
inline constexpr T norm(const PointT<T, sDim, tDim> &x) {
  return sqrt(x * x);
}

template<typename T, int sDim, int tDim>
inline constexpr T normST(const PointT<T, sDim, tDim> &x) {
  T sum{};
  for (int n = 0; n < sDim + tDim; ++n)
    sum += x[n] * x[n];
  return sqrt(sum);
}

template<typename T = double, int sDim = SpaceDimension, int tDim = TimeDimension>
class PointT {
public:
  template<typename dummy = T>
  static typename std::enable_if<std::is_same_v<dummy, double>, PointT>::type
  Random(double leftBound, double rightBound) {
    PointT p;
    for (int i = 0; i < sDim + tDim; ++i) {
      p[i] = leftBound + (rightBound - leftBound) * (double(rand()) / RAND_MAX);
    }
    return p;
  }
private:
  T z[sDim + tDim]{};
  mutable size_t hash{0};
public:
  constexpr PointT() {
    for (int n = 0; n < sDim + tDim; ++n)
      z[n] = T{};
  }

  template<typename TT>
  constexpr PointT(const PointT<TT, sDim, tDim> &y) {
    for (int n = 0; n < sDim + tDim; ++n)
      z[n] = y[n];
  }

  constexpr PointT(const T *a) {
    for (int n = 0; n < sDim + tDim; ++n)
      z[n] = a[n];
  }

  constexpr PointT(std::vector<T> a) {
    for (int n = 0; n < std::min(sDim + tDim, (int)a.size()); ++n)
      z[n] = a[n];
  }

  constexpr PointT(T a, T b = T{}, T c = T{}, T d = T{}) {
    if constexpr (sDim + tDim > 0) z[0] = a;
    if constexpr (sDim + tDim > 1) z[1] = b;
    if constexpr (sDim + tDim > 2) z[2] = c;
    if constexpr (sDim + tDim > 3) z[3] = d;
  }

  constexpr PointT(const PointT &y, const T &t) {
    memcpy(z, y.z, sDim * sizeof(T));
    if constexpr (tDim > 0) z[sDim] = t;
  }

  const T &operator[](const unsigned int i) const { return z[i]; }

  T &operator[](const unsigned int i) {
    hash = 0;
    return z[i];
  }

  const T *operator()() const {
    hash = 0;
    return z;
  }

  constexpr T t() const {
    if constexpr (tDim > 0) {
      return z[sDim];
    } else {
      return 0.0;
    }
  }

  PointT<T, sDim, tDim> &WithT(T t) {
    if constexpr (tDim > 0) {
      hash = 0;
      z[sDim] = t;
      return *this;
    }
    THROW("Point not time dependent!")
  }

  PointT<T, sDim, tDim> CopyWithT(T t) const {
    if constexpr (tDim > 0) { return PointT<T, sDim, tDim>(*this).WithT(t); }
    THROW("Point not time dependent!")
  }

  constexpr PointT<T, sDim, tDim> GetSpatialPoint() {
    if constexpr (sDim == 1) {
      return PointT<T, sDim, tDim>(z[0]);
    } else if constexpr (sDim == 2) {
      return PointT<T, sDim, tDim>(z[0], z[1]);
    } else if constexpr (sDim == 3) {
      return PointT<T, sDim, tDim>(z[0], z[1], z[2]);
    }
    Exit("Not implemented for dim = " + std::to_string(sDim));
  }

  constexpr bool isTimeDep() const { return tDim == 1; }

  constexpr int SpaceDim() const { return sDim; }

  constexpr int TimeDim() const { return tDim; }

  constexpr PointT &operator=(const T &a) {
    hash = 0;
    for (int n = 0; n < sDim + tDim; ++n)
      z[n] = a;
    return *this;
  }

  template<typename TT>
  constexpr PointT &operator=(const PointT<TT, sDim, tDim> &y) {
    hash = 0;
    for (int n = 0; n < sDim + tDim; ++n)
      z[n] = y[n];
    return *this;
  }

  template<typename TT>
  constexpr PointT &operator+=(const PointT<TT, sDim, tDim> &y) {
    hash = 0;
    for (int n = 0; n < sDim + tDim; ++n)
      z[n] += y[n];
    return *this;
  }

  template<typename TT>
  constexpr PointT &operator-=(const PointT<TT, sDim, tDim> &y) {
    hash = 0;
    for (int n = 0; n < sDim + tDim; ++n)
      z[n] -= y[n];
    return *this;
  }

  template<typename TT>
  constexpr PointT &operator*=(const TT &a) {
    hash = 0;
    for (int n = 0; n < sDim + tDim; ++n)
      z[n] *= a;
    return *this;
  }

  template<typename TT>
  constexpr PointT &operator/=(const TT &a) {
    hash = 0;
    for (int n = 0; n < sDim + tDim; ++n)
      z[n] /= a;
    return *this;
  }

  constexpr bool operator==(const PointT &other) const {
    for (int n = 0; n < sDim + tDim; ++n) {
      if (!mpp_geometry::isNear(z[n], other.z[n])) return false;
    }

    return true;
  }

  constexpr bool operator!=(const PointT &y) const { return !(*this == y); }

  constexpr bool isNear(const PointT &other, double tol) const {
    for (int n = 0; n < sDim + tDim; ++n) {
      if (!mpp_geometry::isNear(z[n], other.z[n], tol)) return false;
    }

    return true;
  }

  constexpr bool operator<(const PointT &other) const {
    if constexpr (tDim > 0) {
      if (z[sDim] < other.z[sDim] - GeometricTolerance) return true;
      if (z[sDim] > other.z[sDim] + GeometricTolerance) return false;
    }

    for (int n = 0; n < sDim; ++n) {
      if (z[n] < other.z[n] - GeometricTolerance) return true;
      if (z[n] > other.z[n] + GeometricTolerance) return false;
    }

    return false;
  }

  constexpr bool operator>(const PointT &other) const {
    if constexpr (tDim > 0) {
      if (z[sDim] > other.z[sDim] + GeometricTolerance) return true;
      if (z[sDim] < other.z[sDim] - GeometricTolerance) return false;
    }

    for (int n = 0; n < sDim; ++n) {
      if (z[n] > other.z[n] + GeometricTolerance) return true;
      if (z[n] < other.z[n] - GeometricTolerance) return false;
    }

    return false;
  }

  template<typename OSTREAM>
  constexpr OSTREAM &print(OSTREAM &os) const {
    if constexpr (std::is_same<OSTREAM, MasterLogging>::value) os << beginD;
    os << "(" << z[0];
    for (int i = 1; i < sDim; ++i)
      os << ", " << z[i];
    if (tDim == 1) os << " | " << z[sDim];
    os << ")";
    if constexpr (std::is_same<OSTREAM, MasterLogging>::value) os << endD;
    return os;
  }

  Saver &save(Saver &saver) const {
    saver << (sDim + tDim);
    for (int n = 0; n < sDim + tDim; ++n)
      saver << z[n];
    return saver;
  }

  Loader &load(Loader &loader) {
    hash = 0;
    int dim;
    loader >> dim;
    if (sDim + tDim != dim) THROW("Point Dimension does not fit")
    for (int n = 0; n < sDim + tDim; ++n)
      loader >> z[n];
    return loader;
  }

  size_t Hash() const;

  void Abs() {
    for (int n = 0; n < sDim + tDim; ++n) {
      z[n] = abs(z[n]);
    }
  }

  T Sum() const { return std::accumulate(std::begin(z), std::end(z), T(0.0)); }

  constexpr explicit operator double() const { return norm(*this); }
};

template<typename T, int sDim, int tDim>
size_t PointT<T, sDim, tDim>::Hash() const {
  if (z[0] == infty) return 1;
  if (hash == 0) {
    double hash_sum = 115421.1 * z[0] + 156421541.5;
    if constexpr (sDim + tDim >= 2) hash_sum += 124323.3 * z[1];
    if constexpr (sDim + tDim >= 3) hash_sum += 221313.7 * z[2];
    if constexpr (sDim + tDim >= 4) hash_sum += 3.7332 * z[3];
    hash = size_t(std::abs(hash_sum) + 345453 * GeometricTolerance);
  }
  return hash;
}

namespace std {

template<typename T, int sDim, int tDim>
struct hash<PointT<T, sDim, tDim>> {
  size_t operator()(const PointT<T, sDim, tDim> &P) const { return P.Hash(); }
};
} // namespace std

template<typename T, int sDim, int tDim>
inline constexpr PointT<T, sDim, tDim> operator+(const PointT<T, sDim, tDim> &x,
                                                 const PointT<T, sDim, tDim> &y) {
  PointT<T, sDim, tDim> z(x);
  return z += y;
}

template<typename T, int sDim, int tDim>
inline constexpr PointT<T, sDim, tDim> operator-(const PointT<T, sDim, tDim> &x,
                                                 const PointT<T, sDim, tDim> &y) {
  PointT<T, sDim, tDim> z(x);
  return z -= y;
}

template<typename T, int sDim, int tDim>
inline constexpr PointT<T, sDim, tDim> operator-(const PointT<T, sDim, tDim> &x) {
  PointT<T, sDim, tDim> z(x);
  return z *= T(-1.0);
}

template<typename T, typename TT, int sDim, int tDim>
inline constexpr PointT<T, sDim, tDim> operator*(const TT &b, const PointT<T, sDim, tDim> &x) {
  PointT<T, sDim, tDim> z(x);
  return z *= b;
}

template<typename T, typename TT, int sDim, int tDim>
inline constexpr PointT<T, sDim, tDim> operator*(const PointT<T, sDim, tDim> &x, const TT &b) {
  PointT<T, sDim, tDim> z(x);
  return z *= b;
}

template<typename T, typename TT, int sDim, int tDim>
inline constexpr PointT<T, sDim, tDim> operator/(const PointT<T, sDim, tDim> &x, const TT &b) {
  PointT<T, sDim, tDim> z(x);
  return z /= b;
}

template<typename T, int sDim, int tDim>
inline constexpr T operator*(const PointT<T, sDim, tDim> &x, const PointT<T, sDim, tDim> &y) {
  T sum{};
  for (int n = 0; n < sDim; ++n)
    sum += x[n] * y[n];
  return sum;
}

template<typename T, int sDim, int tDim>
inline constexpr T dist(const PointT<T, sDim, tDim> &x, const PointT<T, sDim, tDim> &y) {
  return norm(x - y);
}

template<typename T, int sDim, int tDim>
inline constexpr T distInfty(const PointT<T, sDim, tDim> &x, const PointT<T, sDim, tDim> &y) {
  T m = std::abs(x[0] - y[0]);
  for (int i = 1; i < sDim; i++) {
    m = std::max(m, std::abs(x[i] - y[i]));
  }
  return m;
}

template double distInfty(const PointT<double, SpaceDimension, TimeDimension> &x,
                          const PointT<double, SpaceDimension, TimeDimension> &y);

template<typename T, int sDim, int tDim>
inline constexpr PointT<T, sDim, tDim> operator^(const PointT<T, sDim, tDim> &x,
                                                 const PointT<T, sDim, tDim> &y) {
  if constexpr (sDim == 3)
    return PointT<T, 3, tDim>(x[1] * y[2] - x[2] * y[1], x[2] * y[0] - x[0] * y[2],
                              x[0] * y[1] - x[1] * y[0]);
  else THROW("Cross product not implemented for points in dimension other than 3!")
}

template<typename T, int sDim, int tDim>
inline constexpr PointT<T, sDim, tDim> curl(const PointT<T, sDim, tDim> &x,
                                            const PointT<T, sDim, tDim> &y) {
  return x ^ y;
}

template<typename T, int sDim, int tDim>
inline constexpr T det(const PointT<T, sDim, tDim> &x, const PointT<T, sDim, tDim> &y,
                       const PointT<T, sDim, tDim> &z) {
  if constexpr (sDim == 3) return x * curl(y, z);
  else THROW("Determinant not implemented for points in dimension other than 3!")
}

template<int sDim, int tDim>
inline constexpr PointT<double, sDim, tDim> mid(const PointT<double, sDim, tDim> &x) {
  return x;
}

template<int sDim, int tDim>
inline constexpr PointT<double, sDim, tDim> abs(const PointT<double, sDim, tDim> &x) {
  auto y = x;
  y.Abs();
  return y;
}

template<int sDim, int tDim>
inline constexpr PointT<double, sDim, tDim> sum(const PointT<double, sDim, tDim> &x) {
  return x.Sum();
}

template<typename T, int sDim, int tDim>
inline constexpr bool LessPoint_x(const PointT<T, sDim, tDim> &x, const PointT<T, sDim, tDim> &y) {
  if (x[0] < x[0] - GeometricTolerance) return true;
  if (x[0] > x[0] + GeometricTolerance) return false;
  if (x[1] < x[1] - GeometricTolerance) return true;
  if (x[1] > x[1] + GeometricTolerance) return false;
  if (x[2] < x[2] - GeometricTolerance) return true;
  if constexpr (tDim > 0) {
    if (x[sDim] < y[sDim] - GeometricTolerance) return true;
    if (x[sDim] > y[sDim] + GeometricTolerance) return false;
  }
  return false;
}

template<typename T, int sDim, int tDim>
inline constexpr bool LessPoint_y(const PointT<T, sDim, tDim> &x, const PointT<T, sDim, tDim> &y) {
  if (x[1] < y[1] - GeometricTolerance) return true;
  if (x[1] > y[1] + GeometricTolerance) return false;
  if (x[2] < y[2] - GeometricTolerance) return true;
  if (x[2] > y[2] + GeometricTolerance) return false;
  if (x[0] < y[0] - GeometricTolerance) return true;
  if constexpr (tDim > 0) {
    if (x[sDim] < y[sDim] - GeometricTolerance) return true;
    if (x[sDim] > y[sDim] + GeometricTolerance) return false;
  }
  return false;
}

template<typename T, int sDim, int tDim>
inline constexpr bool LessPoint_z(const PointT<T, sDim, tDim> &x, const PointT<T, sDim, tDim> &y) {
  if (x[2] < y[2] - GeometricTolerance) return true;
  if (x[2] > y[2] + GeometricTolerance) return false;
  if (x[0] < y[0] - GeometricTolerance) return true;
  if (x[0] > y[0] + GeometricTolerance) return false;
  if (x[1] < y[1] - GeometricTolerance) return true;
  if (x[1] < y[1] - GeometricTolerance) return false;
  if constexpr (tDim > 0) {
    if (x[sDim] < y[sDim] - GeometricTolerance) return true;
    if (x[sDim] > y[sDim] + GeometricTolerance) return false;
  }
  return false;
}

template<typename T, int sDim, int tDim>
inline Saver &operator<<(Saver &saver, const PointT<T, sDim, tDim> &x) {
  return x.save(saver);
}

template<typename T, int sDim, int tDim>
inline Loader &operator>>(Loader &loader, PointT<T, sDim, tDim> &x) {
  return x.load(loader);
}

/// templates cannot be used for the stream since Saver has the same stucture
template<typename T, int sDim, int tDim>
constexpr inline std::ostream &operator<<(std::ostream &os, const PointT<T, sDim, tDim> &z) {
  return z.print(os);
}

using Point = PointT<double, SpaceDimension, TimeDimension>;

template PointT<> &Point::operator/= <size_t>(const size_t &a);

constexpr Point Origin{};

constexpr Point Infty(infty, infty, infty, infty);

#ifdef BUILD_IA

typedef PointT<IAInterval, SpaceDimension, TimeDimension> IAPoint;

template<int sDim, int tDim>
PointT<double, sDim, tDim> mid(const PointT<IAInterval, sDim, tDim> &IAx) {
  PointT<double, sDim, tDim> x;
  for (int n = 0; n < sDim + tDim; ++n)
    x[n] = mid(IAx[n]);
  return x;
}

template<int sDim, int tDim>
PointT<double, sDim, tDim> sup(const PointT<IAInterval, sDim, tDim> &IAx) {
  PointT<double, sDim, tDim> x;
  for (int n = 0; n < sDim + tDim; ++n)
    x[n] = sup(IAx[n]);
  return x;
}

template<int sDim, int tDim>
PointT<double, sDim, tDim> inf(const PointT<IAInterval, sDim, tDim> &IAx) {
  PointT<double, sDim, tDim> x;
  for (int n = 0; n < sDim + tDim; ++n)
    x[n] = inf(IAx[n]);
  return x;
}

template<int sDim, int tDim>
inline PointT<IAInterval, sDim, tDim> operator+(const PointT<double, sDim, tDim> &x,
                                                const PointT<IAInterval, sDim, tDim> &y) {
  PointT<IAInterval, sDim, tDim> IAx(x);
  return IAx += y;
}

template<int sDim, int tDim>
inline PointT<IAInterval, sDim, tDim> operator+(const PointT<IAInterval, sDim, tDim> &x,
                                                const PointT<double, sDim, tDim> &y) {
  PointT<IAInterval, sDim, tDim> IAx(x);
  return IAx += y;
}

template<int sDim, int tDim>
inline PointT<IAInterval, sDim, tDim> operator-(const PointT<double, sDim, tDim> &x,
                                                const PointT<IAInterval, sDim, tDim> &y) {
  PointT<IAInterval, sDim, tDim> IAx(x);
  return IAx -= y;
}

template<int sDim, int tDim>
inline PointT<IAInterval, sDim, tDim> operator-(const PointT<IAInterval, sDim, tDim> &x,
                                                const PointT<double, sDim, tDim> &y) {
  PointT<IAInterval, sDim, tDim> IAx(x);
  return IAx -= y;
}

template<int sDim, int tDim>
inline PointT<IAInterval, sDim, tDim> operator*(const PointT<double, sDim, tDim> &x,
                                                const IAInterval &b) {
  PointT<IAInterval, sDim, tDim> IAx(x);
  return IAx *= b;
}

template<int sDim, int tDim>
inline PointT<IAInterval, sDim, tDim> operator*(const IAInterval &b,
                                                const PointT<double, sDim, tDim> &x) {
  PointT<IAInterval, sDim, tDim> IAx(x);
  return IAx *= b;
}

template<int sDim, int tDim>
inline PointT<IAInterval, sDim, tDim> operator/(const PointT<double, sDim, tDim> &x,
                                                const IAInterval &b) {
  PointT<IAInterval, sDim, tDim> IAx(x);
  return IAx /= b;
}

template<int sDim, int tDim>
inline IAInterval operator*(const PointT<IAInterval, sDim, tDim> &x,
                            const PointT<double, sDim, tDim> &y) {
  PointT<IAInterval, sDim, tDim> IAy(y);
  return x * IAy;
}

template<int sDim, int tDim>
inline IAInterval operator*(const PointT<double, sDim, tDim> &x,
                            const PointT<IAInterval, sDim, tDim> &y) {
  PointT<IAInterval, sDim, tDim> IAx(x);
  return IAx * y;
}

#endif // BUILD_IA

/**
 * Points is used in the initialization of the cells
 */
class Points : public std::vector<Point> {
public:
  Points(int n) : std::vector<Point>(n) {}

  Points(int n, const Point *x) : std::vector<Point>(n) {
    for (int i = 0; i < n; ++i)
      (*this)[i] = x[i];
  }

  Points(int n1, const Point *x1, int n2, const Point *x2) : std::vector<Point>(n1 + n2) {
    for (int i = 0; i < n1; ++i)
      (*this)[i] = x1[i];
    for (int i = 0; i < n2; ++i)
      (*this)[n1 + i] = x2[i];
  }

  Points(int n1, const Point *x1, int n2, const Point *x2, int n3, const Point *x3, int n4,
         const Point *x4) : std::vector<Point>(n1 + n2 + n3 + n4) {
    for (int i = 0; i < n1; ++i)
      (*this)[i] = x1[i];
    for (int i = 0; i < n2; ++i)
      (*this)[n1 + i] = x2[i];
    for (int i = 0; i < n3; ++i)
      (*this)[n1 + n2 + i] = x3[i];
    for (int i = 0; i < n4; ++i)
      (*this)[n1 + n2 + n3 + i] = x4[i];
  }
};

//==================================================================================================
// LEGACY CODE
//==================================================================================================

void Point_2d();

void Point_3d();

void Point_1d_t();

void Point_2d_t();

void Point_3d_t();


#endif // of #ifndef _POINT_H_
