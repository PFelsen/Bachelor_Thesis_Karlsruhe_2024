#ifndef CVECTOR_H
#define CVECTOR_H

#include "RVector.hpp"

#include <complex>

#ifdef BUILD_IA

#include "IACInterval.hpp"

template<typename T>
using COMPLEX_TYPE = typename std::conditional<
    std::is_same_v<T, double>, std::complex<double>,
    typename std::conditional<std::is_same_v<T, IAInterval>, IACInterval, void>::type>::type;

#else

template<typename T>
using COMPLEX_TYPE =
    typename std::conditional<std::is_same_v<T, double>, std::complex<double>, void>::type;

#endif

namespace std {
inline std::complex<double> operator*(const std::complex<double> &a, int b) {
  return a * double(b);
}

inline std::complex<double> operator*(int a, const std::complex<double> &b) {
  return double(a) * b;
}

inline std::complex<double> operator/(const std::complex<double> &a, int b) {
  return a / double(b);
}
} // namespace std

template<typename REAL = double>
class CVectorT {
  using COMPLEX = COMPLEX_TYPE<REAL>;
protected:
  std::vector<COMPLEX> z;
public:
  CVectorT() = default;

  explicit CVectorT(int length);

  template<typename REAL1>
  explicit CVectorT(const CVectorT<REAL1> &v) : z(v.size()) {
    for (int i = 0; i < z.size(); ++i)
      z[i] = v[i];
  }

  template<typename REAL1>
  explicit CVectorT(const RVectorT<REAL1> &v) : z(v.size()) {
    for (int i = 0; i < z.size(); ++i)
      z[i] = v[i];
  }

  template<typename COMPLEX1>
  CVectorT(const std::vector<COMPLEX1> &v) : z(v.begin(), v.end()) {}

  template<typename T>
  constexpr CVectorT(const T &a, int length) : z(length) {
    for (int i = 0; i < z.size(); ++i)
      z[i] = a;
  }

  template<typename T>
  CVectorT &operator=(const T &a) {
    for (int i = 0; i < z.size(); ++i)
      z[i] = a;
    return *this;
  }

  template<typename REAL1>
  CVectorT &operator=(const CVectorT<REAL1> &v) {
    z.resize(v.size());
    for (int i = 0; i < z.size(); ++i)
      z[i] = v[i];
    return *this;
  }

  template<typename REAL1>
  CVectorT &operator=(const RVectorT<REAL1> &v) {
    z.resize(v.size());
    for (int i = 0; i < z.size(); ++i)
      z[i] = v[i];
    return *this;
  }

  template<typename T>
  CVectorT &operator+=(const T &a) {
    for (int i = 0; i < z.size(); ++i)
      z[i] += a;
    return *this;
  }

  template<typename REAL1>
  CVectorT &operator+=(const CVectorT<REAL1> &v) {
    if (z.size() != v.size()) THROW("Size does not fit!")
    for (int i = 0; i < z.size(); ++i)
      z[i] += v[i];
    return *this;
  }

  template<typename REAL1>
  CVectorT &operator+=(const RVectorT<REAL1> &v) {
    if (z.size() != v.size()) THROW("Size does not fit!")
    for (int i = 0; i < z.size(); ++i)
      z[i] += v[i];
    return *this;
  }

  template<typename T>
  CVectorT &operator-=(const T &a) {
    for (int i = 0; i < z.size(); ++i)
      z[i] -= a;
    return *this;
  }

  template<typename REAL1>
  CVectorT &operator-=(const CVectorT<REAL1> &v) {
    if (z.size() != v.size()) THROW("Size does not fit!")
    for (int i = 0; i < z.size(); ++i)
      z[i] -= v[i];
    return *this;
  }

  template<typename REAL1>
  CVectorT &operator-=(const RVectorT<REAL1> &v) {
    if (z.size() != v.size()) THROW("Size does not fit!")
    for (int i = 0; i < z.size(); ++i)
      z[i] -= v[i];
    return *this;
  }

  template<typename T>
  CVectorT &operator*=(const T &a) {
    for (int i = 0; i < z.size(); ++i)
      z[i] *= a;
    return *this;
  }

  template<typename REAL1>
  CVectorT &operator*=(const RVectorT<REAL1> &v) {
    for (int i = 0; i < z.size(); ++i)
      z[i] *= v[i];
    return *this;
  }

  // This a component wise operation
  template<typename REAL1>
  CVectorT &operator*=(const CVectorT<REAL1> &v) {
    for (int i = 0; i < z.size(); ++i)
      z[i] *= v[i];
    return *this;
  }

  template<typename T>
  CVectorT &operator/=(const T &a) {
    for (int i = 0; i < z.size(); ++i)
      z[i] /= a;
    return *this;
  }

  const std::vector<COMPLEX> &asVector() const { return z; }

  std::vector<COMPLEX> &asVector() { return z; }

  int Dim() const { return z.size(); }

  int size() const { return z.size(); }

  void resize(int N) { z.resize(N); }

  COMPLEX &operator[](int i) { return z[i]; }

  const COMPLEX &operator[](int i) const { return z[i]; }

  auto begin() const { return z.begin(); }

  auto begin() { return z.begin(); }

  auto rbegin() const { return z.rbegin(); }

  auto end() const { return z.end(); }

  auto end() { return z.end(); }

  auto rend() const { return z.rend(); }

  auto push_back(const COMPLEX &a) { return z.push_back(a); }

  void clear() { z.clear(); }

  auto append(const CVectorT &v) { return z.insert(z.end(), v.begin(), v.end()); }

  CVectorT &conj();

  RVectorT<REAL> real() const;

  RVectorT<REAL> imag() const;

  REAL normSqr() const;

  REAL norm() const;

  COMPLEX Mean() const;

  REAL Variance() const;

  void Accumulate(int commSplit = 0);

  Saver &save(Saver &saver) const;

  Loader &load(Loader &loader);
};

template<typename REAL>
bool operator==(const CVectorT<REAL> &x, const CVectorT<REAL> &y) {
  if (x.size() != y.size()) return false;
  return (x.real() == y.real()) && (x.imag() == y.imag());
}

template<typename REAL, typename REAL1>
bool operator==(const CVectorT<REAL> &x, const RVectorT<REAL1> &y) {
  return x == CVectorT<REAL>(y);
}

template<typename REAL, typename REAL1>
bool operator==(const RVectorT<REAL1> &x, const CVectorT<REAL> &y) {
  return CVectorT<REAL>(x) == y;
}

template<typename REAL>
bool operator!=(const CVectorT<REAL> &x, const CVectorT<REAL> &y) {
  return !(x == y);
}

template<typename REAL, typename REAL1>
bool operator!=(const CVectorT<REAL> &x, const RVectorT<REAL1> &y) {
  return !(x == y);
}

template<typename REAL, typename REAL1>
bool operator!=(const RVectorT<REAL1> &x, const CVectorT<REAL> &y) {
  return !(x == y);
}

template<typename REAL>
inline CVectorT<REAL> operator-(const CVectorT<REAL> &v) {
  CVectorT<REAL> z(v);
  return z *= -1.0;
}

template<typename REAL>
inline CVectorT<REAL> operator+(const CVectorT<REAL> &v, const CVectorT<REAL> &w) {
  CVectorT<REAL> z(v);
  return z += w;
}

template<typename REAL>
inline CVectorT<REAL> operator+(const CVectorT<REAL> &v, const RVectorT<REAL> &w) {
  CVectorT<REAL> z(v);
  return z += w;
}

template<typename REAL>
inline CVectorT<REAL> operator+(const RVectorT<REAL> &v, const CVectorT<REAL> &w) {
  CVectorT<REAL> z(v);
  return z += w;
}

template<typename REAL>
inline CVectorT<REAL> operator-(const CVectorT<REAL> &v, const CVectorT<REAL> &w) {
  CVectorT<REAL> z(v);
  return z -= w;
}

template<typename REAL>
inline CVectorT<REAL> operator-(const CVectorT<REAL> &v, const RVectorT<REAL> &w) {
  CVectorT<REAL> z(v);
  return z -= w;
}

template<typename REAL>
inline CVectorT<REAL> operator-(const RVectorT<REAL> &v, const CVectorT<REAL> &w) {
  CVectorT<REAL> z(v);
  return z -= w;
}

template<typename REAL, typename REAL1>
inline CVectorT<REAL> operator*(const REAL1 &b, const CVectorT<REAL> &v) {
  CVectorT<REAL> z(v);
  return z *= b;
}

template<typename REAL, typename REAL1>
inline CVectorT<REAL> operator*(const CVectorT<REAL> &v, const REAL1 &b) {
  CVectorT<REAL> z(v);
  return z *= b;
}

template<typename REAL, typename REAL1>
inline CVectorT<REAL> operator/(const CVectorT<REAL> &v, const REAL1 &b) {
  CVectorT<REAL> z(v);
  return z /= b;
}

template<typename REAL>
COMPLEX_TYPE<REAL> operator*(const CVectorT<REAL> &v, const CVectorT<REAL> &w) {
  if (v.size() != w.size()) THROW("Size does not fit!")
  COMPLEX_TYPE<REAL> sum{};
  for (int i = 0; i < v.size(); ++i)
    sum += v[i] * conj(w[i]);
  return sum;
}

template<typename REAL>
COMPLEX_TYPE<REAL> operator*(const CVectorT<REAL> &v, const RVectorT<REAL> &w) {
  if (v.size() != w.size()) THROW("Size does not fit!")
  COMPLEX_TYPE<REAL> sum{};
  for (int i = 0; i < v.size(); ++i)
    sum += v[i] * w[i];
  return sum;
}

template<typename REAL>
COMPLEX_TYPE<REAL> operator*(const RVectorT<REAL> &v, const CVectorT<REAL> &w) {
  COMPLEX_TYPE<REAL> sum{};
  for (int i = 0; i < v.size(); ++i)
    sum += v[i] * conj(w[i]);
  return sum;
}

template<typename REAL>
inline CVectorT<REAL> conj(const CVectorT<REAL> &v) {
  CVectorT<REAL> w(v);
  return w.conj();
}

template<typename REAL>
inline RVectorT<REAL> real(const CVectorT<REAL> &v) {
  return v.real();
}

template<typename REAL>
inline RVectorT<REAL> imag(const CVectorT<REAL> &v) {
  return v.imag();
}

template<typename REAL>
inline REAL normSqr(const CVectorT<REAL> &v) {
  return v.normSqr();
}

template<typename REAL>
inline REAL norm(const CVectorT<REAL> &v) {
  return v.norm();
}

template<typename REAL>
inline Saver &operator<<(Saver &saver, const CVectorT<REAL> &v) {
  return v.save(saver);
}

template<typename REAL>
inline Loader &operator>>(Loader &loader, CVectorT<REAL> &v) {
  return v.load(loader);
}

template<typename REAL>
std::ostream &operator<<(std::ostream &os, const CVectorT<REAL> &v) {
  for (int i = 0; i < v.size(); ++i) {
    if (i == 0) os << v[i];
    else os << "\n" << v[i];
  }
  return os;
}

using CVector = CVectorT<>;

#ifdef BUILD_IA

using IACVector = CVectorT<IAInterval>;

CVector mid(const IACVector &);

inline IACVector operator+(const IACVector &v, const CVector &w) {
  IACVector z(v);
  return z += w;
}

inline IACVector operator+(const CVector &v, const IACVector &w) {
  IACVector z(v);
  return z += w;
}

inline IACVector operator+(const IACVector &v, const RVector &w) {
  IACVector z(v);
  return z += w;
}

inline IACVector operator+(const RVector &v, const IACVector &w) {
  IACVector z(v);
  return z += w;
}

inline IACVector operator-(const IACVector &v, const CVector &w) {
  IACVector z(v);
  return z -= w;
}

inline IACVector operator-(const CVector &v, const IACVector &w) {
  IACVector z(v);
  return z -= w;
}

inline IACVector operator-(const IACVector &v, const RVector &w) {
  IACVector z(v);
  return z -= w;
}

inline IACVector operator-(const RVector &v, const IACVector &w) {
  IACVector z(v);
  return z -= w;
}

inline IACInterval operator*(const IACVector &v, const CVector &w) {
  if (v.size() != w.size()) THROW("Size does not fit!")
  IACInterval sum;
  for (int i = 0; i < v.size(); ++i)
    sum += v[i] * conj(w[i]);
  return sum;
}

inline IACInterval operator*(const CVector &v, const IACVector &w) {
  if (v.size() != w.size()) THROW("Size does not fit!")
  IACInterval sum;
  for (int i = 0; i < v.size(); ++i)
    sum += v[i] * conj(w[i]);
  return sum;
}

inline IACInterval operator*(const CVector &v, const IARVector &w) {
  if (v.size() != w.size()) THROW("Size does not fit!")
  IACInterval sum;
  for (int i = 0; i < v.size(); ++i)
    sum += IACInterval(v[i]) * w[i];
  return sum;
}

inline IACInterval operator*(const IACVector &v, const RVector &w) {
  if (v.size() != w.size()) THROW("Size does not fit!")
  IACInterval sum;
  for (int i = 0; i < v.size(); ++i)
    sum += v[i] * w[i];
  return sum;
}

inline IACInterval operator*(const RVector &v, const IACVector &w) {
  if (v.size() != w.size()) THROW("Size does not fit!")
  IACInterval sum;
  for (int i = 0; i < v.size(); ++i)
    sum += v[i] * conj(w[i]);
  return sum;
}

#endif // BUILD_IA

#endif
