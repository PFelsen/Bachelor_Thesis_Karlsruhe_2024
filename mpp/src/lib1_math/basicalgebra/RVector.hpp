#ifndef RVector_H
#define RVector_H

#include <vector>

#include "Assertion.hpp"
#include "Logging.hpp"


#ifdef BUILD_IA

#include "IAInterval.hpp"

#endif

namespace mpp_ba {
void SetTolerance(double tol = 1e-12);

double GetTolerance();

bool isNear(double a, double b);

bool isNearZero(double a);

#ifdef BUILD_IA

bool isNear(const IAInterval &a, const IAInterval &b);

bool isNearZero(const IAInterval &a);

#endif

template<typename REAL, typename REAL1>
using enable_if_not_same =
    typename std::enable_if<!std::is_same<typename std::decay<REAL1>::type, REAL>::value>::type;

template<typename REAL, typename REAL1, typename RETURN_TYPE>
using enable_if_constructible =
    typename std::enable_if<std::is_constructible<REAL, REAL1>::value, RETURN_TYPE>::type;
} // namespace mpp_ba

class Saver;

class Loader;

template<typename REAL = double>
class RVectorT {
protected:
  std::vector<REAL> z{};
public:
  RVectorT() {}

  explicit RVectorT(int length) : z(length) {}

  RVectorT(const RVectorT &v) : z(v.z) {}

  RVectorT(RVectorT &&v) : z(std::move(v.z)) {}

  /// RVector of other type
  template<typename REAL1, typename Constraint = mpp_ba::enable_if_not_same<REAL, REAL1>>
  RVectorT(const RVectorT<REAL1> &v) : z(v.size()) {
    for (int i = 0; i < v.size(); ++i)
      z[i] = v[i];
  }

  /// selects part of v (indices included)
  RVectorT(const RVectorT &v, int startIdx, int endIdx) :
      z(v.begin() + startIdx, v.begin() + endIdx) {}

  template<typename REAL1, typename Constraint = mpp_ba::enable_if_constructible<REAL, REAL1, REAL>>
  RVectorT(const REAL1 &a, int length) : z(length, a) {}

  template<typename REAL1, typename Constraint = mpp_ba::enable_if_constructible<REAL, REAL1, REAL>>
  RVectorT(const std::initializer_list<REAL1> &v) : z(v.begin(), v.end()) {}

  template<typename REAL1, typename Constraint = mpp_ba::enable_if_constructible<REAL, REAL1, REAL>>
  RVectorT(const std::vector<REAL1> &v) : z(v.begin(), v.end()) {}

  template<typename REAL1>
  mpp_ba::enable_if_constructible<REAL, REAL1, RVectorT> &operator=(const REAL1 &b) {
    std::fill(z.begin(), z.end(), b);
    return *this;
  }

  RVectorT &operator=(const RVectorT &v) {
    if (this != &v) z = v.z;
    return *this;
  }

  RVectorT &operator=(RVectorT &&v) {
    z = std::move(v.z);
    return *this;
  }

  /// RVector of other type
  template<typename REAL1, typename Constraint = mpp_ba::enable_if_not_same<REAL, REAL1>>
  RVectorT &operator=(const RVectorT<REAL1> &v) {
    z.resize(v.size());
    for (int i = 0; i < z.size(); ++i)
      z[i] = v[i];
    return *this;
  }

  template<typename REAL1>
  mpp_ba::enable_if_constructible<REAL, REAL1, RVectorT> &operator+=(const REAL1 &a) {
    for (int i = 0; i < z.size(); ++i)
      z[i] += a;
    return *this;
  }

  template<typename REAL1>
  RVectorT &operator+=(const RVectorT<REAL1> &v) {
    if constexpr (DebugLevel > 0) {
      if (z.size() != v.size()) THROW("Size does not fit!")
    }
    for (int i = 0; i < z.size(); ++i)
      z[i] += v[i];
    return *this;
  }

  template<typename REAL1>
  mpp_ba::enable_if_constructible<REAL, REAL1, RVectorT> &operator-=(const REAL1 &a) {
    for (int i = 0; i < z.size(); ++i)
      z[i] -= a;
    return *this;
  }

  template<typename REAL1>
  RVectorT &operator-=(const RVectorT<REAL1> &v) {
    if constexpr (DebugLevel > 0) {
      if (z.size() != v.size()) THROW("Size does not fit!")
    }
    for (int i = 0; i < z.size(); ++i)
      z[i] -= v[i];
    return *this;
  }

  template<typename REAL1>
  mpp_ba::enable_if_constructible<REAL, REAL1, RVectorT> &operator*=(const REAL1 &a) {
    if constexpr (std::is_same_v<REAL, double>) {
      if (a == 1.0) return *this;
    }
    for (int i = 0; i < z.size(); ++i)
      z[i] *= a;
    return *this;
  }

  // This a component wise operation
  template<typename REAL1>
  RVectorT &operator*=(const RVectorT<REAL1> &v) {
    if constexpr (DebugLevel > 0) {
      if (z.size() != v.size()) THROW("Size does not fit!")
    }
    for (int i = 0; i < z.size(); ++i)
      z[i] *= v[i];
    return *this;
  }

  // This a component wise operation
  template<typename REAL1>
  RVectorT &operator/=(const RVectorT<REAL1> &v) {
    if constexpr (DebugLevel > 0) {
      if (z.size() != v.size()) THROW("Size does not fit!")
    }
    for (int i = 0; i < z.size(); ++i)
      z[i] /= v[i];
    return *this;
  }

  template<typename REAL1>
  mpp_ba::enable_if_constructible<REAL, REAL1, RVectorT> &operator/=(const REAL1 &a) {
    if constexpr (std::is_same_v<REAL, double>) {
      if (a == 1.0) return *this;
    }
    for (int i = 0; i < z.size(); ++i)
      z[i] /= a;
    return *this;
  }

  void removeLast() { z.pop_back(); }

  const std::vector<REAL> &asVector() const { return z; }

  std::vector<REAL> &asVector() { return z; }

  int Dim() const { return z.size(); }

  int size() const { return z.size(); }

  void resize(int N) { z.resize(N); }

  REAL &operator[](int i) { return z[i]; }

  const REAL &operator[](int i) const { return z[i]; }

  const REAL *operator()() const { return z.data(); }

  REAL *operator()() { return z.data(); }

  auto begin() const { return z.begin(); }

  auto begin() { return z.begin(); }

  auto rbegin() const { return z.rbegin(); }

  auto end() const { return z.end(); }

  auto end() { return z.end(); }

  auto rend() const { return z.rend(); }

  auto front() const { return z.front(); }

  REAL &back() { return z.back(); }

  const REAL &back() const { return z.back(); }

  auto push_back(const REAL &a) { return z.push_back(a); }

  void clear() { z.clear(); }

  template<typename REAL1>
  auto append(const RVectorT<REAL1> &v) {
    return z.insert(z.end(), v.begin(), v.end());
  }

  template<typename REAL1>
  void insert(const RVectorT<REAL1> &v, int start) {
    if (start + v.size() > z.size()) { THROW("Cannot insert row: size of v too large") }
    for (int j = 0; j < v.size(); ++j)
      z[start + j] = v[j];
  }

  const std::vector<REAL> &Data() const { return z; }

  REAL normSqr() const;

  REAL norm() const;

  REAL normAccumulated(int commSplit = 0) const;

  REAL normOnes() const;

  REAL Mean() const;

  REAL Variance() const;

  REAL Min() const;

  REAL MinAccumulated(int commSplit = 0) const;

  REAL Max() const;

  REAL MaxAccumulated(int commSplit = 0) const;

  void Accumulate(int commSplit = 0);

  void CleanZero();

  Saver &save(Saver &saver) const;

  Loader &load(Loader &loader);

  template<typename REAL1>
  void Multiply(const REAL &b, const RVectorT<REAL1> &v) {
    if constexpr (DebugLevel > 0) {
      if (z.size() != v.size()) THROW("Size does not fit!")
    }
    for (int i = 0; i < z.size(); ++i)
      z[i] = b * v[i];
  }

  template<typename REAL1>
  void MultiplyPlus(const REAL &b, const RVectorT<REAL1> &v) {
    if constexpr (DebugLevel > 0) {
      if (z.size() != v.size()) THROW("Size does not fit!")
    }
    for (int i = 0; i < z.size(); ++i)
      z[i] += b * v[i];
  }

  template<typename REAL1>
  void MultiplyMinus(const REAL &b, const RVectorT<REAL1> &v) {
    if constexpr (DebugLevel > 0) {
      if (z.size() != v.size()) THROW("Size does not fit!")
    }
    for (int i = 0; i < z.size(); ++i)
      z[i] -= b * v[i];
  }
};

template<typename REAL>
bool operator==(const RVectorT<REAL> &x, const RVectorT<REAL> &y) {
  if (x.size() != y.size()) return false;
  for (int i = 0; i < x.size(); ++i)
    if (!mpp_ba::isNear(x[i], y[i])) return false;
  return true;
}

template<typename REAL>
bool operator!=(const RVectorT<REAL> &x, const RVectorT<REAL> &y) {
  return !(x == y);
}

template<typename REAL>
inline RVectorT<REAL> operator-(const RVectorT<REAL> &v) {
  RVectorT<REAL> z(v);
  return z *= -1.0;
}

template<typename REAL>
inline RVectorT<REAL> operator+(const RVectorT<REAL> &v, const RVectorT<REAL> &w) {
  RVectorT<REAL> z(v);
  return z += w;
}

template<typename REAL>
inline RVectorT<REAL> operator-(const RVectorT<REAL> &v, const RVectorT<REAL> &w) {
  RVectorT<REAL> z(v);
  return z -= w;
}

template<typename REAL, typename REAL1>
inline mpp_ba::enable_if_constructible<REAL, REAL1, RVectorT<REAL>>
operator*(const REAL1 &b, const RVectorT<REAL> &v) {
  RVectorT<REAL> z(v);
  return z *= b;
}

template<typename REAL, typename REAL1>
inline mpp_ba::enable_if_constructible<REAL, REAL1, RVectorT<REAL>>
operator*(const RVectorT<REAL> &v, const REAL1 &b) {
  RVectorT<REAL> z(v);
  return z *= b;
}

template<typename REAL, typename REAL1>
inline mpp_ba::enable_if_constructible<REAL, REAL1, RVectorT<REAL>>
operator/(const RVectorT<REAL> &v, const REAL1 &b) {
  RVectorT<REAL> z(v);
  return z /= b;
}

template<typename REAL>
REAL operator*(const RVectorT<REAL> &v, const RVectorT<REAL> &w) {
  if constexpr (DebugLevel > 0) {
    if (v.size() != w.size()) THROW("Size does not fit!")
  }
  REAL sum{};
  for (int i = 0; i < v.size(); ++i)
    sum += v[i] * w[i];
  return sum;
}

template<typename REAL>
inline REAL normSqr(const RVectorT<REAL> &v) {
  return v.normSqr();
}

template<typename REAL>
inline REAL norm(const RVectorT<REAL> &v) {
  return v.norm();
}

template<typename REAL>
inline Saver &operator<<(Saver &saver, const RVectorT<REAL> &v) {
  return v.save(saver);
}

template<typename REAL>
inline Loader &operator>>(Loader &loader, RVectorT<REAL> &v) {
  return v.load(loader);
}

template<typename REAL>
std::ostream &operator<<(std::ostream &os, const RVectorT<REAL> &v) {
  if constexpr (std::is_same_v<REAL, double>) os << beginD;
  for (int i = 0; i < v.size(); ++i) {
    if (i == 0) os << v[i];
    else {
      os << "\n" << v[i];
      if constexpr (std::is_same_v<REAL, double>) os << beginD;
    }
  }
  if constexpr (std::is_same_v<REAL, double>) os << endD;
  return os;
}

using RVector = RVectorT<double>;

/// Data vector for fem vectors and matrices
using BasicVector = RVectorT<double>;

#ifdef BUILD_IA

using IARVector = RVectorT<IAInterval>;

RVector mid(const IARVector &);

RVector sup(const IARVector &);

RVector inf(const IARVector &);

inline IARVector operator+(const IARVector &a, const RVector &b) {
  IARVector w(a);
  return w += b;
}

inline IARVector operator+(const RVector &a, const IARVector &b) {
  IARVector w(a);
  return w += b;
}

inline IARVector operator-(const IARVector &a, const RVector &b) {
  IARVector w(a);
  return w -= b;
}

inline IARVector operator-(const RVector &a, const IARVector &b) {
  IARVector w(a);
  return w -= b;
}

inline IAInterval operator*(const IARVector &a, const RVector &b) {
  if (a.size() != b.size()) THROW("Size does not fit!")
  IAInterval sum;
  for (int i = 0; i < a.size(); ++i)
    sum += a[i] * b[i];
  return sum;
}

inline IAInterval operator*(const RVector &a, const IARVector &b) {
  if (a.size() != b.size()) THROW("Size does not fit!")
  IAInterval sum;
  for (int i = 0; i < a.size(); ++i)
    sum += a[i] * b[i];
  return sum;
}

#endif // BUILD_IA

#endif