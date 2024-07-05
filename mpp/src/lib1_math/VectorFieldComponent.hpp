#ifndef VECTORFIELDCOMPONENT_HPP
#define VECTORFIELDCOMPONENT_HPP

#include "Point.hpp"

template<typename T = double>
class VectorFieldComponentT {
  T value{};
  int component{};
public:
  constexpr VectorFieldComponentT() {}

  constexpr VectorFieldComponentT(int k, const T &value) : value(value), component(k) {}

  template<typename TT>
  constexpr VectorFieldComponentT(const VectorFieldComponentT<TT> &y) {
    component = y.Component();
    value = y.Value();
  }

  constexpr VectorFieldComponentT &operator=(const T &a) {
    value = a;
    return *this;
  }

  constexpr const T &Value() const { return value; }

  constexpr int Component() const { return component; }

  constexpr T operator[](int k) const { return (k == component) * value; }

  template<typename TT>
  constexpr VectorFieldComponentT &operator=(const VectorFieldComponentT<TT> &y) {
    value = y.Value();
    component = y.Component();
    return *this;
  }

  template<typename TT>
  constexpr VectorFieldComponentT &operator*=(const TT &a) {
    value *= a;
    return *this;
  }

  template<typename TT>
  constexpr VectorFieldComponentT &operator/=(const TT &a) {
    value /= a;
    return *this;
  }

  constexpr bool operator==(const VectorFieldComponentT<T> &other) const {
    if (component != other.component) return false;
    if (mpp_geometry::isNear(value, other.value)) return true;
    return false;
  }

  constexpr bool operator!=(const VectorFieldComponentT<T> &other) const {
    return !((*this) == other);
  }

  Saver &save(Saver &saver) const { return saver << component << value; }

  Loader &load(Loader &loader) { return loader >> component >> value; }
};

template<typename T>
inline constexpr VectorFieldComponentT<T> operator-(const VectorFieldComponentT<T> &x) {
  VectorFieldComponentT<T> z(x);
  return z *= T(-1.0);
}

template<typename T>
constexpr T operator*(const VectorFieldComponentT<T> &x, const VectorFieldComponentT<T> &y) {
  return (x.Component() == y.Component()) * (x.Value() * y.Value());
}

template<typename T, typename TT>
inline constexpr std::enable_if_t<std::is_scalar<TT>::value, VectorFieldComponentT<T>>
operator*(const TT &b, const VectorFieldComponentT<T> &x) {
  VectorFieldComponentT<T> z(x);
  return z *= b;
}

template<typename T, typename TT>
inline constexpr std::enable_if_t<std::is_scalar<TT>::value, VectorFieldComponentT<T>>
operator*(const VectorFieldComponentT<T> &x, const TT &b) {
  VectorFieldComponentT<T> z(x);
  return z *= b;
}

template<typename T, typename TT>
inline constexpr std::enable_if_t<std::is_scalar<TT>::value, VectorFieldComponentT<T>>
operator/(const VectorFieldComponentT<T> &x, const TT &b) {
  VectorFieldComponentT<T> z(x);
  return z /= b;
}

template<typename T>
inline constexpr T norm(const VectorFieldComponentT<T> &x) {
  return abs(x.Value());
}

template<typename T>
inline Saver &operator<<(Saver &saver, const VectorFieldComponentT<T> &x) {
  return x.save(saver);
}

template<typename T>
inline Loader &operator>>(Loader &loader, VectorFieldComponentT<T> &x) {
  return x.load(loader);
}

template<typename T = double>
using VelocityComponentT = VectorFieldComponentT<T>;

using VectorFieldComponent = VectorFieldComponentT<double>;

using VelocityComponent = VectorFieldComponent;

#ifdef BUILD_IA

using IAVectorFieldComponent = VectorFieldComponentT<IAInterval>;

using IAVelocityComponent = IAVectorFieldComponent;

VelocityComponentT<double> mid(const VelocityComponentT<IAInterval> &IAx);

VelocityComponentT<double> sup(const VelocityComponentT<IAInterval> &IAx);

VelocityComponentT<double> inf(const VelocityComponentT<IAInterval> &IAx);

inline VelocityComponentT<IAInterval> operator*(const VelocityComponentT<double> &x,
                                                const IAInterval &b) {
  VelocityComponentT<IAInterval> IAx(x);
  return IAx *= b;
}

inline VelocityComponentT<IAInterval> operator*(const IAInterval &b,
                                                const VelocityComponentT<double> &x) {
  VelocityComponentT<IAInterval> IAx(x);
  return IAx *= b;
}

inline VelocityComponentT<IAInterval> operator/(const VelocityComponentT<double> &x,
                                                const IAInterval &b) {
  VelocityComponentT<IAInterval> IAx(x);
  return IAx /= b;
}

inline IAInterval operator*(const VelocityComponentT<IAInterval> &x,
                            const VelocityComponentT<double> &y) {
  VelocityComponentT<IAInterval> IAy(y);
  return x * IAy;
}

inline IAInterval operator*(const VelocityComponentT<double> &x,
                            const VelocityComponentT<IAInterval> &y) {
  VelocityComponentT<IAInterval> IAx(x);
  return IAx * y;
}

#endif // BUILD_IA

#endif // VECTORFIELDCOMPONENT_HPP
