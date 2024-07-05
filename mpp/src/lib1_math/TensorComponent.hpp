#ifndef TENSORCOMPONENT_HPP
#define TENSORCOMPONENT_HPP
#include <type_traits>
#include <utility>
#include "SaveLoad.hpp"
#include "VectorField.hpp"
#include "VectorFieldComponent.hpp"

template<typename T = double>
class TensorComponentT {
  T value{};
  std::pair<int, int> component{};
public:
  constexpr TensorComponentT() {}

  constexpr TensorComponentT(int k, int l, const T &value) : value(value), component(k, l) {}

  template<typename TT>
  constexpr TensorComponentT(const TensorComponentT<TT> &y) :
      value(y.Value()), component(y.Component()) {}

  constexpr TensorComponentT &operator=(const T &a) {
    value = a;
    return *this;
  }

  constexpr const T &Value() const { return value; }

  constexpr int Row() const { return component.first; }

  constexpr int Column() const { return component.second; }

  constexpr T Entry(int k, int l) const {
    return (k == component.first && l == component.second) * value;
  }

  template<typename TT>
  constexpr TensorComponentT &operator=(const TensorComponentT<TT> &y) {
    value = y.Value();
    component = y.Component();
    return *this;
  }

  template<typename TT>
  constexpr TensorComponentT &operator*=(const TT &a) {
    value *= a;
    return *this;
  }

  template<typename TT>
  constexpr TensorComponentT &operator/=(const TT &a) {
    value /= a;
    return *this;
  }

  template<int dim, int tDim>
  constexpr inline VectorFieldComponentT<T> multiply(const VectorFieldT<T, dim, tDim> &v) const {
    return VectorFieldComponentT<T>(Row(), value * v[Column()]);
  }

  constexpr inline VectorFieldComponentT<T> multiply(const VectorFieldComponentT<T> &v) const {
    return VectorFieldComponentT<T>(Row(), (v.Component() == Column()) * v.Value() * value);
  }

  template<int dim, int tDim>
  constexpr inline VectorFieldComponentT<T>
  applyTransposed(const VectorFieldT<T, dim, tDim> &v) const {
    return VectorFieldComponentT<T>(Column(), value * v[Row()]);
  }

  constexpr inline VectorFieldComponentT<T>
  applyTransposed(const VectorFieldComponentT<T> &v) const {
    return VectorFieldComponentT<T>(Column(), (Row() == v.Component()) * value * v.Value());
  }

  constexpr bool operator==(const TensorComponentT<T> &other) const {
    if (component.first != other.component.first || component.second != other.component.second)
      return false;
    if (mpp_geometry::isNear(value, other.value)) return true;
    return false;
  }

  constexpr bool operator!=(const TensorComponentT<T> &other) const { return !((*this) == other); }

  Saver &save(Saver &saver) const { return saver << component.first << component.second << value; }

  Loader &load(Loader &loader) { return loader >> component.first >> component.second >> value; }
};

template<typename T>
inline constexpr TensorComponentT<T> operator-(const TensorComponentT<T> &x) {
  TensorComponentT<T> z(x);
  return z *= T(-1.0);
}

template<typename T, typename TT>
inline constexpr std::enable_if_t<std::is_scalar<TT>::value, TensorComponentT<T>>
operator*(const TT &b, const TensorComponentT<T> &x) {
  TensorComponentT<T> z(x);
  return z *= b;
}

template<typename T, typename TT>
inline constexpr std::enable_if_t<std::is_scalar<TT>::value, TensorComponentT<T>>
operator*(const TensorComponentT<T> &x, const TT &b) {
  TensorComponentT<T> z(x);
  return z *= b;
}

template<typename T, typename TT>
inline constexpr std::enable_if_t<std::is_scalar<TT>::value, TensorComponentT<T>>
operator/(const TensorComponentT<T> &x, const TT &b) {
  TensorComponentT<T> z(x);
  return z /= b;
}

template<typename T, int dim>
constexpr inline T trace(const TensorComponentT<T> &S) {
  return S.trace();
}

template<typename T, int sDim, int tDim>
constexpr inline VectorFieldComponentT<T> operator*(const TensorComponentT<T> &S,
                                                    const VectorFieldT<T, sDim, tDim> &v) {
  return S.multiply(v);
}

template<typename T>
constexpr inline VectorFieldComponentT<T> operator*(const TensorComponentT<T> &S,
                                                    const VectorFieldComponentT<T> &v) {
  return S.multiply(v);
}

template<typename T, int sDim, int tDim>
constexpr inline VectorFieldComponentT<T> operator*(const VectorFieldT<T, sDim, tDim> &v,
                                                    const TensorComponentT<T> &S) {
  return S.applyTransposed(v);
}

template<typename T>
constexpr inline VectorFieldComponentT<T> operator*(const VectorFieldComponentT<T> &v,
                                                    const TensorComponentT<T> &S) {
  return S.applyTransposed(v);
}

template<typename T>
constexpr inline TensorComponentT<T> operator*(const TensorComponentT<T> &R,
                                               const TensorComponentT<T> &S) {
  return TensorComponentT<T>(R.Row(), S.Column(), (R.Column() == S.Row()) * R.Value() * S.Value());
}

template<typename T>
inline constexpr T norm(const TensorComponentT<T> &x) {
  return abs(x.Value());
}

template<typename T>
inline constexpr T trace(const TensorComponentT<T> &x) {
  return (x.Row() == x.Column()) * x.Value();
}

template<typename T>
inline constexpr T Frobenius(const TensorComponentT<T> &x, const TensorComponentT<T> &y) {
  return (x.Row() == y.Row() && x.Column() == y.Column()) * x.Value() * y.Value();
}

template<typename T>
inline constexpr TensorComponentT<T> Product(const VectorFieldComponentT<T> &x,
                                             const VectorFieldComponentT<T> &y) {
  return TensorComponentT<T>(x.Component(), y.Component(), x.Value() * y.Value());
}

template<typename T>
inline Saver &operator<<(Saver &saver, const TensorComponentT<T> &x) {
  return x.save(saver);
}

template<typename T>
inline Loader &operator>>(Loader &loader, TensorComponentT<T> &x) {
  return x.load(loader);
}

using TensorComponent = TensorComponentT<>;


#endif // TENSORCOMPONENT_HPP
