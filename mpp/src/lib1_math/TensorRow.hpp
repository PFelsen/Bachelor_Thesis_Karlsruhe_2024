#ifndef TENSORROW_HPP
#define TENSORROW_HPP

#include "VectorField.hpp"

template<typename T = double, int dim = SpaceDimension>
class TensorRowT {
protected:
  /// Stores value, i.e., a single row
  VectorFieldT<T, dim> value{};
  int component{};
public:
  template<typename TT>
  constexpr TensorRowT(int i, const VectorFieldT<TT, dim> &U) {
    component = i;
    value = U;
  }

  constexpr TensorRowT() {}

  template<typename TT>
  constexpr TensorRowT(const VectorFieldT<TT, dim> &U, int i) : TensorRowT(i, U) {}

  template<typename TT>
  constexpr TensorRowT(const TensorRowT<TT, dim> &S) {
    component = S.Component();
    value = S.Value();
  }

  template<typename T1, typename T2>
  constexpr TensorRowT(const TensorRowT<T1, dim> &R, const TensorRowT<T2, dim> &S) {
    component = R.Component();
    if constexpr (dim == 1) {
      value[0] = R.Value() * S.Value();
    } else if constexpr (dim == 2) {
      value[0] = R.Value()[S.Component()] * S.Value()[0];
      value[1] = R.Value()[S.Component()] * S.Value()[1];
    } else if constexpr (dim == 3) {
      value[0] = R.Value()[S.Component()] * S.Value()[0];
      value[1] = R.Value()[S.Component()] * S.Value()[1];
      value[2] = R.Value()[S.Component()] * S.Value()[2];
    } else {
      THROW("Not implemented for dim higher than 3")
    }
  }

  constexpr int Dim() const { return dim; }

  constexpr const VectorFieldT<T, dim> &Value() const { return value; }

  constexpr int Component() const { return component; }

  constexpr VectorFieldT<T, dim> operator[](int k) const { return (k == component) * value; }

  template<typename TT>
  constexpr TensorRowT &operator=(const TensorRowT<TT, dim> &S) {
    component = S.Component();
    value = S.Value();
    return *this;
  }

  template<typename TT>
  constexpr TensorRowT &operator*=(const TT &a) {
    value *= a;
    return *this;
  }

  template<typename TT>
  constexpr TensorRowT &operator/=(const TT &a) {
    value /= a;
    return *this;
  }

  template<int tDim>
  constexpr inline VectorFieldComponentT<T> multiply(const VectorFieldT<T, dim, tDim> &v) const {
    if constexpr (dim == 1) return VectorFieldComponentT<T>(component, value[0] * v[0]);
    if constexpr (dim == 2)
      return VectorFieldComponentT<T>(component, value[0] * v[0] + value[1] * v[1]);
    if constexpr (dim == 3)
      return VectorFieldComponentT<T>(component,
                                      value[0] * v[0] + value[1] * v[1] + value[2] * v[2]);
    THROW("not implemented")
  }

  constexpr inline VectorFieldComponentT<T> multiply(const VectorFieldComponentT<T> &v) const {
    return VectorFieldComponentT<T>(component, v.Value() * value[v.Component()]);
  }

  constexpr inline VectorFieldT<T, dim, 0> applyTransposed(const VectorFieldT<T, dim, 0> &v) const {
    return value * v[component];
  }

  constexpr inline VectorFieldT<T, dim, 1> applyTransposed(const VectorFieldT<T, dim, 1> &v) const {
    return value * v[component];
  }

  constexpr inline VectorFieldT<T, dim> applyTransposed(const VectorFieldComponentT<T> &v) const {
    return ((component == v.Component()) * v.Value()) * value;
  }

  constexpr T trace() const { return value[component]; }

  constexpr T norm() const { return sqrt(value * value); }

  constexpr T maxnorm() const {
    // maxnorm = max_i \sum_{j=1}^m |t_ij|
    T sum{};
    for (int j = 0; j < dim; ++j)
      sum += abs(value[j]);
    return sum;
  }

  friend bool operator==(const TensorRowT<T, dim> &R, const TensorRowT<T, dim> &S) {
    return (R.component == S.component) && (R.value == S.value);
  }

  friend bool operator!=(const TensorRowT<T, dim> &R, const TensorRowT<T, dim> &S) {
    return !(R == S);
  }

  Saver &save(Saver &saver) const { return saver << dim << component << value; }

  Loader &load(Loader &loader) {
    int loadedDim;
    loader >> loadedDim;
    if (loadedDim != dim) THROW("Tensor Dimension does not fit")
    return loader >> component >> value;
  }
};

template<typename T, int dim>
constexpr inline TensorRowT<T, dim> operator-(const TensorRowT<T, dim> &S) {
  TensorRowT<T, dim> bS(S);
  return bS *= T(-1.0);
}

template<typename T, typename TT, int dim>
constexpr inline TensorRowT<T, dim> operator*(const TT &b, const TensorRowT<T, dim> &S) {
  TensorRowT<T, dim> bS(S);
  return bS *= b;
}

template<typename T, typename TT, int dim>
constexpr inline TensorRowT<T, dim> operator*(const TensorRowT<T, dim> &S, const TT &b) {
  TensorRowT<T, dim> bS(S);
  return bS *= b;
}

template<typename T, typename TT, int dim>
constexpr inline TensorRowT<T, dim> operator/(const TensorRowT<T, dim> &S, const TT &b) {
  TensorRowT<T, dim> bS(S);
  return bS /= b;
}

template<typename T, int dim>
constexpr inline T trace(const TensorRowT<T, dim> &S) {
  return S.trace();
}

template<typename T, int sDim, int tDim>
constexpr inline VectorFieldComponentT<T> operator*(const TensorRowT<T, sDim> &S,
                                                    const VectorFieldT<T, sDim, tDim> &v) {
  return S.multiply(v);
}

template<typename T, int sDim>
constexpr inline VectorFieldComponentT<T> operator*(const TensorRowT<T, sDim> &S,
                                                    const VectorFieldComponentT<T> &v) {
  return S.multiply(v);
}

template<typename T, int sDim>
constexpr inline VectorFieldT<T, sDim, 0> operator*(const VectorFieldT<T, sDim, 0> &v,
                                                    const TensorRowT<T, sDim> &S) {
  return S.applyTransposed(v);
}

template<typename T, int sDim>
constexpr inline VectorFieldT<T, sDim, 1> operator*(const VectorFieldT<T, sDim, 1> &v,
                                                    const TensorRowT<T, sDim> &S) {
  return S.applyTransposed(v);
}

template<typename T, int sDim>
constexpr inline VectorFieldT<T, sDim> operator*(const VectorFieldComponentT<T> &v,
                                                 const TensorRowT<T, sDim> &S) {
  return S.applyTransposed(v);
}

template<typename T, int dim>
constexpr inline TensorRowT<T, dim> operator*(const TensorRowT<T, dim> &R,
                                              const TensorRowT<T, dim> &S) {
  return TensorRowT<T, dim>(R, S);
}

template<typename T, int dim>
T Frobenius(const TensorRowT<T, dim> &R, const TensorRowT<T, dim> &S) {
  return (R.Component() == S.Component()) * (R.Value() * S.Value());
}

template<typename T, int dim>
TensorRowT<T, dim> Product(const VectorFieldComponentT<T> &R, const VectorFieldT<T, dim> &S) {
  return TensorRowT(R.Component(), R.Value() * S);
}

template<typename T, int dim>
constexpr inline T norm(const TensorRowT<T, dim> &S) {
  return S.norm();
}

template<int dim>
constexpr inline double maxnorm(const TensorRowT<double, dim> &S) {
  return S.maxnorm();
}

template<typename T, int dim>
constexpr inline Saver &operator<<(Saver &saver, const TensorRowT<T, dim> &S) {
  return S.save(saver);
}

template<typename T, int dim>
constexpr inline Loader &operator>>(Loader &loader, TensorRowT<T, dim> &S) {
  return S.load(loader);
}

template<typename T = double, int dim = SpaceDimension>
using VelocityGradientRowT = TensorRowT<T, dim>;

using TensorRow = TensorRowT<>;

using VelocityGradientRow = VelocityGradientRowT<>;

#ifdef BUILD_IA

using IATensorRow = TensorRowT<IAInterval, SpaceDimension>;

using IAVelocityGradientRow = VelocityGradientRowT<IAInterval, SpaceDimension>;

template<int sDim>
TensorRowT<double, sDim> mid(const TensorRowT<IAInterval, sDim> &IA_T) {
  return TensorRowT<double, sDim>(IA_T.Component(), mid(IA_T.Value()));
}

template<int sDim>
TensorRowT<double, sDim> sup(const TensorRowT<IAInterval, sDim> &IA_T) {
  return TensorRowT<double, sDim>(IA_T.Component(), sup(IA_T.Value()));
}

template<int sDim>
TensorRowT<double, sDim> inf(const TensorRowT<IAInterval, sDim> &IA_T) {
  return TensorRowT<double, sDim>(IA_T.Component(), inf(IA_T.Value()));
}

template<int sDim>
inline TensorRowT<IAInterval, sDim> operator*(const TensorRowT<double, sDim> &x,
                                              const IAInterval &b) {
  TensorRowT<IAInterval, sDim> IAx(x);
  return IAx *= b;
}

template<int sDim>
inline TensorRowT<IAInterval, sDim> operator*(const IAInterval &b,
                                              const TensorRowT<double, sDim> &x) {
  TensorRowT<IAInterval, sDim> IAx(x);
  return IAx *= b;
}

template<int sDim>
inline TensorRowT<IAInterval, sDim> operator/(const TensorRowT<double, sDim> &x,
                                              const IAInterval &b) {
  TensorRowT<IAInterval, sDim> IAx(x);
  return IAx /= b;
}

#endif // BUILD_IA

#endif // TENSORROW_HPP
