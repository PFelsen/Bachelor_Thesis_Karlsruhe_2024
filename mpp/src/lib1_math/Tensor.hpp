#ifndef _TENSOR_H_
#define _TENSOR_H_

#include "GlobalDefinitions.hpp"
#include "SymTensor.hpp"
#include "TensorComponent.hpp"
#include "TensorRow.hpp"
#include "VectorField.hpp"

template<typename T = double, int dim = SpaceDimension>
class TensorT {
protected:
  /// Stores values of Tensor in Rows
  std::array<VectorFieldT<T, dim>, dim> t{};
public:
  constexpr TensorT() {
    for (int n = 0; n < dim; ++n)
      t[n] = T{};
  }

  constexpr TensorT(const T &a) {
    for (int n = 0; n < dim; ++n)
      t[n] = a;
  }

  constexpr TensorT(int i, int j, T a = T(1.0)) {
    for (int n = 0; n < dim; ++n)
      t[n] = T{};
    t[i][j] = a;
  }

  template<typename TT>
  constexpr TensorT(int i, const VectorFieldT<TT, dim> &U) {
    for (int n = 0; n < dim; ++n)
      t[n] = T{};
    t[i] = U;
  }

  template<typename TT>
  constexpr TensorT(int i, const VectorFieldComponentT<TT> &U) {
    for (int n = 0; n < dim; ++n)
      t[n] = T{};
    t[i][U.Component()] = U.Value();
  }

  template<typename TT>
  constexpr TensorT(const VectorFieldT<TT, dim> &U, int i) : TensorT(i, U) {}

  template<typename TT>
  constexpr TensorT(const TensorT<TT, dim> &S) {
    for (int n = 0; n < dim; ++n)
      t[n] = S[n];
  }

  template<typename TT>
  constexpr TensorT(const SymTensorT<TT, dim> &S) {
    for (int n = 0; n < dim; ++n)
      for (int m = 0; m < dim; ++m)
        t[n][m] = S(n, m);
  }

  constexpr explicit TensorT(const TensorRowT<T, dim> &S) : TensorT(S.Component(), S.Value()) {}

  constexpr explicit TensorT(const TensorComponentT<T> &S) :
      TensorT(S.Row(), S.Column(), S.Value()) {}

  template<typename TT>
  constexpr TensorT(const VectorFieldT<TT, dim> &U,
                    const VectorFieldT<TT, dim> &V = VectorFieldT<TT, dim>(),
                    const VectorFieldT<TT, dim> &W = VectorFieldT<TT, dim>()) {

    if constexpr (dim > 0) {
      t[0] = U;
    } else {
      THROW("Tensor of Dimension 0 not specified")
    }
    if constexpr (dim > 1) { t[1] = V; }
    if constexpr (dim > 2) { t[2] = W; }
  }

  template<typename T1, typename T2>
  constexpr TensorT(const TensorT<T1, dim> &R, const TensorT<T2, dim> &S) : TensorT() {
    if constexpr (dim == 1) { t[0][0] = R[0][0] * S[0][0]; }
    if constexpr (dim == 2) {
      t[0][0] = R[0][0] * S[0][0] + R[0][1] * S[1][0];
      t[0][1] = R[0][0] * S[0][1] + R[0][1] * S[1][1];
      t[1][0] = R[1][0] * S[0][0] + R[1][1] * S[1][0];
      t[1][1] = R[1][0] * S[0][1] + R[1][1] * S[1][1];
    } else if constexpr (dim == 3) {
      t[0][0] = R[0][0] * S[0][0] + R[0][1] * S[1][0] + R[0][2] * S[2][0];
      t[0][1] = R[0][0] * S[0][1] + R[0][1] * S[1][1] + R[0][2] * S[2][1];
      t[0][2] = R[0][0] * S[0][2] + R[0][1] * S[1][2] + R[0][2] * S[2][2];
      t[1][0] = R[1][0] * S[0][0] + R[1][1] * S[1][0] + R[1][2] * S[2][0];
      t[1][1] = R[1][0] * S[0][1] + R[1][1] * S[1][1] + R[1][2] * S[2][1];
      t[1][2] = R[1][0] * S[0][2] + R[1][1] * S[1][2] + R[1][2] * S[2][2];
      t[2][0] = R[2][0] * S[0][0] + R[2][1] * S[1][0] + R[2][2] * S[2][0];
      t[2][1] = R[2][0] * S[0][1] + R[2][1] * S[1][1] + R[2][2] * S[2][1];
      t[2][2] = R[2][0] * S[0][2] + R[2][1] * S[1][2] + R[2][2] * S[2][2];
    }
  }

  constexpr TensorT(T a00, T a01, T a02, T a10, T a11, T a12, T a20, T a21, T a22) : TensorT() {
    if constexpr (dim > 0) { t[0][0] = a00; }
    if constexpr (dim > 1) {
      t[0][1] = a01;
      t[1][0] = a10;
      t[1][1] = a11;
    }

    if constexpr (dim > 2) {
      t[0][2] = a02;
      t[1][2] = a12;
      t[2][0] = a20;
      t[2][1] = a21;
      t[2][2] = a22;
    }

    if constexpr (dim > 3) { THROW("Constructor for dimension greater 3 not defined") }
  }

  constexpr TensorT(T a00, T a01, T a10, T a11) :
      TensorT(a00, a01, T{}, a10, a11, T{}, T{}, T{}, T(1.0)){};

  constexpr int Dim() const { return dim; }

  constexpr const VectorFieldT<T, dim> &operator[](int k) const { return t[k]; }

  constexpr VectorFieldT<T, dim> &operator[](int k) { return t[k]; }

  constexpr TensorT &operator=(const T &a) {
    for (int n = 0; n < dim; ++n)
      t[n] = a;
    return *this;
  }

  template<typename TT>
  constexpr TensorT &operator=(const TensorT<TT, dim> &S) {
    for (int n = 0; n < dim; ++n)
      t[n] = S[n];
    return *this;
  }

  template<typename TT>
  constexpr TensorT &operator=(const SymTensorT<TT, dim> &S) {
    for (int n = 0; n < dim; ++n)
      for (int m = 0; m < dim; ++m)
        t[n][m] = S(n, m);
    return *this;
  }

  template<typename TT>
  constexpr TensorT &operator=(const TensorRowT<TT, dim> &R) {
    for (int n = 0; n < dim; ++n)
      t[n] = 0.0;
    t[R.Component()] = R.Value();
    return *this;
  }

  template<typename TT>
  constexpr TensorT &operator=(const TensorComponentT<TT> &R) {
    t[R.Row()][R.Column()] = R.Value();
    return *this;
  }

  template<typename TT>
  constexpr TensorT &operator+=(const TensorT<TT, dim> &S) {
    for (int n = 0; n < dim; ++n)
      t[n] += S[n];
    return *this;
  }

  template<typename TT>
  constexpr TensorT &operator+=(const SymTensorT<TT, dim> &S) {
    for (int n = 0; n < dim; ++n)
      for (int m = 0; m < dim; ++m)
        t[n][m] += S(n, m);
    return *this;
  }

  template<typename TT>
  constexpr TensorT &operator+=(const TensorRowT<TT, dim> &S) {
    t[S.Component()] += S.Value();
    return *this;
  }

  template<typename TT>
  constexpr TensorT &operator+=(const TensorComponentT<TT> &C) {
    t[C.Row()][C.Column()] += C.Value();
    return *this;
  }

  template<typename TT>
  constexpr TensorT &operator-=(const TensorT<TT, dim> &S) {
    for (int n = 0; n < dim; ++n)
      t[n] -= S[n];
    return *this;
  }

  template<typename TT>
  constexpr TensorT &operator-=(const SymTensorT<TT, dim> &S) {
    for (int n = 0; n < dim; ++n)
      for (int m = 0; m < dim; ++m)
        t[n][m] -= S(n, m);
    return *this;
  }

  template<typename TT>
  constexpr TensorT &operator-=(const TensorRowT<TT, dim> &S) {
    t[S.Component()] -= S.Value();
    return *this;
  }

  template<typename TT>
  constexpr TensorT &operator-=(const TensorComponentT<TT> &C) {
    t[C.Row()][C.Column()] -= C.Value();
    return *this;
  }

  template<typename TT>
  constexpr TensorT &operator*=(const TT &a) {
    for (int n = 0; n < dim; ++n)
      t[n] *= a;
    return *this;
  }

  constexpr TensorT &operator*=(const TensorT<T, dim> &S) {
    if constexpr (dim == 1) { t[0][0] = t[0][0] * S[0][0]; }
    if constexpr (dim == 2) {
      t[0][0] = t[0][0] * S[0][0] + t[0][1] * S[1][0];
      t[0][1] = t[0][0] * S[0][1] + t[0][1] * S[1][1];
      t[1][0] = t[1][0] * S[0][0] + t[1][1] * S[1][0];
      t[1][1] = t[1][0] * S[0][1] + t[1][1] * S[1][1];
    } else if constexpr (dim == 3) {
      t[0][0] = t[0][0] * S[0][0] + t[0][1] * S[1][0] + t[0][2] * S[2][0];
      t[0][1] = t[0][0] * S[0][1] + t[0][1] * S[1][1] + t[0][2] * S[2][1];
      t[0][2] = t[0][0] * S[0][2] + t[0][1] * S[1][2] + t[0][2] * S[2][2];
      t[1][0] = t[1][0] * S[0][0] + t[1][1] * S[1][0] + t[1][2] * S[2][0];
      t[1][1] = t[1][0] * S[0][1] + t[1][1] * S[1][1] + t[1][2] * S[2][1];
      t[1][2] = t[1][0] * S[0][2] + t[1][1] * S[1][2] + t[1][2] * S[2][2];
      t[2][0] = t[2][0] * S[0][0] + t[2][1] * S[1][0] + t[2][2] * S[2][0];
      t[2][1] = t[2][0] * S[0][1] + t[2][1] * S[1][1] + t[2][2] * S[2][1];
      t[2][2] = t[2][0] * S[0][2] + t[2][1] * S[1][2] + t[2][2] * S[2][2];
    }

    return *this;
  }

  template<typename TT>
  constexpr TensorT &operator/=(const TT &a) {
    for (int n = 0; n < dim; ++n)
      t[n] /= a;
    return *this;
  }

  template<int tDim>
  constexpr inline VectorFieldT<T, dim, tDim> multiply(const VectorFieldT<T, dim, tDim> &v) const {
    if constexpr (dim == 1) return VectorFieldT<T, dim, tDim>(t[0][0] * v);
    if constexpr (dim == 2)
      return VectorFieldT<T, dim, tDim>(t[0][0] * v[0] + t[0][1] * v[1],
                                        t[1][0] * v[0] + t[1][1] * v[1]);
    if constexpr (dim == 3)
      return VectorFieldT<T, dim, tDim>(t[0][0] * v[0] + t[0][1] * v[1] + t[0][2] * v[2],
                                        t[1][0] * v[0] + t[1][1] * v[1] + t[1][2] * v[2],
                                        t[2][0] * v[0] + t[2][1] * v[1] + t[2][2] * v[2]);
    THROW("not implemented")
  }

  constexpr inline VectorFieldT<T, dim> multiply(const VectorFieldComponentT<T> &v) const {
    if constexpr (dim == 1) return VectorFieldT<T, dim>(t[0][0] * v.Value());
    if constexpr (dim == 2)
      return VectorFieldT<T, dim>(t[0][v.Component()] * v.Value(), t[1][v.Component()] * v.Value());
    if constexpr (dim == 3)
      return VectorFieldT<T, dim>(t[0][v.Component()] * v.Value(), t[1][v.Component()] * v.Value(),
                                  t[2][v.Component()] * v.Value());
    THROW("not implemented")
  }

  template<int tDim>
  constexpr inline PointT<T, dim, tDim> multiply(const PointT<T, dim, tDim> &v) const {
    if constexpr (dim == 1) return PointT<T, dim, tDim>(t[0][0] * v);
    if constexpr (dim == 2)
      return PointT<T, dim, tDim>(t[0][0] * v[0] + t[0][1] * v[1], t[1][0] * v[0] + t[1][1] * v[1]);
    if constexpr (dim == 3)
      return PointT<T, dim, tDim>(t[0][0] * v[0] + t[0][1] * v[1] + t[0][2] * v[2],
                                  t[1][0] * v[0] + t[1][1] * v[1] + t[1][2] * v[2],
                                  t[2][0] * v[0] + t[2][1] * v[1] + t[2][2] * v[2]);
    THROW("not implemented")
  }

  constexpr TensorT &transpose() {
    if constexpr (dim == 1) {
      return *this;
    } else if constexpr (dim == 2) {
      T tmp = t[1][0];
      t[1][0] = t[0][1];
      t[0][1] = tmp;
      return *this;
    } else if constexpr (dim == 3) {
      T tmp = t[1][0];
      t[1][0] = t[0][1];
      t[0][1] = tmp;
      tmp = t[2][0];
      t[2][0] = t[0][2];
      t[0][2] = tmp;
      tmp = t[2][1];
      t[2][1] = t[1][2];
      t[1][2] = tmp;
      return *this;
    } else THROW("Transpose not defined for dimension greater than 3!")
  }

  constexpr T trace() const {
    T sum{};
    for (int n = 0; n < dim; ++n)
      sum += t[n][n];
    return sum;
  }

  constexpr T det() const {
    if constexpr (dim == 1) {
      return t[0][0];
    } else if constexpr (dim == 2) {
      return t[0][0] * t[1][1] - t[1][0] * t[0][1];
    } else if constexpr (dim == 3) {
      return t[0][0] * (t[1][1] * t[2][2] - t[2][1] * t[1][2])
             + t[1][0] * (t[2][1] * t[0][2] - t[0][1] * t[2][2])
             + t[2][0] * (t[0][1] * t[1][2] - t[1][1] * t[0][2]);
    } else THROW("Determinant not defined for dimension greater than 3!")
  }

  constexpr TensorT &Invert() {
    if constexpr (dim == 1) {
      t[0][0] = 1.0 / t[0][0];
    } else if constexpr (dim == 2) {
      T Det = det();
      if (abs(Det) < Eps) THROW("Error in Tensor<2> Inversion: det<Eps;")
      T tmp = t[0][0];
      t[0][0] = t[1][1] / Det;
      t[1][1] = tmp / Det;
      t[0][1] = -t[0][1] / Det;
      t[1][0] = -t[1][0] / Det;
      return *this;
    } else if constexpr (dim == 3) {
      T Det = det();
      if (abs(Det) < Eps) THROW("Error in Tensor<3> Inversion: det<Eps;")
      TensorT<T, 3> R(*this);
      t[0][0] = (R[1][1] * R[2][2] - R[2][1] * R[1][2]) / Det;
      t[0][1] = (-R[0][1] * R[2][2] + R[2][1] * R[0][2]) / Det;
      t[0][2] = (R[0][1] * R[1][2] - R[1][1] * R[0][2]) / Det;
      t[1][0] = (-R[1][0] * R[2][2] + R[2][0] * R[1][2]) / Det;
      t[1][1] = (R[0][0] * R[2][2] - R[2][0] * R[0][2]) / Det;
      t[1][2] = (-R[0][0] * R[1][2] + R[1][0] * R[0][2]) / Det;
      t[2][0] = (R[1][0] * R[2][1] - R[2][0] * R[1][1]) / Det;
      t[2][1] = (-R[0][0] * R[2][1] + R[2][0] * R[0][1]) / Det;
      t[2][2] = (R[0][0] * R[1][1] - R[1][0] * R[0][1]) / Det;
    } else THROW("Invert not defined for dimension greater than 3!")
    return *this;
  }

  constexpr TensorT &transposeInvert() {
    if constexpr (dim == 1) {
      t[0][0] = 1.0 / t[0][0];
    } else if constexpr (dim == 2) {
      T Det = det();
      if (abs(Det) < Eps) THROW("Error in Tensor<2> Inversion: det<Eps;")
      T tmp = t[0][0];
      t[0][0] = t[1][1] / Det;
      t[1][1] = tmp / Det;
      tmp = t[0][1];
      t[0][1] = -t[1][0] / Det;
      t[1][0] = -tmp / Det;
    } else if constexpr (dim == 3) {
      T Det = det();
      if (abs(Det) < Eps) THROW("Error in Tensor<2> Inversion: det<Eps;")
      TensorT<T, 3> R(*this);
      t[0][0] = (R[1][1] * R[2][2] - R[2][1] * R[1][2]) / Det;
      t[1][0] = (-R[0][1] * R[2][2] + R[2][1] * R[0][2]) / Det;
      t[2][0] = (R[0][1] * R[1][2] - R[1][1] * R[0][2]) / Det;
      t[0][1] = (-R[1][0] * R[2][2] + R[2][0] * R[1][2]) / Det;
      t[1][1] = (R[0][0] * R[2][2] - R[2][0] * R[0][2]) / Det;
      t[2][1] = (-R[0][0] * R[1][2] + R[1][0] * R[0][2]) / Det;
      t[0][2] = (R[1][0] * R[2][1] - R[2][0] * R[1][1]) / Det;
      t[1][2] = (-R[0][0] * R[2][1] + R[2][0] * R[0][1]) / Det;
      t[2][2] = (R[0][0] * R[1][1] - R[1][0] * R[0][1]) / Det;
    } else THROW("TransposeInvert not defined for dimension greater than 3!")
    return *this;
  }

  constexpr TensorT &sym() {
    if constexpr (dim == 1) {
      return *this;
    } else if constexpr (dim == 2) {
      t[0][1] = t[1][0] = (t[1][0] + t[0][1]) / 2;
      return *this;
    } else if constexpr (dim == 3) {
      t[0][1] = t[1][0] = (t[1][0] + t[0][1]) / 2;
      t[0][2] = t[2][0] = (t[2][0] + t[0][2]) / 2;
      t[1][2] = t[2][1] = (t[1][2] + t[2][1]) / 2;
      return *this;
    } else THROW("Sym not defined for dimension greater than 3!")
  }

  constexpr TensorT &skew() {
    if constexpr (dim == 1) {
      t[0][0] = T{};
      return *this;
    } else if constexpr (dim == 2) {
      t[0][0] = t[1][1] = T{};
      t[0][1] = (t[0][1] - t[1][0]) / 2;
      t[1][0] = -t[0][1];
      return *this;
    } else if constexpr (dim == 3) {
      t[0][0] = t[1][1] = t[2][2] = T{};
      t[0][1] = (t[0][1] - t[1][0]) / 2;
      t[0][2] = (t[0][2] - t[2][0]) / 2;
      t[1][2] = (t[1][2] - t[2][1]) / 2;
      t[1][0] = -t[0][1];
      t[2][0] = -t[0][2];
      t[2][1] = -t[1][2];
      return *this;
    } else THROW("Skew not defined for dimension greater than 3!")
  }

  constexpr T norm() const {
    T sum{};
    for (int n = 0; n < dim; ++n)
      sum += t[n] * t[n];
    return sqrt(sum);
  }

  constexpr T maxnorm() const {
    // maxnorm = max_i \sum_{j=1}^m |t_ij|
    T n{};
    for (int i = 0; i < dim; ++i) {
      T sum{};
      for (int j = 0; j < dim; ++j)
        sum += abs(t[i][j]);
      n = max(n, sum);
    }
    return n;
  }

  const VectorFieldT<T, dim> &Row(int i) const { return t[i]; }

  VectorFieldT<T, dim> Col(int i) const {
    if constexpr (dim == 1) return VectorFieldT<T, dim>(t[0][i]);
    if constexpr (dim == 2) return VectorFieldT<T, dim>(t[0][i], t[1][i]);
    if constexpr (dim == 3) return VectorFieldT<T, dim>(t[0][i], t[1][i], t[2][i]);
    else THROW("Not definet for dimensions higher than 3.")
  }

  template<int tDim>
  VectorFieldT<T, dim, tDim> applyTransposed(const VectorFieldT<T, dim, tDim> &z) const {
    if constexpr (dim == 1) {
      return VectorFieldT<T, dim>(t[0][0] * z[0]);
    } else if constexpr (dim == 2) {
      return VectorFieldT<T, dim>(t[0][0] * z[0] + t[1][0] * z[1], t[0][1] * z[0] + t[1][1] * z[1]);
    } else if constexpr (dim == 3) {
      return VectorFieldT<T, dim>(t[0][0] * z[0] + t[1][0] * z[1] + t[2][0] * z[2],
                                  t[0][1] * z[0] + t[1][1] * z[1] + t[2][1] * z[2],
                                  t[0][2] * z[0] + t[1][2] * z[1] + t[2][2] * z[2]);
    } else {
      THROW("Not definet for dimensions higher than 3.")
    }
  }

  VectorFieldT<T, dim> applyTransposed(const VectorFieldComponentT<T> &z) const {
    return t[z.Component()] * z.Value();
  }

  template<int tDim>
  PointT<T, dim, tDim> applyTransposed(const PointT<T, dim, tDim> &z) const {
    if constexpr (dim == 1) {
      return PointT<T, dim>(t[0][0] * z[0]);
    } else if constexpr (dim == 2) {
      return PointT<T, dim>(t[0][0] * z[0] + t[1][0] * z[1], t[0][1] * z[0] + t[1][1] * z[1]);
    } else if constexpr (dim == 3) {
      return PointT<T, dim>(t[0][0] * z[0] + t[1][0] * z[1] + t[2][0] * z[2],
                            t[0][1] * z[0] + t[1][1] * z[1] + t[2][1] * z[2],
                            t[0][2] * z[0] + t[1][2] * z[1] + t[2][2] * z[2]);
    } else {
      THROW("Not definet for dimensions higher than 3.")
    }
  }

  friend bool operator==(const TensorT<T, dim> &R, const TensorT<T, dim> &S) {
    bool equal = true;
    for (int i = 0; i < dim; ++i)
      equal = equal && R[i] == S[i];
    return equal;
  }

  friend bool operator!=(const TensorT<T, dim> &R, const TensorT<T, dim> &S) { return !(R == S); }

  template<typename OSTREAM>
  constexpr OSTREAM &print(OSTREAM &os) const {
    os << t[0];
    for (int n = 1; n < dim; ++n)
      os << "\n" << t[n];
    return os;
  }

  Saver &save(Saver &saver) const {
    saver << dim;
    for (int n = 0; n < dim; ++n)
      saver << t[n];
    return saver;
  }

  Loader &load(Loader &loader) {
    int loadedDim;
    loader >> loadedDim;
    if (loadedDim != dim) THROW("Tensor Dimension does not fit")
    for (int n = 0; n < dim; ++n)
      loader >> t[n];
    return loader;
  }
};

template<typename T, int dim>
constexpr inline TensorT<T, dim> operator+(const TensorT<T, dim> &R, const TensorT<T, dim> &S) {
  TensorT<T, dim> RS(R);
  return RS += S;
}

template<typename T, int dim>
constexpr inline TensorT<T, dim> operator+(const SymTensorT<T, dim> &R, const TensorT<T, dim> &S) {
  TensorT<T, dim> RS(R);
  return RS += S;
}

template<typename T, int dim>
constexpr inline TensorT<T, dim> operator+(const TensorT<T, dim> &R, const SymTensorT<T, dim> &S) {
  TensorT<T, dim> RS(R);
  return RS += S;
}

template<typename T, int dim>
constexpr inline TensorT<T, dim> operator+(const TensorT<T, dim> &R, const TensorRowT<T, dim> &S) {
  TensorT<T, dim> RS(R);
  return RS += S;
}

template<typename T, int dim>
constexpr inline TensorT<T, dim> operator+(const TensorRowT<T, dim> &R, const TensorT<T, dim> &S) {
  TensorT<T, dim> RS(S);
  return RS += R;
}

template<typename T, int dim>
constexpr inline TensorT<T, dim> operator+(const TensorT<T, dim> &R, const TensorComponentT<T> &S) {
  TensorT<T, dim> RS(R);
  return RS += S;
}

template<typename T, int dim>
constexpr inline TensorT<T, dim> operator+(const TensorComponentT<T> &R, const TensorT<T, dim> &S) {
  TensorT<T, dim> RS(S);
  return RS += R;
}

template<typename T, int dim>
constexpr inline TensorT<T, dim> operator-(const TensorT<T, dim> &R, const TensorT<T, dim> &S) {
  TensorT<T, dim> RS(R);
  return RS -= S;
}

template<typename T, int dim>
constexpr inline TensorT<T, dim> operator-(const SymTensorT<T, dim> &R, const TensorT<T, dim> &S) {
  TensorT<T, dim> RS(R);
  return RS -= S;
}

template<typename T, int dim>
constexpr inline TensorT<T, dim> operator-(const TensorT<T, dim> &R, const SymTensorT<T, dim> &S) {
  TensorT<T, dim> RS(R);
  return RS -= S;
}

template<typename T, int dim>
constexpr inline TensorT<T, dim> operator-(const TensorT<T, dim> &R, const TensorRowT<T, dim> &S) {
  TensorT<T, dim> RS(R);
  return RS -= S;
}

template<typename T, int dim>
constexpr inline TensorT<T, dim> operator-(const TensorRowT<T, dim> &R, const TensorT<T, dim> &S) {
  TensorT<T, dim> RS(R);
  return RS -= S;
}

template<typename T, int dim>
constexpr inline TensorT<T, dim> operator-(const TensorT<T, dim> &R, const TensorComponentT<T> &S) {
  TensorT<T, dim> RS(R);
  return RS -= S;
}

template<typename T, int dim>
constexpr inline TensorT<T, dim> operator-(const TensorComponentT<T> &R, const TensorT<T, dim> &S) {
  TensorT<T, dim> RS(-S);
  return RS += R;
}

template<typename T, int dim>
constexpr inline TensorT<T, dim> operator-(const TensorT<T, dim> &S) {
  TensorT<T, dim> R(S);
  return R *= T(-1.0);
}

template<typename T, typename TT, int dim>
constexpr inline TensorT<T, dim> operator*(const TT &b, const TensorT<T, dim> &S) {
  TensorT<T, dim> bS(S);
  return bS *= b;
}

template<typename T, typename TT, int dim>
constexpr inline TensorT<T, dim> operator*(const TensorT<T, dim> &S, const TT &b) {
  TensorT<T, dim> bS(S);
  return bS *= b;
}

template<typename T, typename TT, int dim>
constexpr inline TensorT<T, dim> operator/(const TensorT<T, dim> &S, const TT &b) {
  TensorT<T, dim> bS(S);
  return bS /= b;
}

template<typename T, int dim>
constexpr inline TensorT<T, dim> transpose(const TensorT<T, dim> &S) {
  TensorT<T, dim> R(S);
  return R.transpose();
}

template<typename T, int dim>
constexpr inline T trace(const TensorT<T, dim> &S) {
  return S.trace();
}

template<typename T, int dim>
constexpr inline T det(const TensorT<T, dim> &S) {
  return S.det();
}

template<typename T, int dim>
constexpr inline TensorT<T, dim> Invert(const TensorT<T, dim> &S) {
  TensorT<T, dim> R(S);
  return R.Invert();
}

template<typename T, int dim>
constexpr inline TensorT<T, dim> transposeInvert(const TensorT<T, dim> &S) {
  TensorT<T, dim> R(S);
  return R.transposeInvert();
}

template<typename T, int dim>
constexpr inline TensorT<T, dim> sym(const TensorT<T, dim> &S) {
  TensorT<T, dim> R(S);
  return R.sym();
}

template<typename T, int dim>
constexpr inline TensorT<T, dim> skew(const TensorT<T, dim> &S) {
  TensorT<T, dim> R(S);
  return R.skew();
}

template<typename T, int sDim, int tDim>
constexpr inline VectorFieldT<T, sDim, tDim> operator*(const TensorT<T, sDim> &S,
                                                       const VectorFieldT<T, sDim, tDim> &v) {
  return S.multiply(v);
}

template<typename T, int sDim>
constexpr inline VectorFieldT<T, sDim> operator*(const TensorT<T, sDim> &S,
                                                 const VectorFieldComponentT<T> &v) {
  return S.multiply(v);
}

template<typename T, int sDim, int tDim>
constexpr inline PointT<T, sDim, tDim> operator*(const TensorT<T, sDim> &S,
                                                 const PointT<T, sDim, tDim> &D) {
  return S.multiply(D);
}

template<typename T, int dim>
constexpr inline VectorFieldT<T, dim> operator*(const VectorFieldT<T, dim> &v,
                                                const TensorT<T, dim> &S) {
  return S.applyTransposed(v);
}

template<typename T, int dim>
constexpr inline VectorFieldT<T, dim> operator*(const VectorFieldComponentT<T> &v,
                                                const TensorT<T, dim> &S) {
  return S.applyTransposed(v);
}

template<typename T, int dim>
constexpr inline TensorT<T, dim> operator*(const TensorT<T, dim> &R, const TensorT<T, dim> &S) {
  return TensorT<T, dim>(R, S);
}

template<typename T, int dim>
T Frobenius(const TensorT<T, dim> &R, const TensorT<T, dim> &S) {
  T sum{};
  for (int n = 0; n < dim; ++n)
    sum += R[n] * S[n];
  return sum;
}

template<typename T, int dim>
constexpr inline T norm(const TensorT<T, dim> &S) {
  return S.norm();
}

template<int dim>
constexpr inline double maxnorm(const TensorT<double, dim> &S) {
  return S.maxnorm();
}

/// templates cannot be used for the stream since Saver has the same stucture
template<typename T, int dim>
constexpr std::ostream &operator<<(std::ostream &os, const TensorT<T, dim> &S) {
  return S.print(os);
}

template<typename T, int dim>
constexpr inline Saver &operator<<(Saver &saver, const TensorT<T, dim> &S) {
  return S.save(saver);
}

template<typename T, int dim>
constexpr inline Loader &operator>>(Loader &loader, TensorT<T, dim> &S) {
  return S.load(loader);
}

template<int dim>
constexpr inline double maxelem(const TensorT<double, dim> &S) {
  double m = 0.0;
  for (int i = 0; i < dim; ++i)
    for (int j = 0; j < dim; ++j)
      if (abs(S[i][j]) > m) m = abs(S[i][j]);
  return m;
}

constexpr inline TensorT<double, 3> MatrixGauss(const TensorT<double, 3> &D,
                                                const TensorT<double, 3> &N) {
  // solves DF = N bei Gaussian elimination for F
  double c0 = D[1][0] / D[0][0];
  double c1 = D[2][0] / D[0][0];
  double c2 = (D[2][1] - c1 * D[0][1]) / (D[1][1] - c0 * D[0][1]);

  TensorT<double, 3> F;
  for (int j = 0; j < 3; ++j) {
    F[2][j] = (N[2][j] - c1 * N[0][j] - c2 * (N[1][j] - c0 * N[0][j]))
              / (D[2][2] - c1 * D[0][2] - c2 * (D[1][2] - c0 * D[0][2]));
    F[1][j] =
        (N[1][j] - c0 * N[0][j] - (D[1][2] - c0 * D[0][2]) * F[2][j]) / (D[1][1] - c0 * D[0][1]);
    F[0][j] = (N[0][j] - D[0][2] * F[2][j] - D[0][1] * F[1][j]) / D[0][0];
  }
  return F;
}

inline TensorT<double, 3> exp(const TensorT<double, 3> &T) {
  // exp(T) by scaling and squaring algorithm using Pade approximation
  // see Golub: 'Matrix Computations' pg.573
  double nt = maxnorm(T);
  double j = max(0.0, 1.0 + floor(log(nt) / log(2.0))); // log_2(|T|_infty)
  // scaling
  TensorT<double, 3> A = (1.0 / pow(2.0, j)) * T;
  double delta = Eps + Eps * norm(T);

  int q = 1; // corresponds to delta >= 1/6
  if (delta >= 1.0 / 1440) {
    q = 2;
    //        dout(10) << "q = 2: eps = " << 1.0 / 1440 << "\n";
  } else if (delta >= 1.0 / 806400) {
    q = 3;
    //        dout(10) << "q = 3: eps = " << 1.0 / 806400 << "\n";
  } else if (delta >= 1.0 / 812851200) {
    q = 4;
    //        dout(10) << "q = 4: eps = " << 1.0 / 812851200 << "\n";
  } else if (delta >= 7.76665066513e-13) {
    q = 5;
    //        dout(10) << "q = 5: eps = " << 7.76665066513e-13 << "\n";
  } else if (delta >= 3.39451515084e-16) {
    q = 6;
    //        dout(10) << "q = 6: eps = " << 3.39451515084e-16 << "\n";
  } else if (delta >= 1.08798562527e-19) {
    q = 7;
    //        dout(10) << "q = 6: eps = " << 1.08798562527e-19 << "\n";
  } else THROW("exp: Selected accuracy unrealistic; file: Tensor.h\n")

  TensorT<double, 3> D(1, 0, 0, 0, 1, 0, 0, 0, 1);
  TensorT<double, 3> N(1, 0, 0, 0, 1, 0, 0, 0, 1);
  TensorT<double, 3> X(1, 0, 0, 0, 1, 0, 0, 0, 1);
  double c = 1;
  int altern_sign = 1;
  for (int k = 1; k < q; ++k) {
    c = (q - k + 1) * c / ((2 * q - k + 1) * k);
    X = A * X;
    altern_sign *= -1;
    D += altern_sign * c * X;
    N += c * X;
  }

  // solve DF = N bei Gaussian Elimination for F
  TensorT<double, 3> F = MatrixGauss(D, N);
  // squaring
  for (int k = 1; k <= int(j); ++k)
    F = F * F;

  return F;
}

constexpr inline int JacobiRotation(TensorT<double, 3> &S, TensorT<double, 3> &U, int i, int j) {
  const double EPS = 1e-10;
  double a_ij = S[i][j];
  if (abs(a_ij) < EPS) return 1;
  double a_ii = S[i][i];
  double a_jj = S[j][j];
  TensorT<double, 3> J(1, 0, 0, 0, 1, 0, 0, 0, 1);
  double t = (a_jj - a_ii) / (2.0 * a_ij);
  if (abs(t) < EPS) t = 1;
  else if (t > 0) t = 1 / (t + sqrt(t * t + 1));
  else t = 1 / (t - sqrt(t * t + 1));
  J[i][i] = J[j][j] = 1 / sqrt(1 + t * t);
  J[j][i] = t * J[j][j];
  J[i][j] = -J[j][i];
  S = J * S * transpose(J);
  U = J * U;
  return 0;
}

inline TensorT<double, 3> diagonalize(TensorT<double, 3> &S) {
  TensorT<double, 3> S_old(S);
  const int ITER = 100;
  TensorT<double, 3> U(1, 0, 0, 0, 1, 0, 0, 0, 1);
  for (int i = 0; i <= ITER; ++i) {
    double a_10 = abs(S[1][0]);
    double a_20 = abs(S[2][0]);
    double a_21 = abs(S[2][1]);
    if (a_10 > a_20) {
      if (a_10 > a_21) {
        if (JacobiRotation(S, U, 1, 0)) return U;
      } else if (JacobiRotation(S, U, 2, 1)) return U;
    } else {
      if (a_20 > a_21) {
        if (JacobiRotation(S, U, 2, 0)) return U;
      } else if (JacobiRotation(S, U, 2, 1)) return U;
    }
  }
  THROW("No convergence in diagonalize; file Tensor.hpp")
  //    mout << S_old << endl;
}

inline TensorT<double, 3> log(TensorT<double, 3> E) {
  TensorT<double, 3> U = diagonalize(E);
  E[0][0] = log(E[0][0]);
  E[1][1] = log(E[1][1]);
  E[2][2] = log(E[2][2]);
  return transpose(U) * E * U;
}

inline TensorT<double, 3> pow(const TensorT<double, 3> &T, double x) {
  // determines T^x=U^T*D^x*U for symmetric 3x3 tensors
  TensorT<double, 3> A = T;
  TensorT<double, 3> U = diagonalize(A);
  A[0][0] = pow(A[0][0], x);
  A[1][1] = pow(A[1][1], x);
  A[2][2] = pow(A[2][2], x);
  return transpose(U) * A * U;
}

inline TensorT<double, 3> sqrt(TensorT<double, 3> T) {
  TensorT<double, 3> A = T;
  TensorT<double, 3> U = diagonalize(A);
  A[0][0] = sqrt(A[0][0]);
  A[1][1] = sqrt(A[1][1]);
  A[2][2] = sqrt(A[2][2]);
  return transpose(U) * A * U;
}

template<typename T>
constexpr inline TensorT<T, 3> anti(const VectorFieldT<T, 3> &V) {
  return TensorT<T, 3>(T{}, -V[2], V[1], V[2], T{}, -V[0], -V[1], V[0], T{});
}

template<typename T>
constexpr inline TensorT<T, 3> anti(T a, T b, T c) {
  return anti(VectorFieldT<T, 3>(a, b, c));
}

template<typename T>
constexpr inline VectorFieldT<T, 3> axl(const TensorT<T, 3> &S) {
  return VectorFieldT<T, 3>(S[2][1], S[0][2], S[1][0]);
}

template<typename T, int sDim>
constexpr inline TensorT<T, sDim> Product(const VectorFieldT<T, sDim> &V,
                                          const VectorFieldT<T, sDim> &W) {
  if constexpr (sDim == 1) return TensorT<T, sDim>(V[0] * W[0]);
  else if constexpr (sDim == 2)
    return TensorT<T, sDim>(V[0] * W[0], V[0] * W[1], V[1] * W[0], V[1] * W[1]);
  else if constexpr (sDim == 3)
    return TensorT<T, sDim>(V[0] * W[0], V[0] * W[1], V[0] * W[2], V[1] * W[0], V[1] * W[1],
                            V[1] * W[2], V[2] * W[0], V[2] * W[1], V[2] * W[2]);
  else THROW("Dyadic Product is not defined for dimensions higher than 3")
}

template<typename T, int sDim>
constexpr TensorT<T, sDim> Cross(const TensorT<T, sDim> R, const TensorT<T, sDim> S) {
  if constexpr (sDim != 3) THROW("Cross Product not defined for dimensions other than 3")

  return TensorT<
      T, sDim>(R[1][1] * S[2][2] - R[1][2] * S[2][1] + R[2][2] * S[1][1] - R[2][1] * S[1][2],
               R[1][2] * S[2][0] - R[1][0] * S[2][2] + R[2][0] * S[1][2] - R[2][2] * S[1][0],
               R[1][0] * S[2][1] - R[1][1] * S[2][0] + R[2][1] * S[1][0] - R[2][0] * S[1][1],
               R[0][2] * S[2][1] - R[0][1] * S[2][2] + R[2][1] * S[0][2] - R[2][2] * S[0][1],
               R[2][2] * S[0][0] - R[2][0] * S[0][2] + R[0][0] * S[2][2] - R[0][2] * S[2][0],
               R[2][0] * S[0][1] - R[2][1] * S[0][0] + R[0][1] * S[2][0] - R[0][0] * S[2][1],
               R[0][1] * S[1][2] - R[0][2] * S[1][1] + R[1][2] * S[0][1] - R[1][1] * S[0][2],
               R[0][2] * S[1][0] - R[0][0] * S[1][2] + R[1][0] * S[0][2] - R[1][2] * S[0][0],
               R[0][0] * S[1][1] - R[0][1] * S[1][0] + R[1][1] * S[0][0] - R[1][0] * S[0][1]);
}

template<typename T, int sDim>
constexpr TensorT<T, sDim> Cross(const TensorT<T, sDim> R, const TensorRowT<T, sDim> S) {
  return Cross(R, TensorT(S));
}

template<typename T, int dim>
constexpr inline TensorT<T, dim> diag(const VectorFieldT<T, dim> &D) {
  TensorT<T, dim> S;
  for (int n = 0; n < dim; ++n)
    S[n][n] = D[n];
  return S;
}

template<typename T, int dim>
constexpr inline VectorFieldT<T, dim> diag(const TensorT<T, dim> &S) {
  VectorFieldT<T, dim> D;
  for (int n = 0; n < dim; ++n)
    D[n] = S[n][n];
  return D;
}

/***
 * TensorRow operators
 */
template<typename T, int dim>
constexpr inline TensorRowT<T, dim> operator*(const TensorRowT<T, dim> &R,
                                              const TensorT<T, dim> &S) {
  VectorFieldT<T, dim> t{};
  for (int n = 0; n < dim; ++n) {
    t[n] = R.Value() * S.Col(n);
  }
  return TensorRowT(R.Component(), t);
}

template<typename T, int dim>
constexpr inline TensorT<T, dim> operator*(const TensorT<T, dim> &S, const TensorRowT<T, dim> &R) {
  TensorT<T, dim> C{};
  auto col = S.Col(R.Component());
  for (int n = 0; n < dim; ++n) {
    C[n] = R.Value()[n] * col;
  }
  return C.transpose();
}

template<typename T, int dim>
T Frobenius(const TensorT<T, dim> &S, const TensorRowT<T, dim> &R) {
  return S[R.Component()] * R.Value();
}

template<typename T, int dim>
T Frobenius(const TensorRowT<T, dim> &R, const TensorT<T, dim> &S) {
  return Frobenius(S, R);
}

template<typename T, int dim>
constexpr inline TensorT<T, dim> sym(const TensorRowT<T, dim> &R) {
  VectorFieldT<T, dim> r = 0.5 * R.Value();
  TensorT<T, dim> S(R.Component(), r);
  for (int n = 0; n < dim; ++n) {
    S[n][R.Component()] += r[n];
  }
  return S;
}

// Returns R^T*S
template<typename T, int dim>
constexpr inline TensorT<T, dim> ColRowProduct(const TensorRowT<T, dim> &R,
                                               const TensorRowT<T, dim> &S) {
  return (R.Component() == S.Component()) * Product(R.Value(), S.Value());
}

/***
 * TensorComponent operators
 */
template<typename T, int dim>
constexpr inline TensorRowT<T, dim> operator*(const TensorComponentT<T> &C,
                                              const TensorT<T, dim> &S) {
  return TensorRowT(C.Row(), C.Value() * S.Row(C.Column()));
}

template<typename T, int dim>
constexpr inline TensorT<T, dim> operator*(const TensorT<T, dim> &S, const TensorComponentT<T> &C) {
  TensorT<T, dim> R{};
  auto col = S.Col(C.Row());
  R[C.Column()] = C.Value() * col;
  return R.transpose();
}

template<typename T, int dim>
T Frobenius(const TensorT<T, dim> &S, const TensorComponentT<T> &R) {
  return S[R.Row()][R.Column()] * R.Value();
}

template<typename T, int dim>
T Frobenius(const TensorComponentT<T> &R, const TensorT<T, dim> &S) {
  return Frobenius(S, R);
}

template<typename T, int dim>
constexpr inline TensorT<T, dim> sym(const TensorComponentT<T> &R) {
  double val = 0.5 * R.Value();
  TensorT<T, dim> S(R.Row(), R.Column(), val);
  S[R.Column()][R.Value()] += val;
  return S;
}

/************/

constexpr TensorT<> One = diag(VectorFieldT<>(1.0));

constexpr TensorT<> Zero;

constexpr inline TensorT<double, 3> dev(const TensorT<double, 3> &T) {
  return T - trace(T) / 3.0 * diag(VectorFieldT<double, 3>(1.0));
}

template<typename T = double, int dim = SpaceDimension>
using VelocityGradientT = TensorT<T, dim>;

using Tensor = TensorT<>;

using VelocityGradient = VelocityGradientT<>;

#ifdef BUILD_IA

using IATensor = TensorT<IAInterval, SpaceDimension>;

using IAVelocityGradient = VelocityGradientT<IAInterval, SpaceDimension>;

template<int sDim>
TensorT<double, sDim> mid(const TensorT<IAInterval, sDim> &IA_T) {
  TensorT<double, sDim> T;
  for (int n = 0; n < sDim; ++n)
    T[n] = mid(IA_T[n]);
  return T;
}

template<int sDim>
TensorT<double, sDim> sup(const TensorT<IAInterval, sDim> &IA_T) {
  TensorT<double, sDim> T;
  for (int n = 0; n < sDim; ++n)
    T[n] = sup(IA_T[n]);
  return T;
}

template<int sDim>
TensorT<double, sDim> inf(const TensorT<IAInterval, sDim> &IA_T) {
  TensorT<double, sDim> T;
  for (int n = 0; n < sDim; ++n)
    T[n] = inf(IA_T[n]);
  return T;
}

template<int sDim>
inline TensorT<IAInterval, sDim> operator+(const TensorT<double, sDim> &x,
                                           const TensorT<IAInterval, sDim> &y) {
  TensorT<IAInterval, sDim> IAx(x);
  return IAx += y;
}

template<int sDim>
inline TensorT<IAInterval, sDim> operator+(const TensorT<IAInterval, sDim> &x,
                                           const TensorT<double, sDim> &y) {
  TensorT<IAInterval, sDim> IAx(x);
  return IAx += y;
}

template<int sDim>
inline TensorT<IAInterval, sDim> operator+(const SymTensorT<double, sDim> &x,
                                           const TensorT<IAInterval, sDim> &y) {
  TensorT<IAInterval, sDim> IAx(x);
  return IAx += y;
}

template<int sDim>
inline TensorT<IAInterval, sDim> operator+(const TensorT<IAInterval, sDim> &x,
                                           const SymTensorT<double, sDim> &y) {
  TensorT<IAInterval, sDim> IAx(x);
  return IAx += y;
}

template<int sDim>
inline TensorT<IAInterval, sDim> operator+(const TensorT<double, sDim> &x,
                                           const SymTensorT<IAInterval, sDim> &y) {
  TensorT<IAInterval, sDim> IAx(x);
  return IAx += y;
}

template<int sDim>
inline TensorT<IAInterval, sDim> operator+(const SymTensorT<IAInterval, sDim> &x,
                                           const TensorT<double, sDim> &y) {
  TensorT<IAInterval, sDim> IAx(x);
  return IAx += y;
}

template<int sDim>
inline TensorT<IAInterval, sDim> operator-(const TensorT<double, sDim> &x,
                                           const TensorT<IAInterval, sDim> &y) {
  TensorT<IAInterval, sDim> IAx(x);
  return IAx -= y;
}

template<int sDim>
inline TensorT<IAInterval, sDim> operator-(const TensorT<IAInterval, sDim> &x,
                                           const TensorT<double, sDim> &y) {
  TensorT<IAInterval, sDim> IAx(x);
  return IAx -= y;
}

template<int sDim>
inline TensorT<IAInterval, sDim> operator-(const SymTensorT<double, sDim> &x,
                                           const TensorT<IAInterval, sDim> &y) {
  TensorT<IAInterval, sDim> IAx(x);
  return IAx -= y;
}

template<int sDim>
inline TensorT<IAInterval, sDim> operator-(const TensorT<IAInterval, sDim> &x,
                                           const SymTensorT<double, sDim> &y) {
  TensorT<IAInterval, sDim> IAx(x);
  return IAx -= y;
}

template<int sDim>
inline TensorT<IAInterval, sDim> operator-(const TensorT<double, sDim> &x,
                                           const SymTensorT<IAInterval, sDim> &y) {
  TensorT<IAInterval, sDim> IAx(x);
  return IAx -= y;
}

template<int sDim>
inline TensorT<IAInterval, sDim> operator-(const SymTensorT<IAInterval, sDim> &x,
                                           const TensorT<double, sDim> &y) {
  TensorT<IAInterval, sDim> IAx(x);
  return IAx -= y;
}

template<int sDim>
inline TensorT<IAInterval, sDim> operator*(const TensorT<double, sDim> &x, const IAInterval &b) {
  TensorT<IAInterval, sDim> IAx(x);
  return IAx *= b;
}

template<int sDim>
inline TensorT<IAInterval, sDim> operator*(const IAInterval &b, const TensorT<double, sDim> &x) {
  TensorT<IAInterval, sDim> IAx(x);
  return IAx *= b;
}

template<int sDim>
inline TensorT<IAInterval, sDim> operator/(const TensorT<double, sDim> &x, const IAInterval &b) {
  TensorT<IAInterval, sDim> IAx(x);
  return IAx /= b;
}

template<int sDim>
VectorFieldT<IAInterval, sDim> operator*(const TensorT<double, sDim> &S,
                                         const VectorFieldT<IAInterval, sDim> &v) {
  if constexpr (sDim == 1) return VectorFieldT<IAInterval, sDim>(S[0] * v);
  if constexpr (sDim == 2) return VectorFieldT<IAInterval, sDim>(S[0] * v, S[1] * v);
  if constexpr (sDim == 3) return VectorFieldT<IAInterval, sDim>(S[0] * v, S[1] * v, S[2] * v);
}

template<int sDim>
inline VectorFieldT<IAInterval, sDim> operator*(const TensorT<IAInterval, sDim> &S,
                                                const VectorFieldT<double, sDim> &v) {
  if constexpr (sDim == 1) return VectorFieldT<IAInterval, sDim>(S[0] * v);
  if constexpr (sDim == 2) return VectorFieldT<IAInterval, sDim>(S[0] * v, S[1] * v);
  if constexpr (sDim == 3) return VectorFieldT<IAInterval, sDim>(S[0] * v, S[1] * v, S[2] * v);
}

template<int dim>
constexpr inline TensorT<IAInterval, dim> operator*(const TensorT<IAInterval, dim> &R,
                                                    const TensorT<double, dim> &S) {
  return TensorT<IAInterval, dim>(R, S);
}

template<int dim>
constexpr inline TensorT<IAInterval, dim> operator*(const TensorT<double, dim> &R,
                                                    const TensorT<IAInterval, dim> &S) {
  return TensorT<IAInterval, dim>(R, S);
}

template<int sDim>
constexpr inline TensorT<IAInterval, sDim> Product(const VectorFieldT<IAInterval, sDim> &V,
                                                   const VectorFieldT<double, sDim> &W) {
  if constexpr (sDim == 1) return TensorT<IAInterval, sDim>(V[0] * W[0]);
  else if constexpr (sDim == 2)
    return TensorT<IAInterval, sDim>(V[0] * W[0], V[0] * W[1], V[1] * W[0], V[1] * W[1]);
  else if constexpr (sDim == 3)
    return TensorT<IAInterval, sDim>(V[0] * W[0], V[0] * W[1], V[0] * W[2], V[1] * W[0],
                                     V[1] * W[1], V[1] * W[2], V[2] * W[0], V[2] * W[1],
                                     V[2] * W[2]);
  else THROW("Dyadic Product is not defined for dimensions higher than 3")
}

template<int sDim>
constexpr inline TensorT<IAInterval, sDim> Product(const VectorFieldT<double, sDim> &V,
                                                   const VectorFieldT<IAInterval, sDim> &W) {
  if constexpr (sDim == 1) return TensorT<IAInterval, sDim>(V[0] * W[0]);
  else if constexpr (sDim == 2)
    return TensorT<IAInterval, sDim>(V[0] * W[0], V[0] * W[1], V[1] * W[0], V[1] * W[1]);
  else if constexpr (sDim == 3)
    return TensorT<IAInterval, sDim>(V[0] * W[0], V[0] * W[1], V[0] * W[2], V[1] * W[0],
                                     V[1] * W[1], V[1] * W[2], V[2] * W[0], V[2] * W[1],
                                     V[2] * W[2]);
  else THROW("Dyadic Product is not defined for dimensions higher than 3")
}

#endif // BUILD_IA

#endif // of #ifndef _TENSOR_H_
