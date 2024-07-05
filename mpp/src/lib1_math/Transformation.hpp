#ifndef TRANSFORMATION_H
#define TRANSFORMATION_H

#include "Tensor.hpp"

/**
 * Note: in the following J is the Jacobian of the transformation
 *       T: \hat K \rightarrow K from the reference element \hat K
 *       to K
 * in J[][] is the Jacobian
 * in JTinv[][] is the transposed inverse of the jacobian, i.e. JTinv=J^{-T}
 *
 * JTinv*z             is  J^{-T} * z
 * z*JTinv             is  J^{-1} * z
 * operator() (P)      is  J^T*P
 * ApplyJ              is  J*P
 */

/*template<typename T = double, int dim = SpaceDimension, int tDim = TimeDimension>
class ITransformationT {
  virtual T Det() const = 0;
  virtual VectorFieldT<T, dim, tDim> ApplyJinvTransposed(const VectorFieldT<T, dim, tDim> &z) const
= 0; virtual VectorFieldT<T, dim, tDim> ApplyJinv(const VectorFieldT<T, dim, tDim> &z) const  = 0;
  virtual VectorFieldT<T, dim, tDim> ApplyJ(const VectorFieldT<T, dim, tDim> &z) const = 0;
  virtual PointT<T, dim, tDim> ApplyJ(const PointT<T, dim, tDim> &z) const = 0;
};*/

template<typename TT = double, int dim = SpaceDimension>
class TransformationT {
protected:
  TensorT<TT> space_J;
#if TimeDimension > 0
  TT time_J{};
#endif
  TensorT<TT> JTinv;
  TT space_detJ;
public:
  TransformationT() = default;

  explicit TransformationT(const std::vector<PointT<TT>> &z) {
    for (int n = 0; n < std::min(dim, int(z.size())); ++n) {
      for (int m = 0; m < dim; ++m) {
        space_J[n][m] = z[n][m];
      }
    }
    for (int n = std::min(dim, int(z.size())); n < dim; ++n) {
      for (int m = 0; m < dim; ++m) {
        space_J[n][m] = (n == m) ? 1 : 0;
      }
    }
    space_detJ = det(space_J);
#if TimeDimension > 0
    time_J = 1.0;
#endif
    JTinv = transposeInvert(space_J);
  }

  TransformationT(TransformationT<TT> SpaceT, TT timeJ) :
      space_J(SpaceT.GetJ()), JTinv(SpaceT.GetJTinv()), space_detJ(SpaceT.Det()) {
    if constexpr (TimeDimension == 0) { THROW("Do not call SpaceTransformation with TDim = 0") }
#if TimeDimension > 0
    time_J = timeJ;
#endif
  }

  TT Det() const {
#if TimeDimension > 0
    return space_detJ * time_J;
#else
    return space_detJ;
#endif
  }

  TT DetTime() const {
#if TimeDimension > 0
    return time_J;
#else
    THROW("Do not call DetTime with TimeDimension == 0")
#endif
  }

  TT DetSpace() const { return space_detJ; }

  VectorFieldT<TT> operator[](int i) const { return JTinv[i]; }

  VectorFieldT<TT> operator()(int i) const { return space_J[i]; }

  VectorFieldT<TT> ApplyJinvTransposed(const VectorFieldT<TT> &z) const {
    VectorFieldT<TT> res = JTinv * z;
#if TimeDimension > 0
    res.t() = z.t() / time_J;
#endif
    return res;
  }

  PointT<TT> ApplyJinvTransposed(const PointT<TT> &z) const {
    PointT<TT> p = JTinv * z;
#if TimeDimension > 0
    p.WithT(z.t() / time_J);
#endif
    return p;
  }

  VectorFieldT<TT> operator()(const VectorFieldT<TT> &z) const {
    VectorFieldT<TT> res = space_J.applyTransposed(z);
#if TimeDimension > 0
    res.t() = z.t() * time_J;
#endif
    return res;
  }

  VectorFieldT<TT> ApplyJinv(const VectorFieldT<TT, dim> &z) const {
    VectorFieldT<TT> res = JTinv.applyTransposed(z);
#if TimeDimension > 0
    res.t() = z.t() / time_J;
#endif
    return res;
  }

  PointT<TT> ApplyJinv(const PointT<TT> &z) const {
    PointT<TT> p = JTinv.applyTransposed(z);
#if TimeDimension > 0
    p.WithT(z.t() / time_J);
#endif
    return p;
  }

  PointT<TT> operator()(const PointT<TT> &z) const {
    PointT<TT> p = space_J.applyTransposed(z);
#if TimeDimension > 0
    p.WithT(z.t() * time_J);
#endif
    return p;
  }

  VectorFieldT<TT> ApplyJ(const VectorFieldT<TT> &z) const {
    VectorFieldT<TT> res = space_J * z;
#if TimeDimension > 0
    res.t() = z.t() * time_J;
#endif
    return res;
  }

  PointT<TT> ApplyJ(const PointT<TT> &z) const {
    PointT<TT> p = space_J * z;
#if TimeDimension > 0
    p.WithT(z.t() * time_J);
#endif
    return p;
  }

  std::ostream &print(std::ostream &s) const {
    return s << JTinv << "--------------------" << std::endl
             << space_J << "DET= " << space_detJ << std::endl
             << std::endl;
  }

  TensorT<TT> GetJ() { return space_J; }

  TensorT<TT> GetJTinv() { return JTinv; }
};

template<typename TT>
inline VectorFieldT<TT> operator*(const TransformationT<TT> &T, const VectorFieldT<TT> &z) {
  return T.ApplyJinvTransposed(z);
}

template<typename TT>
inline PointT<TT> operator*(const TransformationT<TT> &T, const PointT<TT> &z) {
  return T.ApplyJinvTransposed(z);
}

template<typename TT>
inline VectorFieldT<TT> operator*(const VectorFieldT<TT> &z, const TransformationT<TT> &T) {
  return T.ApplyJinv(z);
}

template<typename TT>
PointT<TT> operator*(const PointT<TT> &z, const TransformationT<TT> &T) {
  return T.ApplyJinv(z);
}

template<typename TT, int dim>
std::ostream &operator<<(std::ostream &s, const TransformationT<TT, dim> &T) {
  return T.print(s);
}

using Transformation = TransformationT<>;

#ifdef BUILD_IA

using IATransformation = TransformationT<IAInterval, SpaceDimension>;

#endif

#endif // TRANSFORMATION_H
