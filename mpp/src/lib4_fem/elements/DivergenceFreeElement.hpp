#ifndef DIVERGENCEFREEELEMENT_HPP
#define DIVERGENCEFREEELEMENT_HPP

#include "ArgyrisElement.hpp"

template<typename TT = double, int sDim = SpaceDimension, int tDim = TimeDimension>
class DivergenceFreeElementT : public ArgyrisBasicElementT<TT, sDim, tDim> {
  ShapeValues<VectorFieldT<TT, sDim>> gradient;
  ShapeValues<SymTensorT<TT, sDim>> hessian;
public:
  DivergenceFreeElementT(const VectorMatrixBase &base, const Cell &c);

  DivergenceFreeElementT(const VectorMatrixBase &base,
                         const BasicElementT<TT, sDim, tDim> &baseElement);

  VelocityT<TT, sDim> VelocityField(int q, int i, int k) const;

  VelocityT<TT, sDim> VelocityField(int q, const ShapeId &iter) const;

  VelocityT<TT, sDim> VelocityField(const PointT<TT, sDim, tDim> &z, int i, int k) const;

  VelocityT<TT, sDim> VelocityField(const PointT<TT, sDim, tDim> &z, const ShapeId &iter) const;

  VelocityT<TT, sDim> VelocityField(int q, const Vector &u, int m = 0) const;

  VelocityT<TT, sDim> VelocityField(const PointT<TT, sDim, tDim> &z, const Vector &u,
                                    int m = 0) const;

  VelocityGradientT<TT, sDim> VelocityFieldGradient(int q, int i, int k) const;

  VelocityGradientT<TT, sDim> VelocityFieldGradient(int q, const ShapeId &iter) const;

  VelocityGradientT<TT, sDim> VelocityFieldGradient(const PointT<TT, sDim, tDim> &z, int i,
                                                    int k) const;

  VelocityGradientT<TT, sDim> VelocityFieldGradient(const PointT<TT, sDim, tDim> &z,
                                                    const ShapeId &iter) const;

  VelocityGradientT<TT, sDim> VelocityFieldGradient(int q, const Vector &u, int m = 0) const;

  VelocityGradientT<TT, sDim> VelocityFieldGradient(const PointT<TT, sDim, tDim> &z,
                                                    const Vector &u, int m = 0) const;

  friend class TestDivergenceFreeElement;
};

using DivergenceFreeElement = DivergenceFreeElementT<>;

#ifdef BUILD_IA

using IADivergenceFreeElement = DivergenceFreeElementT<IAInterval, SpaceDimension, TimeDimension>;

#endif

template<int sDim = SpaceDimension, int tDim = TimeDimension>
using DivergenceFreeVectorAccessT = VectorAccessGeneral<DivergenceFreeElementT<double, sDim, tDim>>;

using DivergenceFreeVectorAccess = DivergenceFreeVectorAccessT<>;

template<int sDim = SpaceDimension, int tDim = TimeDimension>
using DivergenceFreeMatrixAccessT = MatrixAccessGeneral<DivergenceFreeElementT<double, sDim, tDim>>;

using DivergenceFreeMatrixAccess = DivergenceFreeMatrixAccessT<>;

#endif // DIVERGENCEFREEELEMENT_HPP
