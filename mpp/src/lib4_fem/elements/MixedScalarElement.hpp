#ifndef MIXEDSCALARELEMENT_HPP
#define MIXEDSCALARELEMENT_HPP

#include "Element.hpp"

template<typename TT = double, int sDim = SpaceDimension, int tDim = TimeDimension>
class MixedScalarElementT : public MixedElementT<TT, sDim, tDim> {
protected:
  vector<const ShapeValues<TT> *> shapeValues{};
  vector<ShapeValues<VectorFieldT<TT, sDim>>> shapeGradients{};

  MixedScalarElementT(const VectorMatrixBase &base, const Cell &c, int face,
                      const vector<ShapeValues<TT>> &faceValues);

  void init();
public:
  MixedScalarElementT(const VectorMatrixBase &g, const Cell &c);

  MixedScalarElementT(const VectorMatrixBase &g, const BasicElementT<TT, sDim, tDim> &baseElement);

  TT Value(int q, int n, int i) const;

  TT Value(int q, int n, const ShapeId &iter) const;

  TT Value(const PointT<TT, sDim, tDim> &z, int n, int i) const;

  TT Value(const PointT<TT, sDim, tDim> &z, int n, const ShapeId &iter) const;

  TT Value(int q, const Vector &u, int n, int k = 0) const;

  TT Value(const PointT<TT, sDim, tDim> &z, const Vector &u, int n, int k = 0) const;

  const VectorFieldT<TT, sDim> &Derivative(int q, int n, int i) const;

  const VectorFieldT<TT, sDim> &Derivative(int q, int n, const ShapeId &iter) const;

  VectorFieldT<TT, sDim> Derivative(const PointT<TT, sDim, tDim> &z, int n, int i) const;

  VectorFieldT<TT, sDim> Derivative(const PointT<TT, sDim, tDim> &z, int n,
                                    const ShapeId &iter) const;

  VectorFieldT<TT, sDim> Derivative(int q, const Vector &u, int n, int k = 0) const;

  VectorFieldT<TT, sDim> Derivative(const PointT<TT, sDim, tDim> &z, const Vector &u, int n,
                                    int k = 0) const;
};

using MixedScalarElement = MixedScalarElementT<>;

template<int sDim = SpaceDimension, int tDim = TimeDimension>
using MixedScalarVectorAccessT =
    MixedVectorAccessEquallyDistributed<MixedScalarElementT<double, sDim, tDim>>;

using MixedScalarVectorAccess = MixedScalarVectorAccessT<>;

template<int sDim = SpaceDimension, int tDim = TimeDimension>
using MixedScalarMatrixAccessT =
    MixedMatrixAccessEquallyDistributed<MixedScalarElementT<double, sDim, tDim>>;

using MixedScalarMatrixAccess = MixedScalarMatrixAccessT<>;

#ifdef BUILD_IA

using IAMixedScalarElement = MixedScalarElementT<IAInterval, SpaceDimension, TimeDimension>;

#endif


template<typename TT = double, int sDim = SpaceDimension, int tDim = TimeDimension>
class MixedScalarFaceElementT :
    public IMixedFaceElementT<TT, sDim, tDim>,
    public MixedScalarElementT<TT, sDim, tDim> {
public:
  MixedScalarFaceElementT(const VectorMatrixBase &base, const Cell &c, int face);
};

using MixedScalarFaceElement = MixedScalarFaceElementT<>;

#ifdef BUILD_IA

using IAMixedScalarFaceElement = MixedScalarFaceElementT<IAInterval, SpaceDimension, TimeDimension>;

#endif

#endif // MIXEDSCALARELEMENT_HPP
