#ifndef MIXEDVECTORFIELDELEMENT_HPP
#define MIXEDVECTORFIELDELEMENT_HPP

#include "MixedScalarElement.hpp"

template<typename TT = double, int sDim = SpaceDimension, int tDim = TimeDimension>
class MixedVectorFieldElementT : public MixedScalarElementT<TT, sDim, tDim> {
protected:
  int dim{};

  MixedVectorFieldElementT(const VectorMatrixBase &base, const Cell &c, int face,
                           const vector<ShapeValues<TT>> &faceValues);
public:
  MixedVectorFieldElementT(const VectorMatrixBase &base, const Cell &c);

  MixedVectorFieldElementT(const VectorMatrixBase &base,
                           const BasicElementT<TT, sDim, tDim> &baseElement);

  int Dim() const { return dim; }

  VectorFieldComponentT<TT> VectorComponentValue(int q, int n, int i, int k) const;

  VectorFieldComponentT<TT> VectorComponentValue(int q, int n, const ShapeId &iter, int k) const;

  VectorFieldComponentT<TT> VectorComponentValue(const PointT<TT, sDim, tDim> &z, int n, int i,
                                                 int k) const;

  VectorFieldComponentT<TT> VectorComponentValue(const PointT<TT, sDim, tDim> &z, int n,
                                                 const ShapeId &iter, int k) const;

  VectorFieldT<TT, sDim> VectorValue(int q, const Vector &u, int n) const;

  VectorFieldT<TT, sDim> VectorValue(const PointT<TT, sDim, tDim> &z, const Vector &u, int n) const;

  TensorRowT<TT> VectorRowGradient(int q, int n, int i, int k) const;

  TensorRowT<TT> VectorRowGradient(int q, int n, const ShapeId &iter, int k) const;

  TensorRowT<TT> VectorRowGradient(const PointT<TT, sDim, tDim> &z, int n, int i, int k) const;


  TensorRowT<TT> VectorRowGradient(const PointT<TT, sDim, tDim> &z, int n, const ShapeId &iter,
                                   int k) const;

  TensorT<TT, sDim> VectorGradient(int q, const Vector &u, int n) const;

  TensorT<TT, sDim> VectorGradient(const PointT<TT, sDim, tDim> &z, const Vector &u, int n) const;

  TT Divergence(int q, int n, int i, int k) const;

  TT Divergence(int q, int n, const ShapeId &iter, int k) const;

  TT Divergence(const PointT<TT, sDim, tDim> &z, int n, int i, int k) const;

  TT Divergence(const PointT<TT, sDim, tDim> &z, int n, const ShapeId &iter, int k) const;

  TT Divergence(int q, const Vector &u, int n) const;

  TT Divergence(const PointT<TT, sDim, tDim> &z, const Vector &u, int n) const;
};

using MixedVectorFieldElement = MixedVectorFieldElementT<>;

template<int sDim = SpaceDimension, int tDim = TimeDimension>
using MixedVectorFieldVectorAccessT =
    MixedVectorAccessEquallyDistributed<MixedVectorFieldElementT<double, sDim, tDim>>;

using MixedVectorFieldVectorAccess = MixedVectorFieldVectorAccessT<>;

template<int sDim = SpaceDimension, int tDim = TimeDimension>
using MixedVectorFieldMatrixAccessT =
    MixedMatrixAccessEquallyDistributed<MixedVectorFieldElementT<double, sDim, tDim>>;

using MixedVectorFieldMatrixAccess = MixedVectorFieldMatrixAccessT<>;

#ifdef BUILD_IA

using IAMixedVectorFieldElement =
    MixedVectorFieldElementT<IAInterval, SpaceDimension, TimeDimension>;

#endif

template<typename TT = double, int sDim = SpaceDimension, int tDim = TimeDimension>
class MixedVectorFieldFaceElementT :
    public IMixedFaceElementT<TT, sDim, tDim>,
    public MixedVectorFieldElementT<TT, sDim, tDim> {
public:
  MixedVectorFieldFaceElementT(const VectorMatrixBase &base, const Cell &c, int face);
};

using MixedVectorFieldFaceElement = MixedVectorFieldFaceElementT<>;

#ifdef BUILD_IA

using IAMixedVectorFieldFaceElement =
    MixedVectorFieldFaceElementT<IAInterval, SpaceDimension, TimeDimension>;

#endif

#endif // MIXEDVECTORFIELDELEMENT_HPP
