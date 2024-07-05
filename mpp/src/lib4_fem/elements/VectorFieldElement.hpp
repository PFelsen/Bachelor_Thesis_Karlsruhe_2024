#ifndef VECTORFIELDELEMENT_H
#define VECTORFIELDELEMENT_H

#include "ScalarElement.hpp"

template<typename TT = double, int sDim = SpaceDimension, int tDim = TimeDimension>
class VectorFieldElementT : public ScalarElementT<TT, sDim, tDim> {
protected:
  int dim;

  VectorFieldElementT(const VectorMatrixBase &base, const Cell &c, int face,
                      const ShapeValues<TT> &faceValues);
public:
  VectorFieldElementT(const VectorMatrixBase &base, const Cell &c);

  VectorFieldElementT(const VectorMatrixBase &base,
                      const BasicElementT<TT, sDim, tDim> &baseElement);

  VectorFieldComponentT<TT> VectorComponentValue(int q, int i, int k) const;

  VectorFieldT<TT, sDim, tDim> VectorValue(int q, const Vector &u) const;

  VectorFieldT<TT, sDim, tDim> VectorValue(const PointT<TT, sDim, tDim> &z, const Vector &u) const;

  VectorFieldT<TT, sDim, tDim> VectorValue(const PointT<TT, sDim, tDim> &z, int i) const;

  TensorRowT<TT, sDim> VectorRowGradient(int q, int i, int k) const;

  TensorT<TT, sDim> VectorGradient(int q, const Vector &u) const;

  TensorT<TT, sDim> VectorGradient(const PointT<TT, sDim, tDim> &z, const Vector &u) const;

  TT Divergence(int q, int i, int k) const;

  TT Divergence(int q, const Vector &u) const;

  int Dim() const { return dim; }
};

using VectorFieldElement = VectorFieldElementT<>;

#ifdef BUILD_IA

using IAVectorFieldElement = VectorFieldElementT<IAInterval, SpaceDimension, TimeDimension>;

#endif

template<typename TT = double, int sDim = SpaceDimension, int tDim = TimeDimension>
class VectorFieldFaceElementT :
    public IFaceElementT<TT, sDim, tDim>,
    public VectorFieldElementT<TT, sDim, tDim> {
public:
  VectorFieldFaceElementT(const VectorMatrixBase &base, const Cell &c, int face);
};

using VectorFieldFaceElement = VectorFieldFaceElementT<>;

#ifdef BUILD_IA

using IAVectorFieldFaceElement = VectorFieldFaceElementT<IAInterval, SpaceDimension, TimeDimension>;

#endif

#endif // VECTORFIELDELEMENT_H
