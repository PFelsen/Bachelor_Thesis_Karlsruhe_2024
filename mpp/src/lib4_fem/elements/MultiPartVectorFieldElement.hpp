#ifndef MULTIPARTVECTORFIELDELEMENT_HPP
#define MULTIPARTVECTORFIELDELEMENT_HPP

#include "MultiPartScalarElement.hpp"

template<typename TT = double, int sDim = SpaceDimension, int tDim = TimeDimension>
class MultiPartVectorFieldElementT : public MultiPartScalarElementT<TT, sDim, tDim> {
protected:
  int dim;
public:
  MultiPartVectorFieldElementT(const VectorMatrixBase &base, const Cell &c, int index = 0);

  MultiPartVectorFieldElementT(const VectorMatrixBase &base,
                               const BasicElementT<TT, sDim, tDim> &baseElement, int index = 0);

  VectorFieldComponentT<TT> VectorComponentValue(int q, int i, int k) const;

  VectorFieldT<TT, sDim, tDim> VectorValue(int q, const Vector &u) const;

  VectorFieldT<TT, sDim, tDim> VectorValue(const PointT<TT, sDim, tDim> &z, const Vector &u) const;

  VectorFieldT<TT, sDim, tDim> VectorValue(const PointT<TT, sDim, tDim> &z, int i) const;

  TensorRowT<TT, sDim> VectorRowGradient(int q, int i, int k) const;

  TensorT<TT, sDim> VectorGradient(int q, const Vector &u) const;

  TensorT<TT, sDim> VectorGradient(const PointT<TT, sDim, tDim> &z, const Vector &u) const;

  int Dim() const { return dim; }
};

using MultiPartVectorFieldElement = MultiPartVectorFieldElementT<>;

#ifdef BUILD_IA

using IAMultiPartVectorFieldElement =
    MultiPartVectorFieldElementT<IAInterval, SpaceDimension, TimeDimension>;

#endif

#endif // MULTIPARTVECTORFIELDELEMENT_HPP
