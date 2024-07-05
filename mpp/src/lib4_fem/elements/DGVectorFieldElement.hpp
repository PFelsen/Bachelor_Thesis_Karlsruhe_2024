#ifndef DGVECTORFIELDELEMENT_HPP
#define DGVECTORFIELDELEMENT_HPP

#include "DGElement.hpp"

template<typename TT = double, int sDim = SpaceDimension, int tDim = TimeDimension>
class DGVectorFieldElementT : public DGElementT<TT, sDim, tDim> {
protected:
  int dim{};

  DGVectorFieldElementT(const VectorMatrixBase &base, const Cell &c, int face,
                        const ShapeValues<TT> &faceValues);
public:
  DGVectorFieldElementT(const VectorMatrixBase &base, const Cell &c);

  DGVectorFieldElementT(const VectorMatrixBase &base,
                        const BasicElementT<TT, sDim, tDim> &baseElement);

  VectorFieldComponentT<TT> VectorComponentValue(int q, int i, int k) const;


  VectorFieldT<TT, sDim> VectorValue(int q, const Vector &u) const;

  VectorFieldT<TT, sDim> VectorValue(const PointT<TT, sDim, tDim> &z, const Vector &u) const;

  TensorRowT<TT, sDim> VectorRowGradient(int q, int i, int k) const;

  TensorT<TT, sDim> VectorGradient(int q, const Vector &u) const;

  TensorT<TT, sDim> VectorGradient(const PointT<TT, sDim, tDim> &z, const Vector &u) const;

  TT Divergence(int q, int i, int k) const;

  TT Divergence(int q, const Vector &u) const;

  int Dim() const { return dim; }
};

using DGVectorFieldElement = DGVectorFieldElementT<>;

#ifdef BUILD_IA

using IADGVectorFieldElement = DGVectorFieldElementT<IAInterval, SpaceDimension, TimeDimension>;

#endif

template<typename TT = double, int sDim = SpaceDimension, int tDim = TimeDimension>
class DGVectorFieldFaceElementT :
    public IFaceElementT<TT, sDim, tDim>,
    public DGVectorFieldElementT<TT, sDim, tDim> {
public:
  DGVectorFieldFaceElementT(const VectorMatrixBase &base, const Cell &c, int face);


  int findQPointID(const DGVectorFieldFaceElementT<TT, sDim, tDim> &otherFaceElem,
                   const PointT<TT, sDim, tDim> &Qf_c) const;
};

using DGVectorFieldFaceElement = DGVectorFieldFaceElementT<>;

#ifdef BUILD_IA

using IADGVectorFieldFaceElement =
    DGVectorFieldFaceElementT<IAInterval, SpaceDimension, TimeDimension>;

#endif

#endif // DGVECTORFIELDELEMENT_HPP
