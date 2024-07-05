#ifndef TAYLORHOODELEMENT_HPP
#define TAYLORHOODELEMENT_HPP

#include "MixedVectorFieldElement.hpp"

template<typename TT = double, int sDim = SpaceDimension, int tDim = TimeDimension>
class TaylorHoodElementT : public MixedVectorFieldElementT<TT, sDim, tDim> {
protected:
  TaylorHoodElementT(const VectorMatrixBase &base, const Cell &c, int face,
                     const vector<ShapeValues<TT>> &faceValues);
public:
  TaylorHoodElementT(const VectorMatrixBase &base, const Cell &c) :
      MixedVectorFieldElementT<TT, sDim, tDim>(base, c) {}

  TaylorHoodElementT(const VectorMatrixBase &base,
                     const BasicElementT<TT, sDim, tDim> &baseElement) :
      MixedVectorFieldElementT<TT, sDim, tDim>(base, baseElement) {}

  constexpr int PressureIndex() const { return 1; }

  int PressureSize() const { return this->size(1); }

  constexpr int VelocityIndex() const { return 0; }

  int VelocitySize() const { return this->size(0); }

  TT PressureValue(int q, int i) const;

  TT PressureValue(int q, const ShapeId &iter) const;

  TT PressureValue(int q, const Vector &u) const;

  const VectorFieldT<TT, sDim> &PressureGradient(int q, int i) const;

  const VectorFieldT<TT, sDim> &PressureGradient(int q, const ShapeId &iter) const;

  VectorFieldT<TT, sDim> PressureGradient(int q, const Vector &u) const;

  VectorFieldComponentT<TT> VelocityFieldComponent(int q, int i, int k) const;

  VectorFieldComponentT<TT> VelocityFieldComponent(int q, const ShapeId &iter, int k) const;

  VectorFieldT<TT, sDim> VelocityField(int q, const Vector &u) const;

  TensorRowT<TT> VelocityFieldRowGradient(int q, int i, int k) const;

  TensorRowT<TT> VelocityFieldRowGradient(int q, const ShapeId &iter, int k) const;

  TensorT<TT, sDim> VelocityFieldGradient(int q, const Vector &u) const;

  TT VelocityFieldDivergence(int q, int i, int k) const;

  TT VelocityFieldDivergence(int q, const ShapeId &iter, int k) const;

  TT VelocityFieldDivergence(int q, const Vector &u) const;
};

using TaylorHoodElement = TaylorHoodElementT<>;

template<int sDim = SpaceDimension, int tDim = TimeDimension>
using TaylorHoodVectorAccessT =
    MixedVectorAccessEquallyDistributed<TaylorHoodElementT<double, sDim, tDim>>;

using TaylorHoodVectorAccess = TaylorHoodVectorAccessT<>;

template<int sDim = SpaceDimension, int tDim = TimeDimension>
using TaylorHoodMatrixAccessT =
    MixedMatrixAccessEquallyDistributed<TaylorHoodElementT<double, sDim, tDim>>;

using TaylorHoodMatrixAccess = TaylorHoodMatrixAccessT<>;

#ifdef BUILD_IA

using IATaylorHoodElement = TaylorHoodElementT<IAInterval, SpaceDimension, TimeDimension>;

#endif

template<typename TT = double, int sDim = SpaceDimension, int tDim = TimeDimension>
class TaylorHoodFaceElementT :
    public IMixedFaceElementT<TT, sDim, tDim>,
    public TaylorHoodElementT<TT, sDim, tDim> {
public:
  TaylorHoodFaceElementT(const VectorMatrixBase &base, const Cell &c, int face);
};

using TaylorHoodFaceElement = TaylorHoodFaceElementT<>;

#ifdef BUILD_IA

using IATaylorHoodFaceElement = TaylorHoodFaceElementT<IAInterval, SpaceDimension, TimeDimension>;

#endif

#endif // TAYLORHOODELEMENT_HPP
