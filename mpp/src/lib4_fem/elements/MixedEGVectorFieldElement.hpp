#ifndef MIXEDEGVFELEMENT_HPP
#define MIXEDEGVFELEMENT_HPP

#include "MixedVectorFieldElement.hpp"

template<typename TT = double, int sDim = SpaceDimension, int tDim = TimeDimension>
class MixedEGVectorFieldElementT : public MixedVectorFieldElementT<TT, sDim, tDim> {
protected:
  MixedEGVectorFieldElementT(const VectorMatrixBase &base, const Cell &c, int face,
                             const vector<ShapeValues<TT>> &faceValues);
public:
  MixedEGVectorFieldElementT(const VectorMatrixBase &base, const Cell &c) :
      MixedVectorFieldElementT<TT, sDim, tDim>(base, c) {}

  constexpr int PenaltyIndex() const { return 1; }

  int PenaltySize() const { return this->size(1); }

  constexpr int ValueIndex() const { return 0; }

  int ValueSize() const { return this->size(0); }

  VectorFieldT<TT, sDim> PenaltyVectorValue(int q) const;

  VectorFieldT<TT, sDim> PenaltyVectorValue(int q, const Vector &u) const;

  constexpr TensorT<TT, sDim> PenaltyVectorGradient(int q) const { return One; };

  TensorT<TT, sDim> PenaltyVectorGradient(int q, const Vector &u) const;

  VectorFieldComponentT<TT> VectorComponentValue(int q, int i, int k) const;

  VectorFieldT<TT, sDim> VectorValue(int q, const Vector &u) const;

  VectorFieldT<TT, sDim> VectorValue(const PointT<TT, sDim, tDim> &z, const Vector &u) const;

  TensorRowT<TT, sDim> VectorRowGradient(int q, int i, int k) const;

  TensorT<TT, sDim> VectorGradient(int q, const Vector &u) const;

  TensorT<TT, sDim> VectorGradient(const PointT<TT, sDim, tDim> &z, const Vector &u) const;
};

using MixedEGVectorFieldElement = MixedEGVectorFieldElementT<>;

#ifdef BUILD_IA

using IAMixedEGVectorFieldElement =
    MixedEGVectorFieldElementT<IAInterval, SpaceDimension, TimeDimension>;

#endif

template<typename TT = double, int sDim = SpaceDimension, int tDim = TimeDimension>
class MixedEGVectorFieldFaceElementT :
    public IMixedFaceElementT<TT, sDim, tDim>,
    public MixedEGVectorFieldElementT<TT, sDim, tDim> {
public:
  MixedEGVectorFieldFaceElementT(const VectorMatrixBase &base, const Cell &c, int face);

  int FindQPointID(const MixedEGVectorFieldFaceElementT<TT, sDim, tDim> &otherFaceElem,
                   const PointT<TT, sDim, tDim> &Qf_c) const;
};

using MixedEGVectorFieldFaceElement = MixedEGVectorFieldFaceElementT<>;

#ifdef BUILD_IA

using IAMixedEGVectorFieldFaceElement =
    MixedEGVectorFieldFaceElementT<IAInterval, SpaceDimension, TimeDimension>;

#endif

#endif // MIXEDEGELEMENT_HPP
