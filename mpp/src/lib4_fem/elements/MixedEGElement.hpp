#ifndef MIXEDEGELEMENT_HPP
#define MIXEDEGELEMENT_HPP

#include "MixedScalarElement.hpp"

template<typename TT = double, int sDim = SpaceDimension, int tDim = TimeDimension>
class MixedEGElementT : public MixedScalarElementT<TT, sDim, tDim> {
protected:
  MixedEGElementT(const VectorMatrixBase &base, const Cell &c, int face,
                  const vector<ShapeValues<TT>> &faceValues);
public:
  MixedEGElementT(const VectorMatrixBase &base, const Cell &c);

  constexpr int PenaltyIndex() const { return 1; }

  int PenaltySize() const { return this->size(1); }

  constexpr int ValueIndex() const { return 0; }

  int ValueSize() const { return this->size(0); }

  TT PenaltyValue(int q) const;

  TT PenaltyValue(int q, const Vector &u) const;

  VectorFieldT<TT, sDim> PenaltyGradient(int q) const;

  VectorFieldT<TT, sDim> PenaltyGradient(int q, const Vector &u) const;

  TT Value(int q, int i) const;

  TT Value(int q, const Vector &u) const;

  TT Value(const PointT<TT, sDim, tDim> &z, const Vector &u) const;

  VectorFieldT<TT, sDim> Derivative(int q, int i) const;

  VectorFieldT<TT, sDim> Derivative(int q, const Vector &u) const;

  VectorFieldT<TT, sDim> Derivative(const PointT<TT, sDim, tDim> &z, const Vector &u) const;
};

using MixedEGElement = MixedEGElementT<>;

#ifdef BUILD_IA

using IAMixedEGElement = MixedEGElementT<IAInterval, SpaceDimension, TimeDimension>;

#endif

template<typename TT = double, int sDim = SpaceDimension, int tDim = TimeDimension>
class MixedEGFaceElementT :
    public IMixedFaceElementT<TT, sDim, tDim>,
    public MixedEGElementT<TT, sDim, tDim> {
public:
  MixedEGFaceElementT(const VectorMatrixBase &base, const Cell &c, int face);

  int FindQPointID(const MixedEGFaceElementT<TT, sDim, tDim> &otherFaceElem,
                   const PointT<TT, sDim, tDim> &Qf_c) const;
};

using MixedEGFaceElement = MixedEGFaceElementT<>;

#ifdef BUILD_IA

using IAMixedEGFaceElement = MixedEGFaceElementT<IAInterval, SpaceDimension, TimeDimension>;

#endif

#endif // MIXEDEGELEMENT_HPP
