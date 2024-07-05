#ifndef DGELEMENT_HPP
#define DGELEMENT_HPP

#include "Element.hpp"

template<typename TT = double, int sDim = SpaceDimension, int tDim = TimeDimension>
class DGElementT : public ElementT<TT, sDim, tDim> {
  void init();
protected:
  const ShapeValues<TT> &value;
  ShapeValues<VectorFieldT<TT, sDim>> gradient;

  vector<PointT<TT, sDim, tDim>> nodalPoints;


  DGElementT(const VectorMatrixBase &base, const Cell &c, int face,
             const ShapeValues<TT> &faceValues);
public:
  DGElementT(const VectorMatrixBase &base, const Cell &c);

  DGElementT(const VectorMatrixBase &base, const BasicElementT<TT, sDim, tDim> &baseElement);

  int shape_size() const;

  int NodalPoints() const override;

  const PointT<TT, sDim, tDim> &NodalPoint(int i) const override;

  TT Value(int q, int i) const;

  TT Value(int q, const Vector &u) const;

  // Used in CDD Project
  TT Value(int q, const Vector &u, int k) const;

  VectorFieldT<TT, sDim> Derivative(int q, int i) const;

  VectorFieldT<TT, sDim> Derivative(int q, const Vector &u) const;
};

using DGElement = DGElementT<>;

#ifdef BUILD_IA

using IADGElement = DGElementT<IAInterval, SpaceDimension, TimeDimension>;

#endif

template<typename TT = double, int sDim = SpaceDimension, int tDim = TimeDimension>
class DGFaceElementT : public IFaceElementT<TT, sDim, tDim>, public DGElementT<TT, sDim, tDim> {
public:
  DGFaceElementT(const VectorMatrixBase &base, const Cell &c, int face);

  int findQPointID(const DGFaceElementT<TT, sDim, tDim> &otherFaceElem,
                   const PointT<TT, sDim, tDim> &Qf_c) const;
};

typedef DGFaceElementT<> DGFaceElement;

#ifdef BUILD_IA

using IADGFaceElement = DGFaceElementT<IAInterval, SpaceDimension, TimeDimension>;

#endif

#endif // DGELEMENT_HPP
