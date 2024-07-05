#ifndef RTELEMENT_HPP
#define RTELEMENT_HPP

#include "Element.hpp"

template<typename TT = double, int sDim = SpaceDimension, int tDim = TimeDimension>
class RTElementT : public ElementT<TT, sDim, tDim> {
protected:
  int dim;
  vector<PointT<TT, sDim, tDim>> normal;
  vector<double> sign;
  ShapeValues<VectorFieldT<TT, sDim>> vShapeValues;
  ShapeValues<TT> divShapeValues;
  int nodalPointsOnFacesTotal;

  void init(const VectorMatrixBase &base);

  void initValues(const VectorMatrixBase &base);

  RTElementT(const VectorMatrixBase &base, const Cell &c, int face, int n);
public:
  RTElementT(const VectorMatrixBase &base, const Cell &c);

  RTElementT(const VectorMatrixBase &base, const BasicElementT<TT, sDim, tDim> &baseElement);

  const VectorFieldT<TT, sDim> &VelocityField(int q, int i, int k) const {
    return VelocityField(q, ShapeId{indexingShape(i, k)});
  }

  const VectorFieldT<TT, sDim> &VelocityField(int q, const ShapeId &iter) const {
    return vShapeValues[q][iter.id];
  }

  VectorFieldT<TT, sDim> VelocityField(const PointT<TT, sDim, tDim> &z, int i, int k) const;

  VectorFieldT<TT, sDim> VelocityField(const PointT<TT, sDim, tDim> &z, const ShapeId &iter) const;

  VectorFieldT<TT, sDim> VelocityField(int q, const Vector &u) const;

  VectorFieldT<TT, sDim> VelocityField(int q, const Vector &u, int m) const;

  VectorFieldT<TT, sDim> VelocityField(const PointT<TT, sDim, tDim> &z, const Vector &u) const;

  VectorFieldT<TT, sDim> VelocityField(const PointT<TT, sDim, tDim> &z, const Vector &u,
                                       int m) const;

  TT VelocityFieldDivergence(int q, int i, int k) const {
    return VelocityFieldDivergence(q, ShapeId{indexingShape(i, k)});
  }

  TT VelocityFieldDivergence(int q, const ShapeId &iter) const {
    return divShapeValues[q][iter.id];
  }

  TT VelocityFieldDivergence(const PointT<TT, sDim, tDim> &z, int i, int k) const;

  TT VelocityFieldDivergence(const PointT<TT, sDim, tDim> &z, const ShapeId &iter) const;

  TT VelocityFieldDivergence(int q, const Vector &u) const;

  TT VelocityFieldDivergence(int q, const Vector &u, int m) const;

  TT VelocityFieldDivergence(const PointT<TT, sDim, tDim> &z, const Vector &u) const;

  TT VelocityFieldDivergence(const PointT<TT, sDim, tDim> &z, const Vector &u, int m) const;

  PointT<TT, sDim, tDim> OuterNormal(int face) const { return normal[face]; }

  PointT<TT, sDim, tDim> OrientedNormal(int face) const { return sign[face] * normal[face]; }

  int get_maxk(int i) const {
    if (i < nodalPointsOnFacesTotal) return 1;
    return 2;
  }

  int indexing_k(int i, int k, int m) const {
    if (i < nodalPointsOnFacesTotal) return m;
    return 2 * m + k;
  }
protected:
  int indexingShape(int i, int k) const {
    if (i < nodalPointsOnFacesTotal) return i;
    return 2 * i - nodalPointsOnFacesTotal + k;
  }

  friend class TestRTElement;
};

using RTElement = RTElementT<>;

#ifdef BUILD_IA

using IARTElement = RTElementT<IAInterval, SpaceDimension, TimeDimension>;

#endif

template<int sDim = SpaceDimension, int tDim = TimeDimension>
using RTVectorAccessT = VectorAccessGeneral<RTElementT<double, sDim, tDim>>;

using RTVectorAccess = RTVectorAccessT<>;

template<int sDim = SpaceDimension, int tDim = TimeDimension>
using RTMatrixAccessT = MatrixAccessGeneral<RTElementT<double, sDim, tDim>>;

using RTMatrixAccess = RTMatrixAccessT<>;

template<typename TT = double, int sDim = SpaceDimension, int tDim = TimeDimension>
class RTFaceElementT : public RTElementT<TT, sDim, tDim> {
public:
  RTFaceElementT(const VectorMatrixBase &base, const Cell &c, int face, int n = 0);

  PointT<TT, sDim, tDim> Normal() const { return this->OuterNormal(this->FaceID()); }
};

using RTFaceElement = RTFaceElementT<>;

#ifdef BUILD_IA

using IARTFaceElement = RTFaceElementT<IAInterval, SpaceDimension, TimeDimension>;

#endif

#endif // RTELEMENT_HPP
