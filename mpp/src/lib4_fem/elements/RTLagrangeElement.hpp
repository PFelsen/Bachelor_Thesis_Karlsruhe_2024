#ifndef RTLAGRANGEELEMENT_HPP
#define RTLAGRANGEELEMENT_HPP

#include "Element.hpp"

template<typename TT = double, int sDim = SpaceDimension, int tDim = TimeDimension>
class RTLagrangeElementT : public MixedElementT<TT, sDim, tDim> {
protected:
  int dim;
  vector<PointT<TT, sDim, tDim>> normal;
  vector<double> sign;
  vector<VectorFieldT<TT, sDim>> vectorfield_face;
  const ShapeValues<TT> &pShapeValues;
  ShapeValues<VectorFieldT<TT, sDim>> pShapeGradients;
  ShapeValues<VectorFieldT<TT, sDim>> vShapeValues;
  ShapeValues<TT> divShapeValues;
  int nodalPointsOnFacesTotalRT;
  int velocitySize;
  int lagrangeDegree;

  void init(const VectorMatrixBase &base);

  void initValues(const VectorMatrixBase &base);

  RTLagrangeElementT(const VectorMatrixBase &base, const Cell &c,
                     const ShapeValues<TT> &pShapeValuesFace, int face);
public:
  RTLagrangeElementT(const VectorMatrixBase &base, const Cell &c);

  constexpr int PressureIndex() const { return 1; }

  int PressureSize() const { return this->size(1); }

  constexpr int VelocityIndex() const { return 0; }

  int VelocitySize() const { return this->velocitySize; }

  TT PressureValue(int q, int i) const { return pShapeValues[q][i]; }

  TT PressureValue(int q, const ShapeId &iter) const { return pShapeValues[q][iter.id]; }

  TT PressureValue(int q, const Vector &u, int m = 0) const;

  TT PressureValue(const PointT<TT, sDim, tDim> &z, const Vector &u, int m = 0) const;

  const VectorFieldT<TT, sDim> &PressureGradient(int q, int i) const;

  const VectorFieldT<TT, sDim> &PressureGradient(int q, const ShapeId &iter) const;

  VectorFieldT<TT, sDim> PressureGradient(int q, const Vector &u, int m = 0) const;

  VectorFieldT<TT, sDim> PressureGradient(const PointT<TT, sDim, tDim> &z, const Vector &u,
                                          int m = 0) const;

  const VectorFieldT<TT, sDim> &VelocityField(int q, int i, int k) const {
    return vShapeValues[q][IndexingVelocityShape(i, k)];
  }

  const VectorFieldT<TT, sDim> &VelocityField(int q, const ShapeId &iter) const {
    return vShapeValues[q][iter.id];
  }

  VectorFieldT<TT, sDim> VelocityField(int q, const Vector &u) const;

  VectorFieldT<TT, sDim> VelocityField(int q, const Vector &u, int m) const;

  VectorFieldT<TT, sDim> VelocityField(const PointT<TT, sDim, tDim> &z, const Vector &u) const;

  VectorFieldT<TT, sDim> VelocityField(const PointT<TT, sDim, tDim> &z, const Vector &u,
                                       int m) const;

  VectorFieldT<TT, sDim> CellFlux(const Vector &u) const;

  TT VelocityFieldDivergence(int q, int i, int k) const {
    return divShapeValues[q][IndexingVelocityShape(i, k)];
  }

  TT VelocityFieldDivergence(int q, const ShapeId &iter) const {
    return divShapeValues[q][iter.id];
  }

  TT VelocityFieldDivergence(int q, const Vector &u) const;

  TT VelocityFieldDivergence(int q, const Vector &u, int m) const;

  TT VelocityFieldDivergence(const PointT<TT, sDim, tDim> &z, const Vector &u) const;

  TT VelocityFieldDivergence(const PointT<TT, sDim, tDim> &z, const Vector &u, int m) const;

  TT Sign(int face) const { return sign[face]; }

  PointT<TT, sDim, tDim> OuterNormal(int face) const { return normal[face]; }

  PointT<TT, sDim, tDim> OrientedNormal(int face) const { return sign[face] * normal[face]; }

  int VelocityMaxk(int i) const {
    if (i < nodalPointsOnFacesTotalRT) return 1;
    return 2;
  }

  int IndexingVelocity_k(int i, int k, int m) const {
    if (i < nodalPointsOnFacesTotalRT) return m;
    return 2 * m + k;
  }
protected:
  int IndexingVelocityShape(int i, int k) const {
    if (i < nodalPointsOnFacesTotalRT) return i;
    return 2 * i - nodalPointsOnFacesTotalRT + k;
  }
};

using RTLagrangeElement = RTLagrangeElementT<>;

#ifdef BUILD_IA

using IARTLagrangeElement = RTLagrangeElementT<IAInterval, SpaceDimension, TimeDimension>;

#endif

template<typename TT = double, int sDim = SpaceDimension, int tDim = TimeDimension>
class RTLagrangeFaceElementT :
    public IMixedFaceElementT<TT, sDim, tDim>,
    public RTLagrangeElementT<TT, sDim, tDim> {
public:
  RTLagrangeFaceElementT(const VectorMatrixBase &base, const cell &c, int face) :
      RTLagrangeFaceElementT(base, *c, face){};

  RTLagrangeFaceElementT(const VectorMatrixBase &base, const Cell &c, int face);

  PointT<TT, sDim, tDim> Normal() const { return this->QNormal(0); }

  PointT<TT, sDim, tDim> Normal(int q) const { return this->QNormal(q); }
};

using RTLagrangeFaceElement = RTLagrangeFaceElementT<>;

#ifdef BUILD_IA

using IARTLagrangeFaceElement = RTLagrangeFaceElementT<IAInterval, SpaceDimension, TimeDimension>;

#endif

#endif // RTLAGRANGEELEMENT_HPP
