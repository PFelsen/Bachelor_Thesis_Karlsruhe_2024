#ifndef RTSHAPESTRIANGULAR_HPP
#define RTSHAPESTRIANGULAR_HPP

#include "LagrangeShapesTriangular.hpp"
#include "RTNodalPoints.hpp"
#include "RTShape.hpp"

/*
 * Based on paper by Ervin
 * Computational Bases for RT_k and BDM_k on Triangles
 */

template<uint32_t order, typename T = double, int sDim = SpaceDimension, int tDim = TimeDimension>
class RTTriangle : public RTShapeT<T, sDim, tDim> {
protected:
  static_assert(order < 9, "RTTriangles cannot have a order of more then 6!");

  static constexpr T functions(const T x, const int i) {
    if constexpr (order > 0) return mpp::shape::shapeFunction<order - 1, T>(x, i);
    else THROW("Not implemented!");
  }

  static constexpr T Dfunctions(const T x, const int i) {
    if constexpr (order > 0) return mpp::shape::shapeDerivative<order - 1, T>(x, i);
    else THROW("Not implemented!");
  }

  static constexpr int numNodalPointsFromDegree = (order + 1) * (order + 3);
  static constexpr int numNodelPointsInCell = (order * (order + 1)) / 2;

  VectorFieldT<T, sDim> RT0_VectorField(int i, const PointT<T, sDim, tDim> &z) const {
    switch (i) {
    case 0:
      return VectorFieldT<T, sDim>(z[0], z[1] - 1.0);
    case 1:
      return VectorFieldT<T, sDim>(z[0], z[1]);
    case 2:
      return VectorFieldT<T, sDim>(z[0] - 1.0, z[1]);
    default:
      THROW("Not implemented")
    }
  }

  T RT0_Divergence(int i, const PointT<T, sDim, tDim> &z) const {
    switch (i) {
    case 0:
    case 1:
    case 2:
      return T(2.0);
    default:
      THROW("Not implemented")
    }
  }

  VectorFieldT<T, sDim> Cell_VectorField(int i, const PointT<T, sDim, tDim> &z) const {
    switch (i) {
    case 0:
      return z[1] * VectorFieldT<T, sDim>(z[0], z[1] - 1.0);
    case 1:
      return z[0] * VectorFieldT<T, sDim>(z[0] - 1.0, z[1]);
    default:
      THROW("Not implemented")
    }
  }

  T Cell_Divergence(int i, const PointT<T, sDim, tDim> &z) const {
    switch (i) {
    case 0:
      return 3.0 * z[1] - 1.0;
    case 1:
      return 3.0 * z[0] - 1.0;
    default:
      THROW("Not implemented")
    }
  }

  constexpr T ValueInterval(T x, int i) const {
    return mpp::shape::polynomial<order, T>(this->transformation(x), i);
  }

  constexpr T DerivativeInterval(T x, int i) const {
    return mpp::shape::polynomialDerivative<order, T>(this->transformation(x), i)
           * this->transformationDerivative();
  }

  T ValueTriangle(const PointT<T, sDim, tDim> &z, int i) const {
    T psi[3] = {1 - z[0] - z[1], z[0], z[1]};
    return calculateTriangle<functions, numNodelPointsInCell, T>(psi, i, order - 1);
  }

  VectorFieldT<T, sDim> GradientTriangle(const PointT<T, sDim, tDim> &z, int i) const {
    T psi[3] = {1 - z[0] - z[1], z[0], z[1]};
    return calculateDerivative<functions, Dfunctions, numNodelPointsInCell, sDim, T>(psi, i,
                                                                                     order - 1);
  }

  int nodalPointsOnFaces;
  int nodalPointsInCell;

  T transformation(T x) const { return x + (2 * x - 1.0) / order; }

  T transformationDerivative() const { return 1.0 + T(2.0) / order; }
public:
  RTTriangle() :
      RTShapeT<T, sDim, tDim>(numNodalPointsFromDegree), nodalPointsOnFaces(3 * (order + 1)),
      nodalPointsInCell(numNodelPointsInCell) {}

  void NodalPoints(const Cell &c, vector<PointT<T, sDim, tDim>> &z) const override {
    z.resize(numNodalPointsFromDegree);
    RTNodalPointsTriangle(c, z, order);
  }

  VectorFieldT<T, sDim> LocalVector(const PointT<T, sDim, tDim> &z, int i) const override {
    if constexpr (order == 0) return RT0_VectorField(i, z);

    if (i < nodalPointsOnFaces) {
      int j = i / 3;
      switch (i % 3) {
      case 0:
        return ValueInterval(z[0], j) * RT0_VectorField(0, z);
      case 1:
        return ValueInterval(z[1], j) * RT0_VectorField(1, z);
      case 2:
        return ValueInterval(z[1], order - j) * RT0_VectorField(2, z);
      }
    } else {
      i -= nodalPointsOnFaces;
      return ValueTriangle(z, i / 2) * Cell_VectorField(i % 2, z);
    }
    THROW("") // will never reach this exit, just to prevent warning
  }

  T LocalDiv(const PointT<T, sDim, tDim> &z, int i) const override {
    if (order == 0) return RT0_Divergence(i, z);

    if (i < nodalPointsOnFaces) {
      int j = i / 3;
      switch (i % 3) {
      case 0:
        return ValueInterval(z[0], j) * RT0_Divergence(0, z)
               + DerivativeInterval(z[0], j) * RT0_VectorField(0, z)[0];
      case 1:
        return ValueInterval(z[1], j) * RT0_Divergence(1, z)
               + DerivativeInterval(z[1], j) * RT0_VectorField(1, z)[1];
      case 2:
        return ValueInterval(z[1], order - j) * RT0_Divergence(2, z)
               + DerivativeInterval(z[1], order - j) * RT0_VectorField(2, z)[1];
      }
    } else {
      i -= nodalPointsOnFaces;
      return ValueTriangle(z, i / 2) * Cell_Divergence(i % 2, z)
             + GradientTriangle(z, i / 2) * Cell_VectorField(i % 2, z);
    }
    THROW("") // will never reach this exit, just to prevent warning
  }

  const std::string Name() const override {
    return std::string{'R', 'T', '0' + order, 'T', 'r', 'i'};
  }
};
#endif // RTSHAPESTRIANGULAR_HPP
