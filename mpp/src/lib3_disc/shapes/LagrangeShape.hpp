#ifndef LAGRANGESHAPE_H
#define LAGRANGESHAPE_H

#include "LagrangeNodalPoints.hpp"
#include "PolynomialSystem.hpp"
#include "Polynomials.hpp"
#include "Shape.hpp"

template<typename T = double, int sDim = SpaceDimension, int tDim = TimeDimension>
class LagrangeShapeT : public ShapeT<T, sDim, tDim> {
public:
  LagrangeShapeT(int numNodalPoints) : ShapeT<T, sDim, tDim>(numNodalPoints) {}

  void fillValues(const vector<PointT<T, sDim, tDim>> &localQuadraturePoints) override {
    this->localValues.resize(localQuadraturePoints.size(), this->numNodalPoints);
    this->localGradient.resize(localQuadraturePoints.size(), this->numNodalPoints);
    for (int q = 0; q < localQuadraturePoints.size(); ++q)
      for (int i = 0; i < this->numNodalPoints; ++i) {
        this->localValues[q][i] = (*this)(localQuadraturePoints[q], i);
        this->localGradient[q][i] = this->LocalGradient(localQuadraturePoints[q], i);
      }
  }

  void
  fillFaceValues(const vector<vector<PointT<T, sDim, tDim>>> &localFaceQuadraturePoints) override {
    this->localFaceValues.resize(localFaceQuadraturePoints.size(),
                                 localFaceQuadraturePoints[0].size(), this->numNodalPoints);
    for (int face = 0; face < localFaceQuadraturePoints.size(); ++face)
      for (int q = 0; q < localFaceQuadraturePoints[face].size(); ++q)
        for (int i = 0; i < this->numNodalPoints; ++i) {
          this->localFaceValues[face][q][i] = (*this)(localFaceQuadraturePoints[face][q], i);
        }
  }
};


#endif // LAGRANGESHAPE_H
