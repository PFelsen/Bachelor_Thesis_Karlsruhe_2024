#ifndef RTSHAPE_H
#define RTSHAPE_H

#include "Shape.hpp"

template<typename T = double, int sDim = SpaceDimension, int tDim = TimeDimension>
class RTShapeT : public ShapeT<T, sDim, tDim> {
public:
  RTShapeT(int numNodalPoints) : ShapeT<T, sDim, tDim>(numNodalPoints) {}

  void fillValues(const vector<PointT<T, sDim, tDim>> &localQuadraturePoints) override {
    this->localVector.resize(localQuadraturePoints.size(), this->numNodalPoints);
    this->localDivs.resize(localQuadraturePoints.size(), this->numNodalPoints);
    for (int q = 0; q < localQuadraturePoints.size(); ++q)
      for (int i = 0; i < this->numNodalPoints; ++i) {
        this->localVector[q][i] = this->LocalVector(localQuadraturePoints[q], i);
        this->localDivs[q][i] = this->LocalDiv(localQuadraturePoints[q], i);
      }
  }

  void
  fillFaceValues(const vector<vector<PointT<T, sDim, tDim>>> &localFaceQuadraturePoints) override {
    // TODO
  }
};

#endif // RTSHAPE_H
