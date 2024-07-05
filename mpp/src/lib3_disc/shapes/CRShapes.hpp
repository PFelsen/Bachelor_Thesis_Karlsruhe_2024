#ifndef CRSHAPES_HPP
#define CRSHAPES_HPP

#include "Shape.hpp"

template<typename T = double, int sDim = SpaceDimension, int tDim = TimeDimension>
class CR1TriT : public ShapeT<T, sDim, tDim> {
public:
  const std::string Name() const override { return "CR1Tri"; }

  T operator()(const PointT<T, sDim, tDim> &z, int i) const override {
    switch (i) {
    case 0:
      return 1.0 - 2 * z[1];
    case 1:
      return 2 * (z[0] + z[1]) - 1.0;
    case 2:
      return 1.0 - 2 * z[0];
    default:
      THROW("Not implemented")
    }
  }

  VectorFieldT<T, sDim> LocalGradient(const PointT<T, sDim, tDim> &z, int i) const override {
    switch (i) {
    case 0:
      return VectorFieldT<T, sDim>(0.0, -2.0);
    case 1:
      return VectorFieldT<T, sDim>(2.0, 2.0);
    case 2:
      return VectorFieldT<T, sDim>(-2.0, 0.0);
    default:
      THROW("Not implemented")
    }
  }

  CR1TriT() : ShapeT<T, sDim, tDim>(3) {}
};

using CR1Tri = CR1TriT<>;

template<typename T, int sDim, int tDim>
ShapeT<T, sDim, tDim> *createCRShape(const CELLTYPE cellType, int degree) {
  switch (cellType) {
  case TRIANGLE:
    switch (degree) {
    case 1:
      return new CR1TriT<T, sDim, tDim>();
    default:
      THROW("Shape not implemented")
    }
  default:
    THROW("Shape not implemented")
  }
}

#endif // CRSHAPES_HPP
