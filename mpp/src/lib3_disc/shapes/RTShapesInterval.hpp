#ifndef RTSHAPESINTERVAL_HPP
#define RTSHAPESINTERVAL_HPP

#include "RTShape.hpp"

template<typename T = double, int sDim = SpaceDimension, int tDim = TimeDimension>
class RT0IntT : public RTShapeT<T, sDim, tDim> {
public:
  const std::string Name() const override { return "RT0Int"; }

  VectorFieldT<T, sDim> LocalVector(const PointT<T, sDim, tDim> &z, int i) const override {
    switch (i) {
    case 0:
      return VectorFieldT<T, sDim>(z[0] - 1, T(0.0), T(0.0));
    case 1:
      return VectorFieldT<T, sDim>(z[0], T(0.0), T(0.0));
    default:
      THROW("Not implemented")
    }
  }

  T LocalDiv(const PointT<T, sDim, tDim> &z, int i) const override { return T(1.0); }

  void NodalPoints(const Cell &c, vector<PointT<T, sDim, tDim>> &z) const override {
    z.resize(c.Faces());
    for (int i = 0; i < c.Faces(); ++i)
      z[i] = c.Face(i);
  }

  explicit RT0IntT() : RTShapeT<T, sDim, tDim>(2) {}
};

typedef RT0IntT<> RT0Int;

#endif // RTSHAPESINTERVAL_HPP
