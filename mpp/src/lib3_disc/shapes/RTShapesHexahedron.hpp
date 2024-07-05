#ifndef RTSHAPESHEXAHEDRON_HPP
#define RTSHAPESHEXAHEDRON_HPP

#include "RTShape.hpp"

template<typename T = double, int sDim = SpaceDimension, int tDim = TimeDimension>
class RT0HexT : public RTShapeT<T, sDim, tDim> {
public:
  const std::string Name() const override { return "RT0Hex"; }

  VectorFieldT<T, sDim> LocalVector(const PointT<T, sDim, tDim> &z, int i) const override {
    switch (i) {
    case 0:
      return VectorFieldT<T, sDim>(T(0.0), T(0.0), z[2] - 1.0);
    case 1:
      return VectorFieldT<T, sDim>(T(0.0), z[1] - 1.0, T(0.0));
    case 2:
      return VectorFieldT<T, sDim>(z[0], T(0.0), T(0.0));
    case 3:
      return VectorFieldT<T, sDim>(T(0.0), z[1], T(0.0));
    case 4:
      return VectorFieldT<T, sDim>(z[0] - 1.0, T(0.0), T(0.0));
    case 5:
      return VectorFieldT<T, sDim>(T(0.0), T(0.0), z[2]);
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

  explicit RT0HexT() : RTShapeT<T, sDim, tDim>(6) {}
};

typedef RT0HexT<> RT0Hex;

#endif // RTSHAPESHEXAHEDRON_HPP
