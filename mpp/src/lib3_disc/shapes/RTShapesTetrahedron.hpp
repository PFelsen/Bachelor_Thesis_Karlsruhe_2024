#ifndef RTSHAPESTETRAHEDRON_HPP
#define RTSHAPESTETRAHEDRON_HPP

#include "RTShape.hpp"

template<typename T = double, int sDim = SpaceDimension, int tDim = TimeDimension>
class RT0TetT : public RTShapeT<T, sDim, tDim> {
public:
  const std::string Name() const override { return "RT0Tet"; }

  VectorFieldT<T, sDim> LocalVector(const PointT<T, sDim, tDim> &z,
                                    int i) const override{THROW("Not implemented")}

  T LocalDiv(const PointT<T, sDim, tDim> &z, int i) const override {
    THROW("Not implemented")
  }

  void NodalPoints(const Cell &c, vector<PointT<T, sDim, tDim>> &z) const override {
    z.resize(c.Faces());
    for (int i = 0; i < c.Faces(); ++i)
      z[i] = c.Face(i);
  }

  explicit RT0TetT() : RTShapeT<T, sDim, tDim>(4) {}
};

typedef RT0TetT<> RT0Tet;

#endif // RTSHAPESTETRAHEDRON_HPP
