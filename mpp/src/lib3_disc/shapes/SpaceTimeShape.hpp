#ifndef SPACETIME_SPACETIMESHAPE_HPP
#define SPACETIME_SPACETIMESHAPE_HPP

#include "Shape.hpp"

template<typename T = double, int sDim = SpaceDimension, int tDim = TimeDimension>
class SpaceTimeShapeT : public ShapeT<T, sDim, tDim> {
protected:
  const ShapeT<T, sDim, tDim> &SS;
  const ShapeT<T, sDim, tDim> &TS;

  ShapeValues<Scalar> dt_values;
  ShapeValues<VectorFieldT<T, sDim>> div_values;
  ShapeValues<VectorFieldT<T, sDim, tDim>> gradients;
public:
  using ShapeT<T, sDim, tDim>::operator();

  SpaceTimeShapeT(const ShapeT<T, sDim, tDim> &ss, const ShapeT<T, sDim, tDim> &ts) :
      ShapeT<T, sDim, tDim>(ss.size() * ts.size()), SS(ss), TS(ts) {}

  T operator()(const PointT<T, sDim, tDim> &z, int i) const override {
    int ti = i / SS.size();
    int si = i % SS.size();
    return SS(z, si) * TS(Point(z.t()), ti);
  }

  void NodalPoints(const Cell &c, vector<PointT<T, sDim, tDim>> &z) const override {
    z.resize(this->size());
    std::vector<Point> spoints;
    SS.NodalPoints(c.SpaceCell(), spoints);
    std::vector<Point> tpoints;
    TS.NodalPoints(c.TimeCell(), tpoints);
    for (int ti = 0; ti < tpoints.size(); ti++) {
      for (int si = 0; si < spoints.size(); si++) {
        int index = ti * SS.size() + si;
        z[index] = spoints[si].CopyWithT(tpoints[ti][0]);
      }
    }
  }

  VectorFieldT<T, sDim, tDim> LocalSpaceTimeGradient(int q, int i) const { return gradients[q][i]; }

  void fillValues(const vector<PointT<T, sDim, tDim>> &localQPoints) override {
    this->localValues.resize(localQPoints.size(), this->numNodalPoints);
    dt_values.resize(localQPoints.size(), this->numNodalPoints);
    div_values.resize(localQPoints.size(), this->numNodalPoints);
    gradients.resize(localQPoints.size(), this->numNodalPoints);

    for (int q = 0; q < localQPoints.size(); ++q) {
      PointT<T, sDim, tDim> tp(localQPoints[q].t());
      for (int ti = 0; ti < TS.size(); ti++) {
        for (int si = 0; si < SS.size(); si++) {
          int i = ti * SS.size() + si;
          this->localValues[q][i] = (*this)(localQPoints[q], i);
          gradients[q][i] =
              VectorFieldT<T, sDim, tDim>(SS.LocalGradient(localQPoints[q], si) * TS(tp, ti));
          gradients[q][i].t() = SS(localQPoints[q], si) * TS.LocalGradient(tp, ti)[0];
        }
      }
    }
  };

  void fillFaceValues(const vector<vector<PointT<T, sDim, tDim>>> &localFaceQPoints) override{
      Warning("fillFaceValues called, u should implement it...")};

  const ShapeT<T, sDim, tDim> &GetSpaceShape() const { return SS; }

  const ShapeT<T, sDim, tDim> &GetTimeShape() const { return TS; }

  std::string NameS() const {
    std::stringstream ss;
    ss << "SpaceTimeShape[SpaceShape=" << SS.Name() << ",TimeShape=" << TS.Name() << "]";
    return ss.str();
  }

  const std::string Name() const override { return "SpaceTimeShape"; }
};

using SpaceTimeShape = SpaceTimeShapeT<>;


#endif // SPACETIME_SPACETIMESHAPE_HPP
