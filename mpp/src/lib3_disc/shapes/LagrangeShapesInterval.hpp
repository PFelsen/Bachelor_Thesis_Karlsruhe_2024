#ifndef LAGRANGESHAPESINTERVAL_HPP
#define LAGRANGESHAPESINTERVAL_HPP

#include "LagrangeShape.hpp"

template<uint32_t degree, typename T = double, int sDim = SpaceDimension, int tDim = TimeDimension>
class IntervalShape : public LagrangeShapeT<T, sDim, tDim> {
private:
  static constexpr auto functions = mpp::shape::polynomial<degree, T>;
  static constexpr auto Dfunctions = mpp::shape::polynomialDerivative<degree, T>;
public:
  IntervalShape() : LagrangeShapeT<T, sDim, tDim>(degree + 1) {}

  void NodalPoints(const Cell &c, vector<PointT<T, sDim, tDim>> &z) const override final {
    z.resize(this->numNodalPoints);
    LagrangeNodalPointsInterval(c, z, degree);
  }

  T operator()(const PointT<T, sDim, tDim> &z, int i) const override final {
    return functions(z[0], i);
  }

  VectorFieldT<T, sDim> LocalGradient(const PointT<T, sDim, tDim> &z, int i) const override {
    return VectorFieldT<T, sDim>(Dfunctions(z[0], i), T(0.0));
  }

  const std::string Name() const override { return std::string{'P', '0' + degree, 'I', 'n', 't'}; }
};

#endif // LAGRANGESHAPESINTERVAL_HPP
