#ifndef LAGRANGESHAPESQUADRILATERAL_HPP
#define LAGRANGESHAPESQUADRILATERAL_HPP

#include "LagrangeShape.hpp"

template<uint32_t degree, typename T = double, int sDim = SpaceDimension, int tDim = TimeDimension>
class QuadraticShape : public LagrangeShapeT<T, sDim, tDim> {
private:
  static constexpr auto functions = mpp::shape::polynomial<degree, T>;
  static constexpr auto Dfunctions = mpp::shape::polynomialDerivative<degree, T>;
public:
  QuadraticShape() : LagrangeShapeT<T, sDim, tDim>(pow(degree + 1, 2)) {}

  void NodalPoints(const Cell &c, vector<PointT<T, sDim, tDim>> &z) const override {
    z.resize(this->numNodalPoints);
    LagrangeNodalPointsQuadrilateral(c, z, degree);
  }

  T operator()(const PointT<T, sDim, tDim> &z, int i) const override {
    int N = degree + 1;
    return functions(z[0], i % N) * functions(z[1], i / N);
  }

  VectorFieldT<T, sDim> LocalGradient(const PointT<T, sDim, tDim> &z, int i) const override {
    int N = degree + 1;
    int m = i % N;
    int n = i / N;
    return VectorFieldT<T, sDim>(Dfunctions(z[0], m) * functions(z[1], n),
                                 functions(z[0], m) * Dfunctions(z[1], n));
  }

  const std::string Name() const override {
    return std::string{'P', '0' + degree, 'Q', 'u', 'a', 'd'};
  }
};
#endif