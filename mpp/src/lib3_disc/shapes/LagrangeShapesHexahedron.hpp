#ifndef LAGRANGESHAPESHEXAHEDRON_HPP
#define LAGRANGESHAPESHEXAHEDRON_HPP

#include "LagrangeShape.hpp"

template<uint32_t degree, typename T = double, int sDim = SpaceDimension, int tDim = TimeDimension>
class HexahedronShape : public LagrangeShapeT<T, sDim, tDim> {
private:
  static constexpr auto functions = mpp::shape::polynomial<degree, T>;
  static constexpr auto Dfunctions = mpp::shape::polynomialDerivative<degree, T>;
public:
  HexahedronShape() : LagrangeShapeT<T, sDim, tDim>(pow(degree + 1, 3)) {}

  void NodalPoints(const Cell &c, vector<PointT<T, sDim, tDim>> &z) const override final {
    z.resize(this->numNodalPoints);
    LagrangeNodalPointsHexahedron(c, z, degree);
  }

  T operator()(const PointT<T, sDim, tDim> &z, int i) const override final {
    int N = degree + 1;
    int iDN = i / N;
    return functions(z[0], i % N) * functions(z[1], iDN % N) * functions(z[2], iDN / N);
  }

  VectorFieldT<T, sDim> LocalGradient(const PointT<T, sDim, tDim> &z, int i) const override {
    int N = degree + 1;
    int l = i % N;
    int iDN = i / N;
    int m = iDN % N;
    int n = iDN / N;
    return VectorFieldT<T, sDim>(Dfunctions(z[0], l) * functions(z[1], m) * functions(z[2], n),
                                 functions(z[0], l) * Dfunctions(z[1], m) * functions(z[2], n),
                                 functions(z[0], l) * functions(z[1], m) * Dfunctions(z[2], n));
  }

  const std::string Name() const override { return std::string{'P', '0' + degree, 'H', 'e', 'x'}; }
};

#endif // LAGRANGESHAPESHEXAHEDRON_HPP
