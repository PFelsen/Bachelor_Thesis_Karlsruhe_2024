#ifndef LAGRANGESHAPESTRIANGULAR_HPP
#define LAGRANGESHAPESTRIANGULAR_HPP

#include <string>
#include "LagrangeShape.hpp"

template<auto functions, uint32_t nodalPoints, typename T>
constexpr T calculateTriangle(const T *const psi, const int i, const int currentPolynomialDegree) {

  int index = 0;
  int degree = currentPolynomialDegree;

  for (int callCount = 0; index != nodalPoints; callCount++) {
    int functionNumber1 = degree + callCount;
    int functionNumber2 = callCount;
    const int functionNumber3 = callCount;

    if (degree == 0) {
      if (index == i) {
        return functions(psi[0], functionNumber1) * functions(psi[1], functionNumber2)
               * functions(psi[2], functionNumber3);
      }
      index++;
    }
    for (int n = 0; n < degree; n++) {
      for (int k = 0; k < 3; k++) {
        if (index == i) {
          const auto second = (1 + k) % 3;
          const auto third = (2 + k) % 3;
          return functions(psi[k], functionNumber1) * functions(psi[second], functionNumber2)
                 * functions(psi[third], functionNumber3);
        }
        index++;
      }
      functionNumber2++;
      functionNumber1--;
    }
    degree -= 3;
  }
  THROW("Not implemented")
}

template<auto functions, auto Dfunctions, int sDim, typename T>
constexpr VectorFieldT<T, sDim> calculateSingleField(const T *const psi, const int functionIds[3],
                                                     const int k) {
  const auto firstID = (3 - k) % 3;
  const auto secondID = (4 - k) % 3;
  const auto thirdID = (5 - k) % 3;

  const auto minuendFunction = functions(psi[0], functionIds[firstID]);
  const auto lowerFunction = functions(psi[1], functionIds[secondID]);
  const auto upperFunction = functions(psi[2], functionIds[thirdID]);

  const auto dSubtract = Dfunctions(psi[0], functionIds[firstID]);
  const auto dUpper = Dfunctions(psi[1], functionIds[secondID]);
  const auto dLower = Dfunctions(psi[2], functionIds[thirdID]);

  const auto subtraction = dSubtract * lowerFunction * upperFunction;

  return VectorFieldT<T, sDim>(minuendFunction * dUpper * upperFunction - subtraction,
                               minuendFunction * dLower * lowerFunction - subtraction);
}

template<auto functions, auto Dfunctions, uint32_t nodalPoints, int sDim, typename T>
constexpr VectorFieldT<T, sDim> calculateDerivative(const T *const psi, const int i,
                                                    const int currentPolynomialDegree) {

  int index = 0;
  int degree = currentPolynomialDegree;
  constexpr auto calculateSingle = calculateSingleField<functions, Dfunctions, sDim, T>;
  int functionIds[3];
  for (int callCount = 0; index != nodalPoints; callCount++) {
    functionIds[0] = degree + callCount;
    functionIds[1] = callCount;
    functionIds[2] = callCount;

    if (degree == 0) {
      if (index == i) { return calculateSingle(psi, functionIds, 0); }
      index++;
    }
    for (int n = 0; n < degree; n++) {
      for (int k = 0; k < 3; k++) {
        if (index == i) { return calculateSingle(psi, functionIds, k); }
        index++;
      }
      functionIds[0]--;
      functionIds[1]++;
    }
    degree -= 3;
  }
  THROW("Not implemented")
}

template<uint32_t degree, typename T = double, int sDim = SpaceDimension, int tDim = TimeDimension>
class TriangularShape : public LagrangeShapeT<T, sDim, tDim> {
private:
  static constexpr auto functions = mpp::shape::shapeFunction<degree, T>;
  static constexpr auto functionDerivatives = mpp::shape::shapeDerivative<degree, T>;
  static constexpr int numNodalPointsFromDegree = ((degree + 1) * (degree + 2) / 2);
public:
  TriangularShape() : LagrangeShapeT<T, sDim, tDim>(numNodalPointsFromDegree) {}

  void NodalPoints(const Cell &c, vector<PointT<T, sDim, tDim>> &z) const override {
    z.resize(this->numNodalPoints);
    LagrangeNodalPointsTriangle(c, z, degree);
  }

  T operator()(const PointT<T, sDim, tDim> &z, int i) const override final {
    T psi[3] = {1 - z[0] - z[1], z[0], z[1]};
    return calculateTriangle<functions, numNodalPointsFromDegree, T>(psi, i, degree);
  }

  VectorFieldT<T, sDim> LocalGradient(const PointT<T, sDim, tDim> &z, int i) const override final {
    T psi[3] = {1 - z[0] - z[1], z[0], z[1]};
    return calculateDerivative<functions, functionDerivatives, numNodalPointsFromDegree, sDim,
                               T>(psi, i, degree);
  }

  const std::string Name() const override { return std::string{'P', '0' + degree, 'T', 'r', 'i'}; }
};

#endif // LAGRANGESHAPESTRIANGULAR_HPP
