#ifndef SIMPLEXPOLYNOMIALS_HPP
#define SIMPLEXPOLYNOMIALS_HPP

#include "Assertion.hpp"

namespace poly {
template<typename T = double, int degree>
T shapeFunction(T x, int i) {

  switch (i) {
  case 0:
    return T(1.0);
  case 1:
    return (degree * x);
  case 2:
    return (degree * x) * (degree * x - 1) / T(2.0);
  case 3:
    return (degree * x) * (degree * x - 1) * (degree * x - 2) / T(6.0);
  case 4:
    return (degree * x) * (degree * x - 1) * (degree * x - 2) * (degree * x - 3) / T(24.0);
  case 5:
    return (degree * x) * (degree * x - 1) * (degree * x - 2) * (degree * x - 3) * (degree * x - 4)
           / T(120.0);
  case 6:
    return (degree * x) * (degree * x - 1) * (degree * x - 2) * (degree * x - 3) * (degree * x - 4)
           * (degree * x - 5) / T(720.0);
  default:
    THROW("Polynomial of degree '" + std::to_string(i) + "' not implemented.")
  }
}

template<typename T = double, int degree>
T shapeDerivative(T x, int i) {
  int di[degree]{};
  di[0] = degree;
  T xi[degree]{};
  xi[0] = x;

  for (int power = 1; power < i; ++power) {
    xi[power] = xi[power - 1] * x;
    di[power] = di[power - 1] * degree;
  }

  switch (i) {
  case 0:
    return T(0.0);
  case 1:
    return T(degree);
  case 2:
    return di[1] / T(1.0) * xi[0] - di[0] / T(2.0);
  case 3:
    return di[2] / T(2.0) * xi[1] + (-6 * di[1] * xi[0] + 2 * di[0]) / T(6.0);
  case 4:
    return di[3] / T(6.0) * xi[2]
           + (-18 * di[2] * xi[1] + 22 * di[1] * xi[0] - 6 * di[0]) / T(24.0);
  case 5:
    return di[4] / T(24.0) * xi[3]
           + (-40 * di[3] * xi[2] + 105 * di[2] * xi[1] - 100 * di[1] * xi[0] + 24 * di[0])
                 / T(120.0);
  case 6:
    return di[5] / T(120.0) * xi[4]
           + (-75 * di[4] * xi[3] + 340 * di[3] * xi[2] - 675 * di[2] * xi[1] + 548 * di[1] * xi[0]
              - 120 * di[0])
                 / T(720.0);
  default:
    Exit("Derivatives of order higher than '" + std::to_string(6) + "'are not implemented.")
  }
}
} // namespace poly
#endif // SIMPLEXPOLYNOMIALS_HPP
