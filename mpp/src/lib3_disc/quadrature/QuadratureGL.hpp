#ifndef QUADRATUREINTERVAL_HPP
#define QUADRATUREINTERVAL_HPP

#include "QuadratureT.hpp"

template<typename T = double, int sDim = SpaceDimension, int tDim = TimeDimension>
class QintGaussLobattoT : public QuadratureT<T, sDim, tDim> {
  /**
   * GaussLobatto-Quadrature is an collocation method.
   * The interval boundaries 0 and 1 are guaranteed to be quadrature points.
   * Intermediate points are chosen to reach highest possible order.
   * For n points GL-Quadrature is exact up to polynomial degree 2 * n - 3.
   **/
public:
  QintGaussLobattoT(int npcount);
};

using QintGaussLobatto = QintGaussLobattoT<>;

#endif // QUADRATUREINTERVAL_HPP
