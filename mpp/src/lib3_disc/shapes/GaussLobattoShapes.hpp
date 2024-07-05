#ifndef GAUSSLOBATTOSHAPES_HPP
#define GAUSSLOBATTOSHAPES_HPP

#include "GaussLobattoShape.hpp"

template<typename T, int sDim, int tDim>
ShapeT<T, sDim, tDim> *createGaussLobattoShape(const CELLTYPE cellType, int degree) {
  switch (cellType) {
  case INTERVAL:
    return new GaussLobattoShape<T, sDim, tDim>(degree, 1);
  case QUADRILATERAL:
    return new GaussLobattoShape<T, sDim, tDim>(degree, 2);
  default:
    THROW("Shape not implemented")
  }
}

#endif // GAUSSLOBATTOSHAPES_HPP
