#ifndef CURLSHAPES_HPP
#define CURLSHAPES_HPP

#include "CurlShapesHexahedron.hpp"
#include "CurlShapesTetrahedron.hpp"

template<typename T, int sDim, int tDim>
ShapeT<T, sDim, tDim> *createCurlShape(const CELLTYPE cellType, int degree) {
  switch (cellType) {
  case TETRAHEDRON:
    switch (degree) {
    case 1:
      return new P1TetCurlT<T, sDim, tDim>();
    case 2:
      return new P2TetCurlT<T, sDim, tDim>();
    default:
      THROW("Shape not implemented")
    }
  case HEXAHEDRON:
    switch (degree) {
    case 1:
      return new P1HexCurlT<T, sDim, tDim>();
    case 2:
      return new P2HexCurlT<T, sDim, tDim>();
    default:
      THROW("Shape not implemented")
    }
  default:
    THROW("Shape not implemented")
  }
}

#endif // CURLSHAPES_HPP
