#ifndef RTSHAPES_HPP
#define RTSHAPES_HPP

#include "RTShapesHexahedron.hpp"
#include "RTShapesInterval.hpp"
#include "RTShapesQuadrilateral.hpp"
#include "RTShapesTetrahedron.hpp"
#include "RTShapesTriangular.hpp"

template<typename T, int sDim, int tDim>
ShapeT<T, sDim, tDim> *createRTShape(const CELLTYPE cellType, int degree) {
  switch (cellType) {
  case INTERVAL:
    switch (degree) {
    case 0:
      return new RT0IntT<T, sDim, tDim>();
    default:
      THROW("Shape not implemented")
    }
  case TRIANGLE:
    switch (degree) {
    case 0:
      return new RTTriangle<0, T, sDim, tDim>();
    case 1:
      return new RTTriangle<1, T, sDim, tDim>();
    case 2:
      return new RTTriangle<2, T, sDim, tDim>();
    case 3:
      return new RTTriangle<3, T, sDim, tDim>();
    case 4:
      return new RTTriangle<4, T, sDim, tDim>();
    case 5:
      return new RTTriangle<5, T, sDim, tDim>();
    case 6:
      return new RTTriangle<6, T, sDim, tDim>();
    case 7:
      return new RTTriangle<7, T, sDim, tDim>();
    case 8:
      return new RTTriangle<8, T, sDim, tDim>();
    default:
      THROW("Shape not implemented")
    }
  case QUADRILATERAL:
    switch (degree) {
    case 0:
      return new RT0QuadT<T, sDim, tDim>();
    default:
      THROW("Shape not implemented")
    }
  case TETRAHEDRON:
    switch (degree) {
    case 0:
      return new RT0TetT<T, sDim, tDim>();
    default:
      THROW("Shape not implemented")
    }
  case HEXAHEDRON:
    switch (degree) {
    case 0:
      return new RT0HexT<T, sDim, tDim>();
    default:
      THROW("Shape not implemented")
    }
  default:
    THROW("Shape not implemented")
  }
}

#endif // RTSHAPES_HPP
