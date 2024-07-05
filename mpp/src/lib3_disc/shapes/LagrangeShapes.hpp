#ifndef LAGRANGESHAPES_HPP
#define LAGRANGESHAPES_HPP

#include <memory>
#include "LagrangeShapesHexahedron.hpp"
#include "LagrangeShapesInterval.hpp"
#include "LagrangeShapesQuadrilateral.hpp"
#include "LagrangeShapesTetrahedron.hpp"
#include "LagrangeShapesTriangular.hpp"

template<typename T, int sDim, int tDim>
ShapeT<T, sDim, tDim> *createLagrangeShape(const CELLTYPE cellType, int degree) {
  switch (cellType) {
  case INTERVAL:
    switch (degree) {
    case 0:
      return new IntervalShape<0, T, sDim, tDim>();
    case 1:
      return new IntervalShape<1, T, sDim, tDim>();
    case 2:
      return new IntervalShape<2, T, sDim, tDim>();
    case 3:
      return new IntervalShape<3, T, sDim, tDim>();
    case 4:
      return new IntervalShape<4, T, sDim, tDim>();
    case 5:
      return new IntervalShape<5, T, sDim, tDim>();
    case 6:
      return new IntervalShape<6, T, sDim, tDim>();
    case 7:
      return new IntervalShape<7, T, sDim, tDim>();
    case 8:
      return new IntervalShape<8, T, sDim, tDim>();
    default:
      THROW("Shape not implemented")
    }
  case TRIANGLE:
    switch (degree) {
    case 0:
      return new TriangularShape<0, T, sDim, tDim>();
    case 1:
      return new TriangularShape<1, T, sDim, tDim>();
    case 2:
      return new TriangularShape<2, T, sDim, tDim>();
    case 3:
      return new TriangularShape<3, T, sDim, tDim>();
    case 4:
      return new TriangularShape<4, T, sDim, tDim>();
    case 5:
      return new TriangularShape<5, T, sDim, tDim>();
    case 6:
      return new TriangularShape<6, T, sDim, tDim>();
    case 7:
      WarningOnMaster("P7Tri not fully implemented") return new TriangularShape<7, T, sDim, tDim>();
    case 8:
      WarningOnMaster("P8Tri not fully implemented") return new TriangularShape<8, T, sDim, tDim>();
    default:
      THROW("Shape not implemented")
    }
  case QUADRILATERAL:
    switch (degree) {
    case 0:
      return new QuadraticShape<0, T, sDim, tDim>();
    case 1:
      return new QuadraticShape<1, T, sDim, tDim>();
    case 2:
      return new QuadraticShape<2, T, sDim, tDim>();
    case 3:
      return new QuadraticShape<3, T, sDim, tDim>();
    case 4:
      return new QuadraticShape<4, T, sDim, tDim>();
    case 5:
      return new QuadraticShape<5, T, sDim, tDim>();
    case 6:
      return new QuadraticShape<6, T, sDim, tDim>();
    case 7:
      return new QuadraticShape<7, T, sDim, tDim>();
    case 8:
      return new QuadraticShape<8, T, sDim, tDim>();
    default:
      THROW("Shape not implemented")
    }
  case TETRAHEDRON:
    switch (degree) {
    case 0:
      return new P0TetT<T, sDim, tDim>();
    case 1:
      return new P1TetT<T, sDim, tDim>();
    case 2:
      return new P2TetT<T, sDim, tDim>();
    case 3:
      return new P3TetT<T, sDim, tDim>();
    case 4:
      return new P4TetT<T, sDim, tDim>();
    default:
      THROW("Shape not implemented")
    }
  case HEXAHEDRON:
  case SPACETIME_QUADRILATERAL:
    switch (degree) {
    case 0:
      return new HexahedronShape<0, T, sDim, tDim>();
    case 1:
      return new HexahedronShape<1, T, sDim, tDim>();
    case 2:
      return new HexahedronShape<2, T, sDim, tDim>();
    case 3:
      return new HexahedronShape<3, T, sDim, tDim>();
    case 4:
      return new HexahedronShape<4, T, sDim, tDim>();
    case 5:
      return new HexahedronShape<5, T, sDim, tDim>();
    case 6:
      return new HexahedronShape<6, T, sDim, tDim>();
    case 7:
      return new HexahedronShape<7, T, sDim, tDim>();
    case 8:
      return new HexahedronShape<8, T, sDim, tDim>();
    default:
      THROW("Shape not implemented")
    }
  default:
    THROW("Shape not implemented")
  }
}

template<typename T, int sDim = SpaceDimension, int tDim = TimeDimension>
std::unique_ptr<ShapeT<T, sDim, tDim>> createLagrangeShapeUnique(const CELLTYPE cellType,
                                                                 int degree) {
  return std::unique_ptr<ShapeT<T, sDim, tDim>>(
      createLagrangeShape<T, sDim, tDim>(cellType, degree));
}

#endif // LAGRANGESHAPES_HPP
