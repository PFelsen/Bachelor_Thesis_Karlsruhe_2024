#ifndef CURLSHAPESHEXAHEDRON_HPP
#define CURLSHAPESHEXAHEDRON_HPP

#include "Shape.hpp"

template<typename T = double, int sDim = SpaceDimension, int tDim = TimeDimension>
class P1HexCurlT : public ShapeT<T, sDim, tDim> {
public:
  const std::string Name() const override { return "P1HexCurl"; }

  VectorFieldT<T, sDim> LocalVector(const Point &z, int i) const override {
    switch (i) {
    case 0:
      return VectorFieldT<T, sDim>((1 - z[1]) * (1 - z[2]), 0.0, 0.0);
    case 1:
      return VectorFieldT<T, sDim>(0.0, z[0] * (1 - z[2]), 0.0);
    case 2:
      return VectorFieldT<T, sDim>(-z[1] * (1 - z[2]), 0.0, 0.0);
    case 3:
      return VectorFieldT<T, sDim>(0.0, (z[0] - 1) * (1 - z[2]), 0.0);
    case 4:
      return VectorFieldT<T, sDim>(0.0, 0.0, (1 - z[0]) * (1 - z[1]));
    case 5:
      return VectorFieldT<T, sDim>(0.0, 0.0, z[0] * (1 - z[1]));
    case 6:
      return VectorFieldT<T, sDim>(0.0, 0.0, z[0] * z[1]);
    case 7:
      return VectorFieldT<T, sDim>(0.0, 0.0, z[1] * (1 - z[0]));
    case 8:
      return VectorFieldT<T, sDim>((1 - z[1]) * z[2], 0.0, 0.0);
    case 9:
      return VectorFieldT<T, sDim>(0.0, z[0] * z[2], 0.0);
    case 10:
      return VectorFieldT<T, sDim>(-z[1] * z[2], 0.0, 0.0);
    case 11:
      return VectorFieldT<T, sDim>(0.0, z[2] * (z[0] - 1), 0.0);
    default:
      THROW("Not implemented")
    }
  }

  VectorFieldT<T, sDim> LocalCurl(const Point &z, int i) const override {
    switch (i) {
    case 0:
      return VectorFieldT<T, sDim>(0.0, z[1] - 1, 1 - z[2]);
    case 1:
      return VectorFieldT<T, sDim>(z[0], 0.0, 1 - z[2]);
    case 2:
      return VectorFieldT<T, sDim>(0.0, z[1], 1 - z[2]);
    case 3:
      return VectorFieldT<T, sDim>(z[0] - 1, 0.0, 1 - z[2]);
    case 4:
      return VectorFieldT<T, sDim>(z[0] - 1, 1 - z[1], 0.0);
    case 5:
      return VectorFieldT<T, sDim>(-z[0], z[1] - 1, 0.0);
    case 6:
      return VectorFieldT<T, sDim>(z[0], -z[1], 0.0);
    case 7:
      return VectorFieldT<T, sDim>(1 - z[0], z[1], 0.0);
    case 8:
      return VectorFieldT<T, sDim>(0.0, 1 - z[1], z[2]);
    case 9:
      return VectorFieldT<T, sDim>(-z[0], 0.0, z[2]);
    case 10:
      return VectorFieldT<T, sDim>(0.0, -z[1], z[2]);
    case 11:
      return VectorFieldT<T, sDim>(1 - z[0], 0.0, z[2]);
    default:
      THROW("Not implementd")
    }
  }

  P1HexCurlT() : ShapeT<T, sDim, tDim>(12) {}
};

using P1HexCurl = P1HexCurlT<>;

template<typename T = double, int sDim = SpaceDimension, int tDim = TimeDimension>
class P2HexCurlT : public ShapeT<T, sDim, tDim> {
public:
  const std::string Name() const override { return "P2HexCurl"; }
public:
  VectorField LocalVector(const Point &z, int i) const {
    switch (i) {
      // the lowest order edge-based
    case 0:
      return VectorField((1 - z[1]) * (1 - z[2]), 0.0, 0.0);
    case 1:
      return VectorField(0.0, z[0] * (1 - z[2]), 0.0);
    case 2:
      return VectorField(-z[1] * (1 - z[2]), 0.0, 0.0);
    case 3:
      return VectorField(0.0, (z[0] - 1) * (1 - z[2]), 0.0);
    case 4:
      return VectorField(0.0, 0.0, (1 - z[0]) * (1 - z[1]));
    case 5:
      return VectorField(0.0, 0.0, z[0] * (1 - z[1]));
    case 6:
      return VectorField(0.0, 0.0, z[0] * z[1]);
    case 7:
      return VectorField(0.0, 0.0, z[1] * (1 - z[0]));
    case 8:
      return VectorField((1 - z[1]) * z[2], 0.0, 0.0);
    case 9:
      return VectorField(0.0, z[0] * z[2], 0.0);
    case 10:
      return VectorField(-z[1] * z[2], 0.0, 0.0);
    case 11:
      return VectorField(0.0, z[2] * (z[0] - 1), 0.0);
      // 2-nd order edge-based
    case 12:
      return VectorField(2 * (-1 + 2 * z[0]) * (-1 + z[1]) * (-1 + z[2]),
                         2 * (-1 + z[0]) * z[0] * (-1 + z[2]),
                         2 * (-1 + z[0]) * z[0] * (-1 + z[1]));
    case 13:
      return VectorField(-2 * (-1 + z[1]) * z[1] * (-1 + z[2]),
                         -2 * z[0] * (-1 + 2 * z[1]) * (-1 + z[2]), -2 * z[0] * (-1 + z[1]) * z[1]);
    case 14:
      return VectorField(-2 * (-1 + 2 * z[0]) * z[1] * (-1 + z[2]),
                         -2 * (-1 + z[0]) * z[0] * (-1 + z[2]), -2 * (-1 + z[0]) * z[0] * z[1]);
    case 15:
      return VectorField(2 * (-1 + z[1]) * z[1] * (-1 + z[2]),
                         2 * (-1 + z[0]) * (-1 + 2 * z[1]) * (-1 + z[2]),
                         2 * (-1 + z[0]) * (-1 + z[1]) * z[1]);
    case 16:
      return VectorField(2 * (-1 + z[1]) * (-1 + z[2]) * z[2], 2 * (-1 + z[0]) * (-1 + z[2]) * z[2],
                         2 * (-1 + z[0]) * (-1 + z[1]) * (-1 + 2 * z[2]));
    case 17:
      return VectorField(-2 * (-1 + z[1]) * (-1 + z[2]) * z[2], -2 * z[0] * (-1 + z[2]) * z[2],
                         -2 * z[0] * (-1 + z[1]) * (-1 + 2 * z[2]));
    case 18:
      return VectorField(2 * z[1] * (-1 + z[2]) * z[2], 2 * z[0] * (-1 + z[2]) * z[2],
                         2 * z[0] * z[1] * (-1 + 2 * z[2]));
    case 19:
      return VectorField(-2 * z[1] * (-1 + z[2]) * z[2], -2 * (-1 + z[0]) * (-1 + z[2]) * z[2],
                         -2 * (-1 + z[0]) * z[1] * (-1 + 2 * z[2]));
    case 20:
      return VectorField(-2 * (-1 + 2 * z[0]) * (-1 + z[1]) * z[2], -2 * (-1 + z[0]) * z[0] * z[2],
                         -2 * (-1 + z[0]) * z[0] * (-1 + z[1]));
    case 21:
      return VectorField(2 * (-1 + z[1]) * z[1] * z[2], 2 * z[0] * (-1 + 2 * z[1]) * z[2],
                         2 * z[0] * (-1 + z[1]) * z[1]);
    case 22:
      return VectorField(2 * (-1 + 2 * z[0]) * z[1] * z[2], 2 * (-1 + z[0]) * z[0] * z[2],
                         2 * (-1 + z[0]) * z[0] * z[1]);
    case 23:
      return VectorField(-2 * (-1 + z[1]) * z[1] * z[2], -2 * (-1 + z[0]) * (-1 + 2 * z[1]) * z[2],
                         -2 * (-1 + z[0]) * (-1 + z[1]) * z[1]);
      // 2-nd order face based
      // face (1, 0, 3, 2)
    case 24:
      return VectorField(4 * (z[1] - 1) * z[1] * (z[2] - 1), 0, 0);
    case 25:
      return VectorField(0, -4 * (z[0] - 1) * z[0] * (z[2] - 1), 0);
    case 26:
      return VectorField(-4 * (2 * z[0] - 1) * (z[1] - 1) * z[1] * (z[2] - 1),
                         4 * (z[0] - 1) * z[0] * (2 * z[1] - 1) * (z[2] - 1), 0);
    case 27:
      return VectorField(-4 * (-1 + 2 * z[0]) * (-1 + z[1]) * z[1] * (-1 + z[2]),
                         -4 * (-1 + z[0]) * z[0] * (-1 + 2 * z[1]) * (-1 + z[2]),
                         -4 * (-1 + z[0]) * z[0] * (-1 + z[1]) * z[1]);
      // face (0, 1, 5, 4)
    case 28:
      return VectorField(-4 * (z[1] - 1) * (z[2] - 1) * z[2], 0, 0);
    case 29:
      return VectorField(0, 0, -4 * (z[0] - 1) * z[0] * (z[1] - 1));
    case 30:
      return VectorField(-4 * (2 * z[0] - 1) * (z[1] - 1) * (z[2] - 1) * z[2], 0,
                         4 * (z[0] - 1) * z[0] * (z[1] - 1) * (2 * z[2] - 1));
    case 31:
      return VectorField(-4 * (-1 + 2 * z[0]) * (-1 + z[1]) * (-1 + z[2]) * z[2],
                         -4 * (-1 + z[0]) * z[0] * (-1 + z[2]) * z[2],
                         -4 * (-1 + z[0]) * z[0] * (-1 + z[1]) * (-1 + 2 * z[2]));
      // face (1, 2, 6, 5)
    case 32:
      return VectorField(0, 4 * z[0] * (z[2] - 1) * z[2], 0);
    case 33:
      return VectorField(0, 0, 4 * z[0] * (z[1] - 1) * z[1]);
    case 34:
      return VectorField(0, 4 * z[0] * (2 * z[1] - 1) * (z[2] - 1) * z[2],
                         -4 * z[0] * (z[1] - 1) * z[1] * (2 * z[2] - 1));
    case 35:
      return VectorField(4 * (-1 + z[1]) * z[1] * (-1 + z[2]) * z[2],
                         4 * z[0] * (-1 + 2 * z[1]) * (-1 + z[2]) * z[2],
                         4 * z[0] * (-1 + z[1]) * z[1] * (-1 + 2 * z[2]));
      // face (2, 3, 7, 6)
    case 36:
      return VectorField(-4 * z[1] * (z[2] - 1) * z[2], 0, 0);
    case 37:
      return VectorField(0, 0, 4 * (z[0] - 1) * z[0] * z[1]);
    case 38:
      return VectorField(4 * (2 * z[0] - 1) * z[1] * (z[2] - 1) * z[2], 0,
                         -4 * (z[0] - 1) * z[0] * z[1] * (2 * z[2] - 1));
    case 39:
      return VectorField(4 * (-1 + 2 * z[0]) * z[1] * (-1 + z[2]) * z[2],
                         4 * (-1 + z[0]) * z[0] * (-1 + z[2]) * z[2],
                         4 * (-1 + z[0]) * z[0] * z[1] * (-1 + 2 * z[2]));
      // face (3, 0, 4, 7)
    case 40:
      return VectorField(0, 4 * (z[0] - 1) * (z[2] - 1) * z[2], 0);
    case 41:
      return VectorField(0, 0, -4 * (z[0] - 1) * (z[1] - 1) * z[1]);
    case 42:
      return VectorField(0, -4 * (z[0] - 1) * (2 * z[1] - 1) * (z[2] - 1) * z[2],
                         4 * (z[0] - 1) * (z[1] - 1) * z[1] * (2 * z[2] - 1));
    case 43:
      return VectorField(-4 * (-1 + z[1]) * z[1] * (-1 + z[2]) * z[2],
                         -4 * (-1 + z[0]) * (-1 + 2 * z[1]) * (-1 + z[2]) * z[2],
                         -4 * (-1 + z[0]) * (-1 + z[1]) * z[1] * (-1 + 2 * z[2]));
      // face (4, 5, 6, 7)
    case 44:
      return VectorField(4 * (z[1] - 1) * z[1] * z[2], 0, 0);
    case 45:
      return VectorField(0, 4 * (z[0] - 1) * z[0] * z[2], 0);
    case 46:
      return VectorField(4 * (2 * z[0] - 1) * (z[1] - 1) * z[1] * z[2],
                         -4 * (z[0] - 1) * z[0] * (2 * z[1] - 1) * z[2], 0);
    case 47:
      return VectorField(4 * (-1 + 2 * z[0]) * (-1 + z[1]) * z[1] * z[2],
                         4 * (-1 + z[0]) * z[0] * (-1 + 2 * z[1]) * z[2],
                         4 * (-1 + z[0]) * z[0] * (-1 + z[1]) * z[1]);
      // 2-nd order cell-based
    case 48:
      return VectorField(4 * z[1] * (z[1] - 1) * z[2] * (z[2] - 1), 0, 0);
    case 49:
      return VectorField(0, 4 * z[0] * (z[0] - 1) * z[2] * (z[2] - 1), 0);
    case 50:
      return VectorField(0, 0, 4 * z[0] * (z[0] - 1) * z[1] * (z[1] - 1));
    case 51:
      return VectorField(8 * (2 * z[0] - 1) * z[1] * (z[1] - 1) * z[2] * (z[2] - 1),
                         -8 * (2 * z[1] - 1) * z[0] * (z[0] - 1) * z[2] * (z[2] - 1),
                         8 * (2 * z[2] - 1) * z[0] * (z[0] - 1) * z[1] * (z[1] - 1));
    case 52:
      return VectorField(8 * (2 * z[0] - 1) * z[1] * (z[1] - 1) * z[2] * (z[2] - 1),
                         8 * (2 * z[1] - 1) * z[0] * (z[0] - 1) * z[2] * (z[2] - 1),
                         -8 * (2 * z[2] - 1) * z[0] * (z[0] - 1) * z[1] * (z[1] - 1));
    case 53:
      return VectorField(8 * (-1 + 2 * z[0]) * (-1 + z[1]) * z[1] * (-1 + z[2]) * z[2],
                         8 * (-1 + z[0]) * z[0] * (-1 + 2 * z[1]) * (-1 + z[2]) * z[2],
                         8 * (-1 + z[0]) * z[0] * (-1 + z[1]) * z[1] * (-1 + 2 * z[2]));
    default:
      THROW("Not implemented")
    }
  }

  VectorFieldT<T, sDim> LocalCurl(const Point &z, int i) const {
    switch (i) {
      // the lowest order edge-based
    case 0:
      return VectorFieldT<T, sDim>(0.0, z[1] - 1, 1 - z[2]);
    case 1:
      return VectorFieldT<T, sDim>(z[0], 0.0, 1 - z[2]);
    case 2:
      return VectorFieldT<T, sDim>(0.0, z[1], 1 - z[2]);
    case 3:
      return VectorFieldT<T, sDim>(z[0] - 1, 0.0, 1 - z[2]);
    case 4:
      return VectorFieldT<T, sDim>(z[0] - 1, 1 - z[1], 0.0);
    case 5:
      return VectorFieldT<T, sDim>(-z[0], z[1] - 1, 0.0);
    case 6:
      return VectorFieldT<T, sDim>(z[0], -z[1], 0.0);
    case 7:
      return VectorFieldT<T, sDim>(1 - z[0], z[1], 0.0);
    case 8:
      return VectorFieldT<T, sDim>(0.0, 1 - z[1], z[2]);
    case 9:
      return VectorFieldT<T, sDim>(-z[0], 0.0, z[2]);
    case 10:
      return VectorFieldT<T, sDim>(0.0, -z[1], z[2]);
    case 11:
      return VectorFieldT<T, sDim>(1 - z[0], 0.0, z[2]);
      // 2-nd order edge-based
    case 12:
      return VectorFieldT<T, sDim>(0, 0, 0);
    case 13:
      return VectorFieldT<T, sDim>(0, 0, 0);
    case 14:
      return VectorFieldT<T, sDim>(0, 0, 0);
    case 15:
      return VectorFieldT<T, sDim>(0, 0, 0);
    case 16:
      return VectorFieldT<T, sDim>(0, 0, 0);
    case 17:
      return VectorFieldT<T, sDim>(0, 0, 0);
    case 18:
      return VectorFieldT<T, sDim>(0, 0, 0);
    case 19:
      return VectorFieldT<T, sDim>(0, 0, 0);
    case 20:
      return VectorFieldT<T, sDim>(0, 0, 0);
    case 21:
      return VectorFieldT<T, sDim>(0, 0, 0);
    case 22:
      return VectorFieldT<T, sDim>(0, 0, 0);
    case 23:
      return VectorFieldT<T, sDim>(0, 0, 0);
      // 2-nd order face based
      // face (1, 0, 3, 2)
    case 24:
      return VectorFieldT<T, sDim>(0, 4 * (z[1] - 1) * z[1], -4 * (2 * z[1] - 1) * (z[2] - 1));
    case 25:
      return VectorFieldT<T, sDim>(4 * (z[0] - 1) * z[0], 0, -4 * (2 * z[0] - 1) * (z[2] - 1));
    case 26:
      return VectorFieldT<T, sDim>(-4 * (z[0] - 1) * z[0] * (2 * z[1] - 1),
                                   -4 * (2 * z[0] - 1) * (z[1] - 1) * z[1],
                                   8 * (2 * z[0] - 1) * (2 * z[1] - 1) * (z[2] - 1));
    case 27:
      return VectorFieldT<T, sDim>(0, 0, 0);
      // face (0, 1, 5, 4)
    case 28:
      return VectorFieldT<T, sDim>(0, -4 * (z[1] - 1) * (2 * z[2] - 1), 4 * (z[2] - 1) * z[2]);
    case 29:
      return VectorFieldT<T, sDim>(-4 * (z[0] - 1) * z[0], 4 * (2 * z[0] - 1) * (z[1] - 1), 0);
    case 30:
      return VectorFieldT<T, sDim>(4 * (z[0] - 1) * z[0] * (2 * z[2] - 1),
                                   -8 * (2 * z[0] - 1) * (z[1] - 1) * (2 * z[2] - 1),
                                   4 * (2 * z[0] - 1) * (z[2] - 1) * z[2]);
    case 31:
      return VectorFieldT<T, sDim>(0, 0, 0);
      // face (1, 2, 6, 5)
    case 32:
      return VectorFieldT<T, sDim>(-4 * z[0] * (2 * z[2] - 1), 0, 4 * (z[2] - 1) * z[2]);
    case 33:
      return VectorFieldT<T, sDim>(4 * z[0] * (2 * z[1] - 1), -4 * (z[1] - 1) * z[1], 0);
    case 34:
      return VectorFieldT<T, sDim>(-8 * z[0] * (2 * z[1] - 1) * (2 * z[2] - 1),
                                   4 * (z[1] - 1) * z[1] * (2 * z[2] - 1),
                                   4 * (2 * z[1] - 1) * (z[2] - 1) * z[2]);
    case 35:
      return VectorFieldT<T, sDim>(0, 0, 0);
      // face (2, 3, 7, 6)
    case 36:
      return VectorFieldT<T, sDim>(0, -4 * z[1] * (2 * z[2] - 1), 4 * (z[2] - 1) * z[2]);
    case 37:
      return VectorFieldT<T, sDim>(4 * (z[0] - 1) * z[0], -4 * (2 * z[0] - 1) * z[1], 0);
    case 38:
      return VectorFieldT<T, sDim>(-4 * (z[0] - 1) * z[0] * (2 * z[2] - 1),
                                   8 * (2 * z[0] - 1) * z[1] * (2 * z[2] - 1),
                                   -4 * (2 * z[0] - 1) * (z[2] - 1) * z[2]);
    case 39:
      return VectorFieldT<T, sDim>(0, 0, 0);
      // face (3, 0, 4, 7)
    case 40:
      return VectorFieldT<T, sDim>(-4 * (z[0] - 1) * (2 * z[2] - 1), 0, 4 * (z[2] - 1) * z[2]);
    case 41:
      return VectorFieldT<T, sDim>(-4 * (z[0] - 1) * (2 * z[1] - 1), 4 * (z[1] - 1) * z[1], 0);
    case 42:
      return VectorFieldT<T, sDim>(8 * (z[0] - 1) * (2 * z[1] - 1) * (2 * z[2] - 1),
                                   -4 * (z[1] - 1) * z[1] * (2 * z[2] - 1),
                                   -4 * (2 * z[1] - 1) * (z[2] - 1) * z[2]);
    case 43:
      return VectorFieldT<T, sDim>(0, 0, 0);
      // face (4, 5, 6, 7)
    case 44:
      return VectorFieldT<T, sDim>(0, 4 * (z[1] - 1) * z[1], -4 * (2 * z[1] - 1) * z[2]);
    case 45:
      return VectorFieldT<T, sDim>(-4 * (z[0] - 1) * z[0], 0, 4 * (2 * z[0] - 1) * z[2]);
    case 46:
      return VectorFieldT<T, sDim>(4 * (z[0] - 1) * z[0] * (2 * z[1] - 1),
                                   4 * (2 * z[0] - 1) * (z[1] - 1) * z[1],
                                   -8 * (2 * z[0] - 1) * (2 * z[1] - 1) * z[2]);
    case 47:
      return VectorFieldT<T, sDim>(0, 0, 0);
      // 2-nd order cell-based
    case 48:
      return VectorFieldT<T, sDim>(0, 4 * (z[1] - 1) * z[1] * (2 * z[2] - 1),
                                   -4 * (2 * z[1] - 1) * (z[2] - 1) * z[2]);
    case 49:
      return VectorFieldT<T, sDim>(-4 * (z[0] - 1) * z[0] * (2 * z[2] - 1), 0,
                                   4 * (2 * z[0] - 1) * (z[2] - 1) * z[2]);
    case 50:
      return VectorFieldT<T, sDim>(4 * (z[0] - 1) * z[0] * (2 * z[1] - 1),
                                   -4 * (2 * z[0] - 1) * (z[1] - 1) * z[1], 0);
    case 51:
      return VectorFieldT<T, sDim>(16 * (z[0] - 1) * z[0] * (2 * z[1] - 1) * (2 * z[2] - 1), 0,
                                   -16 * (2 * z[0] - 1) * (2 * z[1] - 1) * (z[2] - 1) * z[2]);
    case 52:
      return VectorFieldT<T, sDim>(-16 * (z[0] - 1) * z[0] * (2 * z[1] - 1) * (2 * z[2] - 1),
                                   16 * (2 * z[0] - 1) * (z[1] - 1) * z[1] * (2 * z[2] - 1), 0);
    case 53:
      return VectorFieldT<T, sDim>(0, 0, 0);
    default:
      THROW("Not implemented")
    }
  }

  P2HexCurlT() : ShapeT<T, sDim, tDim>(54) {}
};

using P2HexCurl = P2HexCurlT<>;

#endif // CURLSHAPESHEXAHEDRON_HPP
