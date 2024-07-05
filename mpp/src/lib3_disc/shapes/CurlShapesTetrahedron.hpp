#ifndef CURLSHAPESTETRAHEDRON_HPP
#define CURLSHAPESTETRAHEDRON_HPP

#include "Shape.hpp"

// TODO: Maybe use P1Tet???
template<typename T = double, int sDim = SpaceDimension, int tDim = TimeDimension>
class P1TetCurlT : public ShapeT<T, sDim, tDim> {
  const std::string Name() const override { return "P1TetCurl"; }

  T shape(const PointT<T, sDim, tDim> &z, int i) const {
    switch (i) {
    case 0:
      return 1 - z[0] - z[1] - z[2];
    case 1:
      return z[0];
    case 2:
      return z[1];
    case 3:
      return z[2];
    default:
      THROW("Not implemented")
    }
  }

  VectorFieldT<T, sDim> gradient(int i) const {
    switch (i) {
    case 0:
      return VectorFieldT<T, sDim>(T(-1.0), T(-1.0), T(-1.0));
    case 1:
      return VectorFieldT<T, sDim>(T(1.0), T(0.0), T(0.0));
    case 2:
      return VectorFieldT<T, sDim>(T(0.0), T(1.0), T(0.0));
    case 3:
      return VectorFieldT<T, sDim>(T(0.0), T(0.0), T(1.0));
    default:
      THROW("Not implemented")
    }
  }
public:
  VectorFieldT<T, sDim> LocalVector(const PointT<T, sDim, tDim> &z, int i) const override {
    switch (i) {
    case 0:
      return shape(z, 0) * gradient(1) - shape(z, 1) * gradient(0);
    case 1:
      return shape(z, 1) * gradient(2) - shape(z, 2) * gradient(1);
    case 2:
      return shape(z, 0) * gradient(2) - shape(z, 2) * gradient(0);
    case 3:
      return shape(z, 0) * gradient(3) - shape(z, 3) * gradient(0);
    case 4:
      return shape(z, 1) * gradient(3) - shape(z, 3) * gradient(1);
    case 5:
      return shape(z, 2) * gradient(3) - shape(z, 3) * gradient(2);
    default:
      THROW("Not implemented")
    }
  }

  VectorFieldT<T, sDim> LocalCurl(const PointT<T, sDim, tDim> &z, int i) const override {
    switch (i) {
    case 0:
      return 2 * curl(gradient(0), gradient(1));
    case 1:
      return 2 * curl(gradient(1), gradient(2));
    case 2:
      return 2 * curl(gradient(0), gradient(2));
    case 3:
      return 2 * curl(gradient(0), gradient(3));
    case 4:
      return 2 * curl(gradient(1), gradient(3));
    case 5:
      return 2 * curl(gradient(2), gradient(3));
    default:
      THROW("Not implemented")
    }
  }

  explicit P1TetCurlT() : ShapeT<T, sDim, tDim>(6) {}
};

using P1TetCurl = P1TetCurlT<>;

template<typename T = double, int sDim = SpaceDimension, int tDim = TimeDimension>
class P2TetCurlT : public ShapeT<T, sDim, tDim> {
public:
  const std::string Name() const override { return "P2TetCurl"; }

  explicit P2TetCurlT() : ShapeT<T, sDim, tDim>(20) {}

  VectorFieldT<T, sDim> LocalCurl(const PointT<T, sDim, tDim> &z, int i) const override {
    switch (i) {
    case 0:
      return VectorFieldT<T, sDim>(0.0,
                                   6.928203230275506 * z[0] + 9.464101615137755 * z[1]
                                       + 9.464101615137753 * z[2] - 6.464101615137753,
                                   -6.92820323027553 * z[0] - 9.464101615137768 * z[1]
                                       - 9.464101615137766 * z[2] + 6.464101615137766);
    case 1:
      return VectorFieldT<T, sDim>(0.0,
                                   -6.928203230275506 * z[0] + 2.535898384862245 * z[1]
                                       + 2.535898384862247 * z[2] + 0.464101615137754,
                                   6.928203230275506 * z[0] - 2.535898384862247 * z[1]
                                       - 2.535898384862249 * z[2] - 0.464101615137752);
    case 2:
      return VectorFieldT<T, sDim>(0.0, 0.0,
                                   9.464101615137755 * z[0] + 2.535898384862247 * z[1] - 3.0);
    case 3:
      return VectorFieldT<T, sDim>(0.0, 0.0,
                                   2.535898384862253 * z[0] + 9.464101615137762 * z[1] - 3.0);
    case 4:
      return VectorFieldT<T, sDim>(-9.464101615137762 * z[0] - 6.928203230275509 * z[1]
                                       - 9.464101615137759 * z[2] + 6.464101615137755,
                                   0.0,
                                   9.464101615137762 * z[0] + 6.928203230275516 * z[1]
                                       + 9.464101615137762 * z[2] - 6.46410161513776);
    case 5:
      return VectorFieldT<T, sDim>(-2.535898384862241 * z[0] + 6.928203230275506 * z[1]
                                       - 2.535898384862247 * z[2] - 0.464101615137753,
                                   0.0,
                                   2.535898384862241 * z[0] - 6.928203230275516 * z[1]
                                       + 2.535898384862241 * z[2] + 0.464101615137758);
    case 6:
      return VectorFieldT<T, sDim>(9.464101615137753 * z[0] + 9.464101615137746 * z[1]
                                       + 6.928203230275506 * z[2] - 6.464101615137751,
                                   -9.46410161513775 * z[0] - 9.464101615137753 * z[1]
                                       - 6.928203230275506 * z[2] + 6.464101615137752,
                                   0.0);
    case 7:
      return VectorFieldT<T, sDim>(2.535898384862247 * z[0] + 2.535898384862244 * z[1]
                                       - 6.928203230275509 * z[2] + 0.464101615137756,
                                   -2.535898384862244 * z[0] - 2.535898384862247 * z[1]
                                       + 6.928203230275509 * z[2] - 0.464101615137754,
                                   0.0);
    case 8:
      return VectorFieldT<T, sDim>(0.0, -9.464101615137753 * z[0] - 2.535898384862244 * z[2] + 3.0,
                                   0.0);
    case 9:
      return VectorFieldT<T, sDim>(0.0, -2.535898384862244 * z[0] - 9.464101615137753 * z[2] + 3.0,
                                   0.0);
    case 10:
      return VectorFieldT<T, sDim>(9.464101615137753 * z[1] + 2.535898384862244 * z[2] - 3.0, 0.0,
                                   0.0);
    case 11:
      return VectorFieldT<T, sDim>(2.535898384862244 * z[1] + 9.464101615137753 * z[2] - 3.0, 0.0,
                                   0.0);
    case 12:
      return VectorFieldT<T, sDim>(0.0, -12.0 * z[1],
                                   12.0 * z[0] + 24.0 * z[1] + 12.0 * z[2] - 12.0);
    case 13:
      return VectorFieldT<T, sDim>(12.0 * z[0], 0.0,
                                   -24.0 * z[0] - 12.0 * z[1] - 12.0 * z[2] + 12.0);
    case 14:
      return VectorFieldT<T, sDim>(-12.0 * z[0], 0.0, 12.0 * z[2]);
    case 15:
      return VectorFieldT<T, sDim>(12.0 * z[0], -12.0 * z[1], 0.0);
    case 16:
      return VectorFieldT<T, sDim>(12.0 * z[0] + 12.0 * z[1] + 24.0 * z[2] - 12.0, 0.0,
                                   -12.0 * z[2]);
    case 17:
      return VectorFieldT<T, sDim>(-12.0 * z[0] - 24.0 * z[1] - 12.0 * z[2] + 12.0, 12.0 * z[1],
                                   0.0);
    case 18:
      return VectorFieldT<T, sDim>(0.0, -12.0 * z[0] - 12.0 * z[1] - 24.0 * z[2] + 12.0,
                                   12.0 * z[2]);
    case 19:
      return VectorFieldT<T, sDim>(-12.0 * z[0], 24.0 * z[0] + 12.0 * z[1] + 12.0 * z[2] - 12.0,
                                   0.0);
    default:
      THROW("Not implemented")
    }
  }

  VectorFieldT<T, sDim> LocalVector(const PointT<T, sDim, tDim> &z, int i) const override {
    switch (i) {
    case 0:
      return VectorFieldT<T, sDim>(2.309401076758505 * z[0] * z[1] + 2.309401076758502 * z[0] * z[2]
                                       - 1.732050807568878 * z[0] + 3.154700538379254 * z[1] * z[1]
                                       + 6.309401076758502 * z[1] * z[2] - 4.520725942163693 * z[1]
                                       + 3.154700538379251 * z[2] * z[2] - 4.52072594216369 * z[2]
                                       + 1.366025403784439,
                                   -2.309401076758505 * z[0] * z[0]
                                       - 3.154700538379254 * z[0] * z[1]
                                       - 3.154700538379251 * z[0] * z[2] + 1.943375672974065 * z[0],
                                   -2.309401076758502 * z[0] * z[0]
                                       - 3.154700538379251 * z[0] * z[1]
                                       - 3.154700538379251 * z[0] * z[2]
                                       + 1.943375672974064 * z[0]);
    case 1:
      return VectorFieldT<T, sDim>(-2.309401076758502 * z[0] * z[1]
                                       - 2.309401076758502 * z[0] * z[2] + 1.732050807568877 * z[0]
                                       + 0.845299461620749 * z[1] * z[1]
                                       + 1.690598923241497 * z[1] * z[2] - 0.47927405783631 * z[1]
                                       + 0.845299461620749 * z[2] * z[2] - 0.47927405783631 * z[2]
                                       - 0.366025403784439,
                                   2.309401076758502 * z[0] * z[0] - 0.845299461620749 * z[0] * z[1]
                                       - 0.845299461620749 * z[0] * z[2] - 0.943375672974063 * z[0],
                                   2.309401076758502 * z[0] * z[0] - 0.845299461620748 * z[0] * z[1]
                                       - 0.845299461620749 * z[0] * z[2]
                                       - 0.943375672974064 * z[0]);
    case 2:
      return VectorFieldT<T, sDim>(0.84529946162075 * z[0] * z[1] + 3.154700538379253 * z[1] * z[1]
                                       - 1.788675134594814 * z[1],
                                   -0.84529946162075 * z[0] * z[0] - 3.154700538379253 * z[0] * z[1]
                                       + 1.211324865405189 * z[0],
                                   0.0);
    case 3:
      return VectorFieldT<T, sDim>(3.154700538379251 * z[0] * z[1] + 0.845299461620749 * z[1] * z[1]
                                       - 1.211324865405188 * z[1],
                                   -3.154700538379251 * z[0] * z[0]
                                       - 0.845299461620749 * z[0] * z[1] + 1.788675134594813 * z[0],
                                   0.0);
    case 4:
      return VectorFieldT<T, sDim>(-3.154700538379253 * z[0] * z[1]
                                       - 2.309401076758505 * z[1] * z[1]
                                       - 3.154700538379251 * z[1] * z[2] + 1.943375672974066 * z[1],
                                   3.154700538379253 * z[0] * z[0] + 2.309401076758505 * z[0] * z[1]
                                       + 6.309401076758501 * z[0] * z[2] - 4.520725942163692 * z[0]
                                       + 2.309401076758502 * z[1] * z[2] - 1.732050807568877 * z[1]
                                       + 3.154700538379251 * z[2] * z[2] - 4.52072594216369 * z[2]
                                       + 1.366025403784439,
                                   -3.15470053837925 * z[0] * z[1] - 2.309401076758502 * z[1] * z[1]
                                       - 3.154700538379251 * z[1] * z[2]
                                       + 1.943375672974064 * z[1]);
    case 5:
      return VectorFieldT<T, sDim>(-0.845299461620747 * z[0] * z[1]
                                       + 2.309401076758504 * z[1] * z[1]
                                       - 0.845299461620749 * z[1] * z[2] - 0.943375672974065 * z[1],
                                   0.845299461620747 * z[0] * z[0] - 2.309401076758504 * z[0] * z[1]
                                       + 1.690598923241497 * z[0] * z[2] - 0.479274057836308 * z[0]
                                       - 2.309401076758502 * z[1] * z[2] + 1.732050807568877 * z[1]
                                       + 0.845299461620749 * z[2] * z[2] - 0.47927405783631 * z[2]
                                       - 0.366025403784439,
                                   -0.845299461620748 * z[0] * z[1]
                                       + 2.309401076758502 * z[1] * z[1]
                                       - 0.845299461620749 * z[1] * z[2]
                                       - 0.943375672974064 * z[1]);
    case 6:
      return VectorFieldT<T, sDim>(-3.154700538379251 * z[0] * z[2]
                                       - 3.154700538379251 * z[1] * z[2]
                                       - 2.309401076758502 * z[2] * z[2] + 1.943375672974064 * z[2],
                                   -3.154700538379251 * z[0] * z[2] - 3.15470053837925 * z[1] * z[2]
                                       - 2.309401076758502 * z[2] * z[2] + 1.943375672974063 * z[2],
                                   3.154700538379251 * z[0] * z[0] + 6.309401076758502 * z[0] * z[1]
                                       + 2.309401076758502 * z[0] * z[2] - 4.52072594216369 * z[0]
                                       + 3.15470053837925 * z[1] * z[1]
                                       + 2.309401076758502 * z[1] * z[2] - 4.520725942163688 * z[1]
                                       - 1.732050807568877 * z[2] + 1.366025403784439);
    case 7:
      return VectorFieldT<T, sDim>(-0.845299461620748 * z[0] * z[2]
                                       - 0.845299461620749 * z[1] * z[2]
                                       + 2.309401076758503 * z[2] * z[2] - 0.943375672974064 * z[2],
                                   -0.845299461620749 * z[0] * z[2]
                                       - 0.845299461620748 * z[1] * z[2]
                                       + 2.309401076758503 * z[2] * z[2] - 0.943375672974064 * z[2],
                                   0.845299461620748 * z[0] * z[0] + 1.690598923241498 * z[0] * z[1]
                                       - 2.309401076758503 * z[0] * z[2] - 0.47927405783631 * z[0]
                                       + 0.845299461620748 * z[1] * z[1]
                                       - 2.309401076758503 * z[1] * z[2] - 0.47927405783631 * z[1]
                                       + 1.732050807568877 * z[2] - 0.366025403784439);
    case 8:
      return VectorFieldT<T, sDim>(0.845299461620748 * z[0] * z[2] + 3.154700538379251 * z[2] * z[2]
                                       - 1.788675134594813 * z[2],
                                   0.0,
                                   -0.845299461620748 * z[0] * z[0]
                                       - 3.154700538379251 * z[0] * z[2]
                                       + 1.211324865405187 * z[0]);
    case 9:
      return VectorFieldT<T, sDim>(3.154700538379251 * z[0] * z[2] + 0.845299461620748 * z[2] * z[2]
                                       - 1.211324865405187 * z[2],
                                   0.0,
                                   -3.154700538379251 * z[0] * z[0]
                                       - 0.845299461620748 * z[0] * z[2]
                                       + 1.788675134594812 * z[0]);
    case 10:
      return VectorFieldT<T, sDim>(0.0,
                                   0.845299461620748 * z[1] * z[2] + 3.154700538379251 * z[2] * z[2]
                                       - 1.788675134594813 * z[2],
                                   -0.845299461620748 * z[1] * z[1]
                                       - 3.154700538379251 * z[1] * z[2]
                                       + 1.211324865405187 * z[1]);
    case 11:
      return VectorFieldT<T, sDim>(0.0,
                                   3.154700538379251 * z[1] * z[2] + 0.845299461620748 * z[2] * z[2]
                                       - 1.211324865405187 * z[2],
                                   -3.154700538379251 * z[1] * z[1]
                                       - 0.845299461620748 * z[1] * z[2]
                                       + 1.788675134594812 * z[1]);
    case 12:
      return VectorFieldT<T, sDim>(8.0 * z[0] * z[1] + 4.0 * z[1] * z[1] + 4.0 * z[1] * z[2]
                                       - 4.0 * z[1],
                                   -8.0 * z[0] * z[0] - 4.0 * z[0] * z[1] - 8.0 * z[0] * z[2]
                                       + 8.0 * z[0],
                                   4.0 * z[0] * z[1]);
    case 13:
      return VectorFieldT<T, sDim>(-4.0 * z[0] * z[1] - 8.0 * z[1] * z[1] - 8.0 * z[1] * z[2]
                                       + 8.0 * z[1],
                                   4.0 * z[0] * z[0] + 8.0 * z[0] * z[1] + 4.0 * z[0] * z[2]
                                       - 4.0 * z[0],
                                   4.0 * z[0] * z[1]);
    case 14:
      return VectorFieldT<T, sDim>(-4.0 * z[1] * z[2], 8.0 * z[0] * z[2], -4.0 * z[0] * z[1]);
    case 15:
      return VectorFieldT<T, sDim>(8.0 * z[1] * z[2], -4.0 * z[0] * z[2], -4.0 * z[0] * z[1]);
    case 16:
      return VectorFieldT<T, sDim>(4.0 * z[1] * z[2],
                                   4.0 * z[0] * z[2] + 8.0 * z[1] * z[2] + 4.0 * z[2] * z[2]
                                       - 4.0 * z[2],
                                   -8.0 * z[0] * z[1] - 8.0 * z[1] * z[1] - 4.0 * z[1] * z[2]
                                       + 8.0 * z[1]);
    case 17:
      return VectorFieldT<T, sDim>(4.0 * z[1] * z[2],
                                   -8.0 * z[0] * z[2] - 4.0 * z[1] * z[2] - 8.0 * z[2] * z[2]
                                       + 8.0 * z[2],
                                   4.0 * z[0] * z[1] + 4.0 * z[1] * z[1] + 8.0 * z[1] * z[2]
                                       - 4.0 * z[1]);
    case 18:
      return VectorFieldT<T, sDim>(8.0 * z[0] * z[2] + 4.0 * z[1] * z[2] + 4.0 * z[2] * z[2]
                                       - 4.0 * z[2],
                                   4.0 * z[0] * z[2],
                                   -8.0 * z[0] * z[0] - 8.0 * z[0] * z[1] - 4.0 * z[0] * z[2]
                                       + 8.0 * z[0]);
    case 19:
      return VectorFieldT<T, sDim>(-4.0 * z[0] * z[2] - 8.0 * z[1] * z[2] - 8.0 * z[2] * z[2]
                                       + 8.0 * z[2],
                                   4.0 * z[0] * z[2],
                                   4.0 * z[0] * z[0] + 4.0 * z[0] * z[1] + 8.0 * z[0] * z[2]
                                       - 4.0 * z[0]);
    default:
      THROW("Not implemented")
    }
  }
};

using P2TetCurl = P2TetCurlT<>;


#endif // CURLSHAPESTETRAHEDRON_HPP
