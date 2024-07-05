#ifndef LAGRANGESHAPESTETRAHEDRON_HPP
#define LAGRANGESHAPESTETRAHEDRON_HPP

#include "LagrangeShape.hpp"
#include "SimplexPolynomials.hpp"

// TODO: recalcNodalPoints, recalc, recalcD

template<typename T = double, int sDim = SpaceDimension, int tDim = TimeDimension>
class TetrahedronShapeT : public LagrangeShapeT<T, sDim, tDim> {
  static int calcNumberOfPoints(int polynomialDegree) {
    return ((polynomialDegree + 1) * (2 + (polynomialDegree + 1) * (polynomialDegree + 4))) / 6;
  }
public:
  int polynomialDegree;

  TetrahedronShapeT(int polynomialDegree) :
      LagrangeShapeT<T, sDim, tDim>(calcNumberOfPoints(polynomialDegree)),
      polynomialDegree(polynomialDegree) {}

  void NodalPoints(const Cell &c, vector<PointT<T, sDim, tDim>> &z) const override {
    z.resize(this->numNodalPoints);
    LagrangeNodalPointsTetrahedron(c, z, polynomialDegree);
  }
};

typedef TetrahedronShapeT<> TetrahedronShape;

template<typename T = double, int sDim = SpaceDimension, int tDim = TimeDimension>
class P0TetT : public TetrahedronShapeT<T, sDim, tDim> {
public:
  const std::string Name() const override { return "P0Tet"; }

  T operator()(const PointT<T, sDim, tDim> &z, int i) const override { return T(1.0); }

  VectorFieldT<T, sDim> LocalGradient(const PointT<T, sDim, tDim> &z, int i) const override {
    return VectorFieldT<T, sDim>();
  }

  explicit P0TetT() : TetrahedronShapeT<T, sDim, tDim>(0) {}
};

typedef P0TetT<> P0Tet;

template<typename T = double, int sDim = SpaceDimension, int tDim = TimeDimension>
class P1TetT : public TetrahedronShapeT<T, sDim, tDim> {
public:
  const std::string Name() const override { return "P1Tet"; }

  T operator()(const PointT<T, sDim, tDim> &z, int i) const override {
    switch (i) {
    case 0:
      return (1 - z[0] - z[1] - z[2]);
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

  VectorFieldT<T, sDim> LocalGradient(const PointT<T, sDim, tDim> &z, int i) const override {
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

  explicit P1TetT() : TetrahedronShapeT<T, sDim, tDim>(1) {}
};

typedef P1TetT<> P1Tet;

template<typename T = double, int sDim = SpaceDimension, int tDim = TimeDimension>
class P2TetT : public TetrahedronShapeT<T, sDim, tDim> {
  T (*p2Fun)(T, int) = &poly::shapeFunction<T, 2>;

  T (*p2Der)(T, int) = &poly::shapeDerivative<T, 2>;
public:
  const std::string Name() const { return "P2Tet"; }

  T operator()(const PointT<T, sDim, tDim> &z, int i) const {
    T psi[4]{1 - z[0] - z[1] - z[2], z[0], z[1], z[2]};
    switch (i) {
    case 0:
      return p2Fun(psi[0], 2);
    case 2:
      return p2Fun(psi[1], 2);
    case 5:
      return p2Fun(psi[2], 2);
    case 9:
      return p2Fun(psi[3], 2);
    case 1:
      return p2Fun(psi[0], 1) * p2Fun(psi[1], 1);
    case 4:
      return p2Fun(psi[1], 1) * p2Fun(psi[2], 1);
    case 3:
      return p2Fun(psi[0], 1) * p2Fun(psi[2], 1);
    case 6:
      return p2Fun(psi[0], 1) * p2Fun(psi[3], 1);
    case 7:
      return p2Fun(psi[1], 1) * p2Fun(psi[3], 1);
    case 8:
      return p2Fun(psi[2], 1) * p2Fun(psi[3], 1);
    default:
      THROW("Not implemented")
    }
  }

  VectorFieldT<T, sDim> LocalGradient(const PointT<T, sDim, tDim> &z, int i) const {
    T psi[4]{1 - z[0] - z[1] - z[2], z[0], z[1], z[2]};
    switch (i) {
    case 0:
      return VectorFieldT<T, sDim>(-p2Der(psi[0], 2), -p2Der(psi[0], 2), -p2Der(psi[0], 2));
    case 2:
      return VectorFieldT<T, sDim>(p2Der(psi[1], 2), T(0.0), T(0.0));
    case 5:
      return VectorFieldT<T, sDim>(T(0.0), p2Der(psi[2], 2), T(0.0));
    case 9:
      return VectorFieldT<T, sDim>(T(0.0), T(0.0), p2Der(psi[3], 2));
    case 1:
      return VectorFieldT<T, sDim>(p2Fun(psi[0], 1) * p2Der(psi[1], 1)
                                       - p2Der(psi[0], 1) * p2Fun(psi[1], 1),
                                   -p2Der(psi[0], 1) * p2Fun(psi[1], 1),
                                   -p2Der(psi[0], 1) * p2Fun(psi[1], 1));
    case 4:
      return VectorFieldT<T, sDim>(p2Der(psi[1], 1) * p2Fun(psi[2], 1),
                                   p2Fun(psi[1], 1) * p2Der(psi[2], 1), T(0.0));
    case 3:
      return VectorFieldT<T, sDim>(-p2Der(psi[0], 1) * p2Fun(psi[2], 1),
                                   p2Fun(psi[0], 1) * p2Der(psi[2], 1)
                                       - p2Der(psi[0], 1) * p2Fun(psi[2], 1),
                                   -p2Der(psi[0], 1) * p2Fun(psi[2], 1));
    case 6:
      return VectorFieldT<T, sDim>(-p2Der(psi[0], 1) * p2Fun(psi[3], 1),
                                   -p2Der(psi[0], 1) * p2Fun(psi[3], 1),
                                   p2Fun(psi[0], 1) * p2Der(psi[3], 1)
                                       - p2Der(psi[0], 1) * p2Fun(psi[3], 1));
    case 7:
      return VectorFieldT<T, sDim>(p2Der(psi[1], 1) * p2Fun(psi[3], 1), T(0.0),
                                   p2Fun(psi[1], 1) * p2Der(psi[3], 1));
    case 8:
      return VectorFieldT<T, sDim>(T(0.0), p2Der(psi[2], 1) * p2Fun(psi[3], 1),
                                   p2Fun(psi[2], 1) * p2Der(psi[3], 1));
    default:
      THROW("Not implemented")
    }
  }

  explicit P2TetT() : TetrahedronShapeT<T, sDim, tDim>(2) {}
};

typedef P2TetT<> P2Tet;

template<typename T = double, int sDim = SpaceDimension, int tDim = TimeDimension>
class P3TetT : public TetrahedronShapeT<T, sDim, tDim> {
  T (*p3Fun)(T, int) = &poly::shapeFunction<T, 3>;

  T (*p3Der)(T, int) = &poly::shapeDerivative<T, 3>;
public:
  const std::string Name() const { return "P3Tet"; }

  T operator()(const PointT<T, sDim, tDim> &z, int i) const {
    T psi[4]{1 - z[0] - z[1] - z[2], z[0], z[1], z[2]};
    switch (i) {
    // Corners
    case 0:
      return p3Fun(psi[0], 3);
    case 3:
      return p3Fun(psi[1], 3);
    case 9:
      return p3Fun(psi[2], 3);
    case 19:
      return p3Fun(psi[3], 3);
      // Edges
    case 1:
      return p3Fun(psi[0], 2) * p3Fun(psi[1], 1);
    case 2:
      return p3Fun(psi[0], 1) * p3Fun(psi[1], 2);
    case 6:
      return p3Fun(psi[1], 2) * p3Fun(psi[2], 1);
    case 8:
      return p3Fun(psi[1], 1) * p3Fun(psi[2], 2);
    case 4:
      return p3Fun(psi[0], 2) * p3Fun(psi[2], 1);
    case 7:
      return p3Fun(psi[0], 1) * p3Fun(psi[2], 2);

    case 10:
      return p3Fun(psi[0], 2) * p3Fun(psi[3], 1);
    case 16:
      return p3Fun(psi[0], 1) * p3Fun(psi[3], 2);
    case 12:
      return p3Fun(psi[1], 2) * p3Fun(psi[3], 1);
    case 17:
      return p3Fun(psi[1], 1) * p3Fun(psi[3], 2);
    case 15:
      return p3Fun(psi[2], 2) * p3Fun(psi[3], 1);
    case 18:
      return p3Fun(psi[2], 1) * p3Fun(psi[3], 2);
      // Faces
    case 5:
      return p3Fun(psi[0], 1) * p3Fun(psi[1], 1) * p3Fun(psi[2], 1);
    case 11:
      return p3Fun(psi[0], 1) * p3Fun(psi[1], 1) * p3Fun(psi[3], 1);
    case 13:
      return p3Fun(psi[0], 1) * p3Fun(psi[2], 1) * p3Fun(psi[3], 1);
    case 14:
      return p3Fun(psi[1], 1) * p3Fun(psi[2], 1) * p3Fun(psi[3], 1);

    default:
      THROW("Not implemented")
    }
  }

  VectorFieldT<T, sDim> LocalGradient(const PointT<T, sDim, tDim> &z, int i) const {
    T psi[4]{1 - z[0] - z[1] - z[2], z[0], z[1], z[2]};
    switch (i) {
    // Corners
    case 0:
      return VectorFieldT<T, sDim>(-p3Der(psi[0], 3), -p3Der(psi[0], 3), -p3Der(psi[0], 3));
    case 3:
      return VectorFieldT<T, sDim>(p3Der(psi[1], 3), T(0.0), T(0.0));
    case 9:
      return VectorFieldT<T, sDim>(T(0.0), p3Der(psi[2], 3), T(0.0));
    case 19:
      return VectorFieldT<T, sDim>(T(0.0), T(0.0), p3Der(psi[3], 3));
      // Edges
    case 1:
      return VectorFieldT<T, sDim>(p3Fun(psi[0], 2) * p3Der(psi[1], 1)
                                       - p3Der(psi[0], 2) * p3Fun(psi[1], 1),
                                   -p3Der(psi[0], 2) * p3Fun(psi[1], 1),
                                   -p3Der(psi[0], 2) * p3Fun(psi[1], 1));
    case 2:
      return VectorFieldT<T, sDim>(p3Fun(psi[0], 1) * p3Der(psi[1], 2)
                                       - p3Der(psi[0], 1) * p3Fun(psi[1], 2),
                                   -p3Der(psi[0], 1) * p3Fun(psi[1], 2),
                                   -p3Der(psi[0], 1) * p3Fun(psi[1], 2));
    case 6:
      return VectorFieldT<T, sDim>(p3Der(psi[1], 2) * p3Fun(psi[2], 1),
                                   p3Fun(psi[1], 2) * p3Der(psi[2], 1), T(0.0));
    case 8:
      return VectorFieldT<T, sDim>(p3Der(psi[1], 1) * p3Fun(psi[2], 2),
                                   p3Fun(psi[1], 1) * p3Der(psi[2], 2), T(0.0));
    case 4:
      return VectorFieldT<T, sDim>(-p3Der(psi[0], 2) * p3Fun(psi[2], 1),
                                   p3Fun(psi[0], 2) * p3Der(psi[2], 1)
                                       - p3Der(psi[0], 2) * p3Fun(psi[2], 1),
                                   -p3Der(psi[0], 2) * p3Fun(psi[2], 1));
    case 7:
      return VectorFieldT<T, sDim>(-p3Der(psi[0], 1) * p3Fun(psi[2], 2),
                                   p3Fun(psi[0], 1) * p3Der(psi[2], 2)
                                       - p3Der(psi[0], 1) * p3Fun(psi[2], 2),
                                   -p3Der(psi[0], 1) * p3Fun(psi[2], 2));

    case 10:
      return VectorFieldT<T, sDim>(-p3Der(psi[0], 2) * p3Fun(psi[3], 1),
                                   -p3Der(psi[0], 2) * p3Fun(psi[3], 1),
                                   p3Fun(psi[0], 2) * p3Der(psi[3], 1)
                                       - p3Der(psi[0], 2) * p3Fun(psi[3], 1));
    case 16:
      return VectorFieldT<T, sDim>(-p3Der(psi[0], 1) * p3Fun(psi[3], 2),
                                   -p3Der(psi[0], 1) * p3Fun(psi[3], 2),
                                   p3Fun(psi[0], 1) * p3Der(psi[3], 2)
                                       - p3Der(psi[0], 1) * p3Fun(psi[3], 2));
    case 12:
      return VectorFieldT<T, sDim>(p3Der(psi[1], 2) * p3Fun(psi[3], 1), T(0.0),
                                   p3Fun(psi[1], 2) * p3Der(psi[3], 1));
    case 17:
      return VectorFieldT<T, sDim>(p3Der(psi[1], 1) * p3Fun(psi[3], 2), T(0.0),
                                   p3Fun(psi[1], 1) * p3Der(psi[3], 2));
    case 15:
      return VectorFieldT<T, sDim>(T(0.0), p3Der(psi[2], 2) * p3Fun(psi[3], 1),
                                   p3Fun(psi[2], 2) * p3Der(psi[3], 1));
    case 18:
      return VectorFieldT<T, sDim>(T(0.0), p3Der(psi[2], 1) * p3Fun(psi[3], 2),
                                   p3Fun(psi[2], 1) * p3Der(psi[3], 2));
      // Faces
    case 5:
      return VectorFieldT<T, sDim>(p3Fun(psi[0], 1) * p3Der(psi[1], 1) * p3Fun(psi[2], 1)
                                       - p3Der(psi[0], 1) * p3Fun(psi[1], 1) * p3Fun(psi[2], 1),
                                   p3Fun(psi[0], 1) * p3Fun(psi[1], 1) * p3Der(psi[2], 1)
                                       - p3Der(psi[0], 1) * p3Fun(psi[1], 1) * p3Fun(psi[2], 1),
                                   -p3Der(psi[0], 1) * p3Fun(psi[1], 1) * p3Fun(psi[2], 1));
    case 11:
      return VectorFieldT<T, sDim>(p3Fun(psi[0], 1) * p3Der(psi[1], 1) * p3Fun(psi[3], 1)
                                       - p3Der(psi[0], 1) * p3Fun(psi[1], 1) * p3Fun(psi[3], 1),
                                   -p3Der(psi[0], 1) * p3Fun(psi[1], 1) * p3Fun(psi[3], 1),
                                   p3Fun(psi[0], 1) * p3Fun(psi[1], 1) * p3Der(psi[3], 1)
                                       - p3Der(psi[0], 1) * p3Fun(psi[1], 1) * p3Fun(psi[3], 1));
    case 13:
      return VectorFieldT<T, sDim>(-p3Der(psi[0], 1) * p3Fun(psi[2], 1) * p3Fun(psi[3], 1),
                                   p3Fun(psi[0], 1) * p3Der(psi[2], 1) * p3Fun(psi[3], 1)
                                       - p3Der(psi[0], 1) * p3Fun(psi[2], 1) * p3Fun(psi[3], 1),
                                   p3Fun(psi[0], 1) * p3Fun(psi[2], 1) * p3Der(psi[3], 1)
                                       - p3Der(psi[0], 1) * p3Fun(psi[2], 1) * p3Fun(psi[3], 1));
    case 14:
      return VectorFieldT<T, sDim>(p3Der(psi[1], 1) * p3Fun(psi[2], 1) * p3Fun(psi[3], 1),
                                   p3Fun(psi[1], 1) * p3Der(psi[2], 1) * p3Fun(psi[3], 1),
                                   p3Fun(psi[1], 1) * p3Fun(psi[2], 1) * p3Der(psi[3], 1));

    default:
      THROW("Not implemented")
    }
  }

  explicit P3TetT() : TetrahedronShapeT<T, sDim, tDim>(3) {}
};

typedef P3TetT<> P3Tet;

template<typename T = double, int sDim = SpaceDimension, int tDim = TimeDimension>
class P4TetT : public TetrahedronShapeT<T, sDim, tDim> {
  T (*p4Fun)(T, int) = &poly::shapeFunction<T, 4>;

  T (*p4Der)(T, int) = &poly::shapeDerivative<T, 4>;
public:
  const std::string Name() const { return "P4Tet"; }

  T operator()(const PointT<T, sDim, tDim> &z, int i) const {
    T psi[4]{1 - z[0] - z[1] - z[2], z[0], z[1], z[2]};
    switch (i) {
    // Corners
    case 0:
      return p4Fun(psi[0], 4);
    case 4:
      return p4Fun(psi[1], 4);
    case 14:
      return p4Fun(psi[2], 4);
    case 34:
      return p4Fun(psi[3], 4);
      // Edges
    case 1:
      return p4Fun(psi[0], 3) * p4Fun(psi[1], 1);
    case 2:
      return p4Fun(psi[0], 2) * p4Fun(psi[1], 2);
    case 3:
      return p4Fun(psi[0], 1) * p4Fun(psi[1], 3);
    case 8:
      return p4Fun(psi[1], 3) * p4Fun(psi[2], 1);
    case 11:
      return p4Fun(psi[1], 2) * p4Fun(psi[2], 2);
    case 13:
      return p4Fun(psi[1], 1) * p4Fun(psi[2], 3);
    case 5:
      return p4Fun(psi[0], 3) * p4Fun(psi[2], 1);
    case 9:
      return p4Fun(psi[0], 2) * p4Fun(psi[2], 2);
    case 12:
      return p4Fun(psi[0], 1) * p4Fun(psi[2], 3);

    case 15:
      return p4Fun(psi[0], 3) * p4Fun(psi[3], 1);
    case 25:
      return p4Fun(psi[0], 2) * p4Fun(psi[3], 2);
    case 31:
      return p4Fun(psi[0], 1) * p4Fun(psi[3], 3);
    case 18:
      return p4Fun(psi[1], 3) * p4Fun(psi[3], 1);
    case 27:
      return p4Fun(psi[1], 2) * p4Fun(psi[3], 2);
    case 32:
      return p4Fun(psi[1], 1) * p4Fun(psi[3], 3);
    case 24:
      return p4Fun(psi[2], 3) * p4Fun(psi[3], 1);
    case 30:
      return p4Fun(psi[2], 2) * p4Fun(psi[3], 2);
    case 33:
      return p4Fun(psi[2], 1) * p4Fun(psi[3], 3);
      // Faces
    case 6:
      return p4Fun(psi[0], 2) * p4Fun(psi[1], 1) * p4Fun(psi[2], 1);
    case 7:
      return p4Fun(psi[0], 1) * p4Fun(psi[1], 2) * p4Fun(psi[2], 1);
    case 10:
      return p4Fun(psi[0], 1) * p4Fun(psi[1], 1) * p4Fun(psi[2], 2);

    case 16:
      return p4Fun(psi[0], 2) * p4Fun(psi[1], 1) * p4Fun(psi[3], 1);
    case 17:
      return p4Fun(psi[0], 1) * p4Fun(psi[1], 2) * p4Fun(psi[3], 1);
    case 26:
      return p4Fun(psi[0], 1) * p4Fun(psi[1], 1) * p4Fun(psi[3], 2);

    case 19:
      return p4Fun(psi[0], 2) * p4Fun(psi[2], 1) * p4Fun(psi[3], 1);
    case 22:
      return p4Fun(psi[0], 1) * p4Fun(psi[2], 2) * p4Fun(psi[3], 1);
    case 28:
      return p4Fun(psi[0], 1) * p4Fun(psi[2], 1) * p4Fun(psi[3], 2);

    case 21:
      return p4Fun(psi[1], 2) * p4Fun(psi[2], 1) * p4Fun(psi[3], 1);
    case 23:
      return p4Fun(psi[1], 1) * p4Fun(psi[2], 2) * p4Fun(psi[3], 1);
    case 29:
      return p4Fun(psi[1], 1) * p4Fun(psi[2], 1) * p4Fun(psi[3], 2);

      // Inner
    case 20:
      return p4Fun(psi[0], 1) * p4Fun(psi[1], 1) * p4Fun(psi[2], 1) * p4Fun(psi[3], 1);

    default:
      THROW("Not implemented")
    }
  }

  VectorFieldT<T, sDim> LocalGradient(const PointT<T, sDim, tDim> &z, int i) const {
    T psi[4]{1 - z[0] - z[1] - z[2], z[0], z[1], z[2]};
    switch (i) {
    // Corners
    case 0:
      return VectorFieldT<T, sDim>(-p4Der(psi[0], 4), -p4Der(psi[0], 4), -p4Der(psi[0], 4));
    case 4:
      return VectorFieldT<T, sDim>(p4Der(psi[1], 4), T(0.0), T(0.0));
    case 14:
      return VectorFieldT<T, sDim>(T(0.0), p4Der(psi[2], 4), T(0.0));
    case 34:
      return VectorFieldT<T, sDim>(T(0.0), T(0.0), p4Der(psi[3], 4));
      // Edges
    case 1:
      return VectorFieldT<T, sDim>(p4Fun(psi[0], 3) * p4Der(psi[1], 1)
                                       - p4Der(psi[0], 3) * p4Fun(psi[1], 1),
                                   -p4Der(psi[0], 3) * p4Fun(psi[1], 1),
                                   -p4Der(psi[0], 3) * p4Fun(psi[1], 1));
    case 2:
      return VectorFieldT<T, sDim>(p4Fun(psi[0], 2) * p4Der(psi[1], 2)
                                       - p4Der(psi[0], 2) * p4Fun(psi[1], 2),
                                   -p4Der(psi[0], 2) * p4Fun(psi[1], 2),
                                   -p4Der(psi[0], 2) * p4Fun(psi[1], 2));
    case 3:
      return VectorFieldT<T, sDim>(p4Fun(psi[0], 1) * p4Der(psi[1], 3)
                                       - p4Der(psi[0], 1) * p4Fun(psi[1], 3),
                                   -p4Der(psi[0], 1) * p4Fun(psi[1], 3),
                                   -p4Der(psi[0], 1) * p4Fun(psi[1], 3));

    case 8:
      return VectorFieldT<T, sDim>(p4Der(psi[1], 3) * p4Fun(psi[2], 1),
                                   p4Fun(psi[1], 3) * p4Der(psi[2], 1), T(0.0));
    case 11:
      return VectorFieldT<T, sDim>(p4Der(psi[1], 2) * p4Fun(psi[2], 2),
                                   p4Fun(psi[1], 2) * p4Der(psi[2], 2), T(0.0));
    case 13:
      return VectorFieldT<T, sDim>(p4Der(psi[1], 1) * p4Fun(psi[2], 3),
                                   p4Fun(psi[1], 1) * p4Der(psi[2], 3), T(0.0));

    case 5:
      return VectorFieldT<T, sDim>(-p4Der(psi[0], 3) * p4Fun(psi[2], 1),
                                   p4Fun(psi[0], 3) * p4Der(psi[2], 1)
                                       - p4Der(psi[0], 3) * p4Fun(psi[2], 1),
                                   -p4Der(psi[0], 3) * p4Fun(psi[2], 1));
    case 9:
      return VectorFieldT<T, sDim>(-p4Der(psi[0], 2) * p4Fun(psi[2], 2),
                                   p4Fun(psi[0], 2) * p4Der(psi[2], 2)
                                       - p4Der(psi[0], 2) * p4Fun(psi[2], 2),
                                   -p4Der(psi[0], 2) * p4Fun(psi[2], 2));
    case 12:
      return VectorFieldT<T, sDim>(-p4Der(psi[0], 1) * p4Fun(psi[2], 3),
                                   p4Fun(psi[0], 1) * p4Der(psi[2], 3)
                                       - p4Der(psi[0], 1) * p4Fun(psi[2], 3),
                                   -p4Der(psi[0], 1) * p4Fun(psi[2], 3));

    case 15:
      return VectorFieldT<T, sDim>(-p4Der(psi[0], 3) * p4Fun(psi[3], 1),
                                   -p4Der(psi[0], 3) * p4Fun(psi[3], 1),
                                   p4Fun(psi[0], 3) * p4Der(psi[3], 1)
                                       - p4Der(psi[0], 3) * p4Fun(psi[3], 1));
    case 25:
      return VectorFieldT<T, sDim>(-p4Der(psi[0], 2) * p4Fun(psi[3], 2),
                                   -p4Der(psi[0], 2) * p4Fun(psi[3], 2),
                                   p4Fun(psi[0], 2) * p4Der(psi[3], 2)
                                       - p4Der(psi[0], 2) * p4Fun(psi[3], 2));
    case 31:
      return VectorFieldT<T, sDim>(-p4Der(psi[0], 1) * p4Fun(psi[3], 3),
                                   -p4Der(psi[0], 1) * p4Fun(psi[3], 3),
                                   p4Fun(psi[0], 1) * p4Der(psi[3], 3)
                                       - p4Der(psi[0], 1) * p4Fun(psi[3], 3));
    case 18:
      return VectorFieldT<T, sDim>(p4Der(psi[1], 3) * p4Fun(psi[3], 1), T(0.0),
                                   p4Fun(psi[1], 3) * p4Der(psi[3], 1));
    case 27:
      return VectorFieldT<T, sDim>(p4Der(psi[1], 2) * p4Fun(psi[3], 2), T(0.0),
                                   p4Fun(psi[1], 2) * p4Der(psi[3], 2));
    case 32:
      return VectorFieldT<T, sDim>(p4Der(psi[1], 1) * p4Fun(psi[3], 3), T(0.0),
                                   p4Fun(psi[1], 1) * p4Der(psi[3], 3));
    case 24:
      return VectorFieldT<T, sDim>(T(0.0), p4Der(psi[2], 3) * p4Fun(psi[3], 1),
                                   p4Fun(psi[2], 3) * p4Der(psi[3], 1));
    case 30:
      return VectorFieldT<T, sDim>(T(0.0), p4Der(psi[2], 2) * p4Fun(psi[3], 2),
                                   p4Fun(psi[2], 2) * p4Der(psi[3], 2));
    case 33:
      return VectorFieldT<T, sDim>(T(0.0), p4Der(psi[2], 1) * p4Fun(psi[3], 3),
                                   p4Fun(psi[2], 1) * p4Der(psi[3], 3));
      // Faces
    case 6:
      return VectorFieldT<T, sDim>(p4Fun(psi[0], 2) * p4Der(psi[1], 1) * p4Fun(psi[2], 1)
                                       - p4Der(psi[0], 2) * p4Fun(psi[1], 1) * p4Fun(psi[2], 1),
                                   p4Fun(psi[0], 2) * p4Fun(psi[1], 1) * p4Der(psi[2], 1)
                                       - p4Der(psi[0], 2) * p4Fun(psi[1], 1) * p4Fun(psi[2], 1),
                                   -p4Der(psi[0], 2) * p4Fun(psi[1], 1) * p4Fun(psi[2], 1));
    case 7:
      return VectorFieldT<T, sDim>(p4Fun(psi[0], 1) * p4Der(psi[1], 2) * p4Fun(psi[2], 1)
                                       - p4Der(psi[0], 1) * p4Fun(psi[1], 2) * p4Fun(psi[2], 1),
                                   p4Fun(psi[0], 1) * p4Fun(psi[1], 2) * p4Der(psi[2], 1)
                                       - p4Der(psi[0], 1) * p4Fun(psi[1], 2) * p4Fun(psi[2], 1),
                                   -p4Der(psi[0], 1) * p4Fun(psi[1], 2) * p4Fun(psi[2], 1));
    case 10:
      return VectorFieldT<T, sDim>(p4Fun(psi[0], 1) * p4Der(psi[1], 1) * p4Fun(psi[2], 2)
                                       - p4Der(psi[0], 1) * p4Fun(psi[1], 1) * p4Fun(psi[2], 2),
                                   p4Fun(psi[0], 1) * p4Fun(psi[1], 1) * p4Der(psi[2], 2)
                                       - p4Der(psi[0], 1) * p4Fun(psi[1], 1) * p4Fun(psi[2], 2),
                                   -p4Der(psi[0], 1) * p4Fun(psi[1], 1) * p4Fun(psi[2], 2));
    case 16:
      return VectorFieldT<T, sDim>(p4Fun(psi[0], 2) * p4Der(psi[1], 1) * p4Fun(psi[3], 1)
                                       - p4Der(psi[0], 2) * p4Fun(psi[1], 1) * p4Fun(psi[3], 1),
                                   -p4Der(psi[0], 2) * p4Fun(psi[1], 1) * p4Fun(psi[3], 1),
                                   p4Fun(psi[0], 2) * p4Fun(psi[1], 1) * p4Der(psi[3], 1)
                                       - p4Der(psi[0], 2) * p4Fun(psi[1], 1) * p4Fun(psi[3], 1));
    case 17:
      return VectorFieldT<T, sDim>(p4Fun(psi[0], 1) * p4Der(psi[1], 2) * p4Fun(psi[3], 1)
                                       - p4Der(psi[0], 1) * p4Fun(psi[1], 2) * p4Fun(psi[3], 1),
                                   -p4Der(psi[0], 1) * p4Fun(psi[1], 2) * p4Fun(psi[3], 1),
                                   p4Fun(psi[0], 1) * p4Fun(psi[1], 2) * p4Der(psi[3], 1)
                                       - p4Der(psi[0], 1) * p4Fun(psi[1], 2) * p4Fun(psi[3], 1));
    case 26:
      return VectorFieldT<T, sDim>(p4Fun(psi[0], 1) * p4Der(psi[1], 1) * p4Fun(psi[3], 2)
                                       - p4Der(psi[0], 1) * p4Fun(psi[1], 1) * p4Fun(psi[3], 2),
                                   -p4Der(psi[0], 1) * p4Fun(psi[1], 1) * p4Fun(psi[3], 2),
                                   p4Fun(psi[0], 1) * p4Fun(psi[1], 1) * p4Der(psi[3], 2)
                                       - p4Der(psi[0], 1) * p4Fun(psi[1], 1) * p4Fun(psi[3], 2));
    case 19:
      return VectorFieldT<T, sDim>(-p4Der(psi[0], 2) * p4Fun(psi[2], 1) * p4Fun(psi[3], 1),
                                   p4Fun(psi[0], 2) * p4Der(psi[2], 1) * p4Fun(psi[3], 1)
                                       - p4Der(psi[0], 2) * p4Fun(psi[2], 1) * p4Fun(psi[3], 1),
                                   p4Fun(psi[0], 2) * p4Fun(psi[2], 1) * p4Der(psi[3], 1)
                                       - p4Der(psi[0], 2) * p4Fun(psi[2], 1) * p4Fun(psi[3], 1));
    case 22:
      return VectorFieldT<T, sDim>(-p4Der(psi[0], 1) * p4Fun(psi[2], 2) * p4Fun(psi[3], 1),
                                   p4Fun(psi[0], 1) * p4Der(psi[2], 2) * p4Fun(psi[3], 1)
                                       - p4Der(psi[0], 1) * p4Fun(psi[2], 2) * p4Fun(psi[3], 1),
                                   p4Fun(psi[0], 1) * p4Fun(psi[2], 2) * p4Der(psi[3], 1)
                                       - p4Der(psi[0], 1) * p4Fun(psi[2], 2) * p4Fun(psi[3], 1));
    case 28:
      return VectorFieldT<T, sDim>(-p4Der(psi[0], 1) * p4Fun(psi[2], 1) * p4Fun(psi[3], 2),
                                   p4Fun(psi[0], 1) * p4Der(psi[2], 1) * p4Fun(psi[3], 2)
                                       - p4Der(psi[0], 1) * p4Fun(psi[2], 1) * p4Fun(psi[3], 2),
                                   p4Fun(psi[0], 1) * p4Fun(psi[2], 1) * p4Der(psi[3], 2)
                                       - p4Der(psi[0], 1) * p4Fun(psi[2], 1) * p4Fun(psi[3], 2));
    case 21:
      return VectorFieldT<T, sDim>(p4Der(psi[1], 2) * p4Fun(psi[2], 1) * p4Fun(psi[3], 1),
                                   p4Fun(psi[1], 2) * p4Der(psi[2], 1) * p4Fun(psi[3], 1),
                                   p4Fun(psi[1], 2) * p4Fun(psi[2], 1) * p4Der(psi[3], 1));
    case 23:
      return VectorFieldT<T, sDim>(p4Der(psi[1], 1) * p4Fun(psi[2], 2) * p4Fun(psi[3], 1),
                                   p4Fun(psi[1], 1) * p4Der(psi[2], 2) * p4Fun(psi[3], 1),
                                   p4Fun(psi[1], 1) * p4Fun(psi[2], 2) * p4Der(psi[3], 1));
    case 29:
      return VectorFieldT<T, sDim>(p4Der(psi[1], 1) * p4Fun(psi[2], 1) * p4Fun(psi[3], 2),
                                   p4Fun(psi[1], 1) * p4Der(psi[2], 1) * p4Fun(psi[3], 2),
                                   p4Fun(psi[1], 1) * p4Fun(psi[2], 1) * p4Der(psi[3], 2));

      // Inner
    case 20:
      return VectorFieldT<
          T, sDim>(p4Fun(psi[0], 1) * p4Der(psi[1], 1) * p4Fun(psi[2], 1) * p4Fun(psi[3], 1)
                       - p4Der(psi[0], 1) * p4Fun(psi[1], 1) * p4Fun(psi[2], 1) * p4Fun(psi[3], 1),
                   p4Fun(psi[0], 1) * p4Fun(psi[1], 1) * p4Der(psi[2], 1) * p4Fun(psi[3], 1)
                       - p4Der(psi[0], 1) * p4Fun(psi[1], 1) * p4Fun(psi[2], 1) * p4Fun(psi[3], 1),
                   p4Fun(psi[0], 1) * p4Fun(psi[1], 1) * p4Fun(psi[2], 1) * p4Der(psi[3], 1)
                       - p4Der(psi[0], 1) * p4Fun(psi[1], 1) * p4Fun(psi[2], 1) * p4Fun(psi[3], 1));
    default:
      THROW("Not implemented")
    }
  }

  explicit P4TetT() : TetrahedronShapeT<T, sDim, tDim>(4) {}
};

typedef P4TetT<> P4Tet;

#endif // LAGRANGESHAPESTETRAHEDRON_HPP
