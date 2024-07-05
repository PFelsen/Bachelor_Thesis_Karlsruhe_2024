#include "Shapes.hpp"

#include "ArgyrisShapes.hpp"
#include "CRShapes.hpp"
#include "CurlShapes.hpp"
#include "LagrangeShapes.hpp"
#include "RTShapes.hpp"
#include "SerendipityShapes.hpp"

class P1rotquad : public Shape {
public:
  const std::string Name() const { return "P1rotquad"; }

  double operator()(const Point &z, int i) const {
    switch (i) {
    case 0:
      return -z[0] * z[0] + z[1] * z[1] + z[0] - 2 * z[1] + 0.75;
    case 1:
      return z[0] * z[0] - z[1] * z[1] + z[1] - 0.25;
    case 2:
      return -z[0] * z[0] + z[1] * z[1] + z[0] - 0.25;
    case 3:
      return z[0] * z[0] - z[1] * z[1] - 2 * z[0] + z[1] + 0.75;
    }
    THROW("Not implemented")
    return 0.0;
  }

  VectorField LocalGradient(const Point &z, int i) const {
    switch (i) {

    case 0:
      return VectorField(-2 * z[0] + 1, 2 * z[1] - 2);
    case 1:
      return VectorField(2 * z[0], -2 * z[1] + 1);
    case 2:
      return VectorField(-2 * z[0] + 1, 2 * z[1]);
    case 3:
      return VectorField(2 * z[0] - 2, -2 * z[1] + 1);
    }
    THROW("Not implemented")
    return zero;
  }

  //  P1rotquad(const Quadrature &Q,
  //            const Quadrature &FaceQ = GetQuadrature("Qint5"))
  P1rotquad() : Shape(4) {
    //    fill();
  }

  void NodalPoints(const Cell &c, Point *z) const {
    for (int i = 0; i < c.Faces(); ++i)
      z[i] = c.Face(i);
  }
};

class P1rotaltquad : public Shape {
public:
  const std::string Name() const { return "P1rotaltquad"; }

  double operator()(const Point &z, int i) const {
    switch (i) {
    case 0:
      return 0.5 * (3.0 / 2 + 3 * z[0] - 5 * z[1] - 3 * z[0] * z[0] + 3 * z[1] * z[1]);
    case 1:
      return 0.5 * (-0.5 - z[0] + 3 * z[1] + 3 * z[0] * z[0] - 3 * z[1] * z[1]);
    case 2:
      return 0.5 * (-0.5 + 3 * z[0] - z[1] - 3 * z[0] * z[0] + 3 * z[1] * z[1]);
    case 3:
      return 0.5 * (3.0 / 2 - 5 * z[0] + 3 * z[1] + 3 * z[0] * z[0] - 3 * z[1] * z[1]);
    }
    THROW("Not implemented")
    return 0.0;
  }

  VectorField LocalGradient(const Point &z, int i) const {
    switch (i) {
    case 0:
      return VectorField(1.5 - 3 * z[0], -2.5 + 3 * z[1]);
    case 1:
      return VectorField(-0.5 + 3 * z[0], 1.5 - 3 * z[1]);
    case 2:
      return VectorField(1.5 - 3 * z[0], -0.5 + 3 * z[1]);
    case 3:
      return VectorField(-2.5 + 3 * z[0], 1.5 - 3 * z[1]);
    }
    THROW("Not implemented")
    return zero;
  }

  //  P1rotaltquad(const Quadrature &Q,
  //               const Quadrature &FaceQ = GetQuadrature("Qint5"))
  P1rotaltquad() : Shape(4) {
    //    fill();
  }

  void NodalPoints(const Cell &c, Point *z) const {
    for (int i = 0; i < c.Faces(); ++i)
      z[i] = c.Face(i);
  }
};
