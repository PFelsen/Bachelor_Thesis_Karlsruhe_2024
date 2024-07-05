#ifndef ARGYRISSHAPES_HPP
#define ARGYRISSHAPES_HPP

#include "Shape.hpp"

template<typename T = double, int sDim = SpaceDimension, int tDim = TimeDimension>
class ArgyrisTriT : public ShapeT<T, sDim, tDim> {
public:
  const std::string Name() const override { return "ArgyrisTri"; }

  T operator()(const PointT<T, sDim, tDim> &z, int i) const override;

  VectorFieldT<T, sDim> LocalGradient(const PointT<T, sDim, tDim> &z, int i) const override {
    return VectorFieldT<T, sDim>(D_x(z, i), D_y(z, i));
  }

  SymTensorT<T, sDim> LocalHessian(const PointT<T, sDim, tDim> &z, int i) const override {
    return SymTensorT<T, sDim>(D_xx(z, i), D_xy(z, i), D_yy(z, i));
  }

  void NodalPoints(const Cell &c, vector<PointT<T, sDim, tDim>> &z) const override {
    z.resize(6);
    for (int i = 0; i < c.Corners(); ++i)
      z[i] = c[i];
    for (int i = 0; i < c.Faces(); ++i)
      z[c.Corners() + i] = c.Face(i);
  }

  void fillValues(const vector<PointT<T, sDim, tDim>> &localQuadraturePoints) override {
    this->localValues.resize(localQuadraturePoints.size(), this->numNodalPoints);
    this->localGradient.resize(localQuadraturePoints.size(), this->numNodalPoints);
    this->localHessian.resize(localQuadraturePoints.size(), this->numNodalPoints);
    for (int q = 0; q < localQuadraturePoints.size(); ++q)
      for (int i = 0; i < this->numNodalPoints; ++i) {
        this->localValues[q][i] = (*this)(localQuadraturePoints[q], i);
        this->localGradient[q][i] = this->LocalGradient(localQuadraturePoints[q], i);
        this->localHessian[q][i] = this->LocalHessian(localQuadraturePoints[q], i);
      }
  }

  void
  fillFaceValues(const vector<vector<PointT<T, sDim, tDim>>> &localFaceQuadraturePoints) override {
    this->localFaceValues.resize(localFaceQuadraturePoints.size(),
                                 localFaceQuadraturePoints[0].size(), this->numNodalPoints);
    for (int face = 0; face < localFaceQuadraturePoints.size(); ++face)
      for (int q = 0; q < localFaceQuadraturePoints[face].size(); ++q)
        for (int i = 0; i < this->numNodalPoints; ++i) {
          this->localFaceValues[face][q][i] = (*this)(localFaceQuadraturePoints[face][q], i);
        }
  }

  explicit ArgyrisTriT() : ShapeT<T, sDim, tDim>(21) {}
private:
  T D_x(const PointT<T, sDim, tDim> &z, int i) const;

  T D_y(const PointT<T, sDim, tDim> &z, int i) const;

  T D_xx(const PointT<T, sDim, tDim> &z, int i) const;

  T D_xy(const PointT<T, sDim, tDim> &z, int i) const;

  T D_yy(const PointT<T, sDim, tDim> &z, int i) const;
};

typedef ArgyrisTriT<> ArgyrisTri;

template<typename T, int sDim, int tDim>
T ArgyrisTriT<T, sDim, tDim>::operator()(const PointT<T, sDim, tDim> &z, int i) const {
  switch (i) {
  case 0:
    return 1 - 10 * (pow(z[0], 3) + pow(z[1], 3))
           + 15 * (z[0] * z[0] - z[1] * z[1]) * (z[0] * z[0] - z[1] * z[1])
           - 6 * (pow(z[0], 5) + pow(z[1], 5)) + 30 * z[0] * z[0] * z[1] * z[1] * (z[0] + z[1]);
  case 1:
    return ((((-3 * z[0] + 8) * z[0] + (z[1] * z[1] - 6)) * z[0] - 10 * z[1] * z[1] * (z[1] - 1))
                * z[0]
            - (4 * z[1] + 1) * (2 * z[1] - 1) * (z[1] - 1) * (z[1] - 1))
           * z[0];
  case 2:
    return ((((-3 * z[1] + 8) * z[1] + (z[0] * z[0] - 6)) * z[1] - 10 * z[0] * z[0] * (z[0] - 1))
                * z[1]
            - (4 * z[0] + 1) * (2 * z[0] - 1) * (z[0] - 1) * (z[0] - 1))
           * z[1];
  case 3:
    return z[0] * z[0] / 2
           * (1 - pow(z[0], 3) + 2 * pow(z[1], 3) - 3 * (1 - z[0]) * (z[0] + z[1] * z[1]));
  case 4:
    return z[0] * z[1] * (1 + (z[0] + z[1]) * (-4 + (z[0] + z[1]) * (5 - 2 * (z[0] + z[1]))));
  case 5:
    return z[1] * z[1] / 2
           * (1 - pow(z[1], 3) + 2 * pow(z[0], 3) - 3 * (1 - z[1]) * (z[1] + z[0] * z[0]));
  case 6:
    return z[0] * z[0]
           * (10 * z[0] + 6 * pow(z[0], 3)
              - 15 * (z[1] * z[1] * (z[1] - 1) + z[0] * (z[0] + z[1] * z[1])));
  case 7:
    return z[0] * z[0]
           * (z[0] * (-4 + z[0] * (7 - 3 * z[0])) - T(7.0) / 2 * z[1] * z[1] * (1 - z[0] - z[1]));
  case 8:
    return (-5 + 2 * z[0] * (7 - 4 * z[0]) + 18.5 * z[1] * (1 - z[0]) - 13.5 * z[1] * z[1]) * z[0]
           * z[0] * z[1];
  case 9:
    return 0.25 * z[0] * z[0]
           * (2 * z[0] * (1 - z[0]) * (1 - z[0]) + z[1] * z[1] * (1 - z[0] - z[1]));
  case 10:
    return (2 + 2 * z[0] * (-3 + 2 * z[0]) - 7 * z[1] * (1 - z[0]) + 5 * z[1] * z[1]) * 0.5 * z[0]
           * z[0] * z[1];
  case 11:
    return 0.25 * z[0] * z[0] * z[1] * z[1] * (5 - 3 * z[0] - 5 * z[1]);
  case 12:
    return z[1] * z[1]
           * (10 * z[1] + 6 * pow(z[1], 3)
              - 15 * (z[0] * z[0] * (z[0] - 1) + z[1] * (z[1] + z[0] * z[0])));
  case 13:
    return (-5 + 2 * z[1] * (7 - 4 * z[1]) + 18.5 * z[0] * (1 - z[1]) - 13.5 * z[0] * z[0]) * z[0]
           * z[1] * z[1];
  case 14:
    return z[1] * z[1]
           * (z[1] * (-4 + z[1] * (7 - 3 * z[1])) - 3.5 * z[0] * z[0] * (1 - z[0] - z[1]));
  case 15:
    return 0.25 * z[0] * z[0] * z[1] * z[1] * (5 - 3 * z[1] - 5 * z[0]);
  case 16:
    return (2 + 2 * z[1] * (-3 + 2 * z[1]) - 7 * z[0] * (1 - z[1]) + 5 * z[0] * z[0]) * 0.5 * z[0]
           * z[1] * z[1];
  case 17:
    return 0.25 * z[1] * z[1]
           * (2 * z[1] * (1 - z[1]) * (1 - z[1]) + z[0] * z[0] * (1 - z[0] - z[1]));
  case 18:
    return -16 * z[0] * z[0] * z[1] * (1 - z[0] - z[1]) * (1 - z[0] - z[1]);
  case 19:
    return -8 * sqrt(2) * z[0] * z[0] * z[1] * z[1] * (1 - z[0] - z[1]);
  case 20:
    return -16 * z[0] * z[1] * z[1] * (1 - z[0] - z[1]) * (1 - z[0] - z[1]);
  default:
    THROW("Not implemented")
  }
}

template<typename T, int sDim, int tDim>
T ArgyrisTriT<T, sDim, tDim>::D_x(const PointT<T, sDim, tDim> &z, int i) const {
  switch (i) {
  case 0:
    return (z[0] * (z[0] * z[0] - 2 * z[0] - (3 * z[1] * z[1] - 1)) - z[1] * z[1] * (2 * z[1] - 2))
           * (-30) * z[0];
  case 1:
    return (((-15 * z[0] + 32) * z[0] + 3 * (z[1] * z[1] - 6)) * z[0]
            + 20 * (1 - z[1]) * z[1] * z[1])
               * z[0]
           + (-8 * z[1] * z[1] + 18 * z[1] - 11) * z[1] * z[1] + 1;
  case 2:
    return 2 * ((-16 * z[0] + 3 * (-5 * z[1] + 9)) * z[0] + ((z[1] + 10) * z[1] - 11)) * z[1]
           * z[0];
  case 3:
    return (((-2.5 * z[0] + 6) * z[0] + 4.5 * (-1 + z[1] * z[1])) * z[0]
            + ((2 * z[1] - 3) * z[1] * z[1] + 1))
           * z[0];
  case 4:
    return ((-8 * z[0] - 18 * z[1] + 15) * z[0] * z[1] + 4 * ((-3 * z[1] + 5) * z[1] - 2) * z[1])
               * z[0]
           + (((-2 * z[1] + 5) * z[1] - 4) * z[1] + 1) * z[1];
  case 5:
    return -3 * (1 - z[0] - z[1]) * z[0] * z[1] * z[1];
  case 6:
    return 15
           * ((2 * (z[0] - 2) * z[0] + (-3 * z[1] * z[1] + 2)) * z[0]
              + 2 * (1 - z[1]) * z[1] * z[1])
           * z[0];
  case 7:
    return (((-15 * z[0] + 28) * z[0] + (-12 + 10.5 * z[1] * z[1])) * z[0]
            - 7 * (1 - z[1]) * z[1] * z[1])
           * z[0];
  case 8:
    return ((-32 * z[0] + (42 - 55.5 * z[1])) * z[0] + ((-27 * z[1] + 37) * z[1] - 10)) * z[0]
           * z[1];
  case 9:
    return (((2.5 * z[0] - 4) * z[0] + (1.5 - 0.75 * z[1] * z[1])) * z[0]
            + 0.5 * (1 - z[1]) * z[1] * z[1])
           * z[0];
  case 10:
    return ((8 * z[0] - 9 + 10.5 * z[1]) * z[0] + ((5 * z[1] - 7) * z[1] + 2)) * z[0] * z[1];
  case 11:
    return (-2.25 * z[0] + 2.5 * (1 - z[1])) * z[0] * z[1] * z[1];
  case 12:
    return 15 * (-3 * z[0] + 2 * (1 - z[1])) * z[0] * z[1] * z[1];
  case 13:
    return ((-40.5 * z[0] + 37 * (1 - z[1])) * z[0] + 2 * (-4 * z[1] + 7) * z[1] - 5) * z[1] * z[1];
  case 14:
    return (10.5 * z[0] - 7 * (1 - z[1])) * z[0] * z[1] * z[1];
  case 15:
    return (2.5 - 3.75 * z[0] - 1.5 * z[1]) * z[0] * z[1] * z[1];
  case 16:
    return ((7.5 * z[0] - 7 * (1 - z[1])) * z[0] + ((2 * z[1] - 3) * z[1] + 1)) * z[1] * z[1];
  case 17:
    return (-0.75 * z[0] + 0.5 * (1 - z[1])) * z[0] * z[1] * z[1];
  case 18:
    return 32 * (-2 * z[0] * z[0] + (3 * z[0] - 1 + z[1]) * (1 - z[1])) * z[0] * z[1];
  case 19:
    return 8 * sqrt(T(2.0)) * (3 * z[0] - 2 * (1 - z[1])) * z[0] * z[1] * z[1];
  case 20:
    return 16 * (-3 * z[0] * z[0] + (4 * z[0] - 1 + z[1]) * (1 - z[1])) * z[1] * z[1];
  default:
    THROW("Not implemented")
  }
}

template<typename T, int sDim, int tDim>
T ArgyrisTriT<T, sDim, tDim>::D_y(const PointT<T, sDim, tDim> &z, int i) const {
  switch (i) {
  case 0:
    return (z[1] * (z[1] * z[1] - 2 * z[1] - (3 * z[0] * z[0] - 1)) - z[0] * z[0] * (2 * z[0] - 2))
           * (-30) * z[1];
  case 1:
    return 2 * ((-16 * z[1] + 3 * (-5 * z[0] + 9)) * z[1] + ((z[0] + 10) * z[0] - 11)) * z[0]
           * z[1];
  case 2:
    return (((-15 * z[1] + 32) * z[1] + 3 * (z[0] * z[0] - 6)) * z[1]
            + 20 * (1 - z[0]) * z[0] * z[0])
               * z[1]
           + (-8 * z[0] * z[0] + 18 * z[0] - 11) * z[0] * z[0] + 1;
  case 3:
    return -3 * (1 - z[0] - z[1]) * z[0] * z[0] * z[1];
  case 4:
    return ((-8 * z[1] - 18 * z[0] + 15) * z[0] * z[1] + 4 * ((-3 * z[0] + 5) * z[0] - 2) * z[0])
               * z[1]
           + (((-2 * z[0] + 5) * z[0] - 4) * z[0] + 1) * z[0];
  case 5:
    return (((-2.5 * z[1] + 6) * z[1] + 4.5 * (-1 + z[0] * z[0])) * z[1]
            + ((2 * z[0] - 3) * z[0] * z[0] + 1))
           * z[1];
  case 6:
    return 15 * (-3 * z[1] + 2 * (1 - z[0])) * z[0] * z[0] * z[1];
  case 7:
    return (10.5 * z[1] - 7 * (1 - z[0])) * z[0] * z[0] * z[1];
  case 8:
    return ((-40.5 * z[1] + 37 * (1 - z[0])) * z[1] + 2 * (-4 * z[0] + 7) * z[0] - 5) * z[0] * z[0];
  case 9:
    return (-0.75 * z[1] + 0.5 * (1 - z[0])) * z[0] * z[0] * z[1];
  case 10:
    return ((7.5 * z[1] - 7 * (1 - z[0])) * z[1] + ((2 * z[0] - 3) * z[0] + 1)) * z[0] * z[0];
  case 11:
    return (2.5 - 3.75 * z[1] - 1.5 * z[0]) * z[0] * z[0] * z[1];
  case 12:
    return 15
           * ((2 * (z[1] - 2) * z[1] + (-3 * z[0] * z[0] + 2)) * z[1]
              + 2 * (1 - z[0]) * z[0] * z[0])
           * z[1];
  case 13:
    return ((-32 * z[1] + (42 - 55.5 * z[0])) * z[1] + ((-27 * z[0] + 37) * z[0] - 10)) * z[0]
           * z[1];
  case 14:
    return (((-15 * z[1] + 28) * z[1] + (-12 + 10.5 * z[0] * z[0])) * z[1]
            - 7 * (1 - z[0]) * z[0] * z[0])
           * z[1];
  case 15:
    return (-2.25 * z[1] + 2.5 * (1 - z[0])) * z[0] * z[0] * z[1];
  case 16:
    return ((8 * z[1] - 9 + 10.5 * z[0]) * z[1] + ((5 * z[0] - 7) * z[0] + 2)) * z[0] * z[1];
  case 17:
    return (((2.5 * z[1] - 4) * z[1] + (1.5 - 0.75 * z[0] * z[0])) * z[1]
            + 0.5 * (1 - z[0]) * z[0] * z[0])
           * z[1];
  case 18:
    return 16 * (-3 * z[1] * z[1] + (4 * z[1] - 1 + z[0]) * (1 - z[0])) * z[0] * z[0];
  case 19:
    return 8 * sqrt(T(2.0)) * (3 * z[1] - 2 * (1 - z[0])) * z[0] * z[0] * z[1];
  case 20:
    return 32 * (-2 * z[1] * z[1] + (3 * z[1] - 1 + z[0]) * (1 - z[0])) * z[0] * z[1];
  default:
    THROW("Not implemented")
  }
}

template<typename T, int sDim, int tDim>
T ArgyrisTriT<T, sDim, tDim>::D_xx(const PointT<T, sDim, tDim> &z, int i) const {
  switch (i) {
  case 0:
    return 60
           * (z[0] * (z[0] * (-2 * z[0] + 3) + (3 * z[1] * z[1] - 1)) + z[1] * z[1] * (z[1] - 1));
  case 1:
    return (12 * (-5 * z[0] + 8) * z[0] + 6 * (z[1] * z[1] - 6)) * z[0]
           + 20 * (1 - z[1]) * z[1] * z[1];
  case 2:
    return (6 * (-8 * z[0] + (-5 * z[1] + 9)) * z[0] + (z[1] + 10) * z[1] - 11) * 2 * z[1];
  case 3:
    return (2 * (-5 * z[0] + 9) * z[0] + 9 * (z[1] * z[1] - 1)) * z[0]
           + (2 * z[1] - 3) * z[1] * z[1] + 1;
  case 4:
    return 6 * (-4 * z[0] - 6 * z[1] + 5) * z[1] * z[0] + (4 * (-3 * z[1] + 5) * z[1] - 8) * z[1];
  case 5:
    return 3 * (z[1] + 2 * z[0] - 1) * z[1] * z[1];
  case 6:
    return 30
           * ((2 * (2 * z[0] - 3) * z[0] + (-3 * z[1] * z[1] + 2)) * z[0]
              + (1 - z[1]) * z[1] * z[1]);
  case 7:
    return (12 * (-5 * z[0] + 7) * z[0] + 3 * (7 * z[1] * z[1] - 8)) * z[0]
           - 7 * (1 - z[1]) * z[1] * z[1];
  case 8:
    return ((-96 * z[0] + (-111 * z[1] + 84)) * z[0] + ((-27 * z[1] + 37) * z[1] - 10)) * z[1];
  case 9:
    return (2 * (5 * z[0] - 6) * z[0] + (3 - 1.5 * z[1] * z[1])) * z[0]
           + 0.5 * (1 - z[1]) * z[1] * z[1];
  case 10:
    return (3 * (8 * z[0] + 7 * z[1] - 6) * z[0] + ((5 * z[1] - 7) * z[1] + 2)) * z[1];
  case 11:
    return (2.5 * (1 - z[1]) - 4.5 * z[0]) * z[1] * z[1];
  case 12:
    return 30 * (1 - z[1] - 3 * z[0]) * z[1] * z[1];
  case 13:
    return (37 * (1 - z[1]) - 81 * z[0]) * z[1] * z[1];
  case 14:
    return 7 * (z[1] + 3 * z[0] - 1) * z[1] * z[1];
  case 15:
    return (2.5 - 1.5 * z[1] - 7.5 * z[0]) * z[1] * z[1];
  case 16:
    return (15 * z[0] - 7 * (1 - z[1])) * z[1] * z[1];
  case 17:
    return (0.5 * (1 - z[1]) - 1.5 * z[0]) * z[1] * z[1];
  case 18:
    return 32 * (-6 * z[0] * z[0] + (6 * z[0] - 1 + z[1]) * (1 - z[1])) * z[1];
  case 19:
    return -16 * sqrt(T(2.0)) * (1 - 3 * z[0] - z[1]) * z[1] * z[1];
  case 20:
    return 32 * (2 * (1 - z[1]) - 3 * z[0]) * z[1] * z[1];
  default:
    THROW("Not implemented")
  }
}

template<typename T, int sDim, int tDim>
T ArgyrisTriT<T, sDim, tDim>::D_xy(const PointT<T, sDim, tDim> &z, int i) const {
  switch (i) {
  case 0:
    return 60 * z[0] * z[1] * (3 * (z[0] + z[1]) - 2);
  case 1:
    return ((3 * z[0] + 20) * z[0] + (-30 * z[0] - 16 * z[1] + 27) * z[1] - 11) * 2 * z[1];
  case 2:
    return ((3 * z[1] + 20) * z[1] + (-30 * z[1] - 16 * z[0] + 27) * z[0] - 11) * 2 * z[0];
  case 3:
    return 3 * (3 * z[0] + 2 * z[1] - 2) * z[0] * z[1];
  case 4:
    return ((-8 * z[0] - 36 * z[1] + 15) * z[0] + 4 * (-9 * z[1] + 10) * z[1] - 8) * z[0]
           + ((-8 * z[1] + 15) * z[1] - 8) * z[1] + 1;
  case 5:
    return 3 * (3 * z[1] + 2 * z[0] - 2) * z[0] * z[1];
  case 6:
    return 30 * (2 - 3 * z[0] - 3 * z[1]) * z[0] * z[1];
  case 7:
    return 7 * (3 * (z[0] + z[1]) - 2) * z[0] * z[1];
  case 8:
    return ((-32 * z[0] + (-111 * z[1] + 42)) * z[0] + ((-81 * z[1] + 74) * z[1] - 10)) * z[0];
  case 9:
    return (1 - 1.5 * (z[0] + z[1])) * z[0] * z[1];
  case 10:
    return ((8 * z[0] + 3 * (7 * z[1] - 3)) * z[0] + ((15 * z[1] - 14) * z[1] + 2)) * z[0];
  case 11:
    return (-4.5 * z[0] + 5 * (1 - 1.5 * z[1])) * z[0] * z[1];
  case 12:
    return 30 * (2 - 3 * z[1] - 3 * z[0]) * z[0] * z[1];
  case 13:
    return ((-32 * z[1] + (-111 * z[0] + 42)) * z[1] + ((-81 * z[0] + 74) * z[0] - 10)) * z[1];
  case 14:
    return 7 * (3 * (z[0] + z[1]) - 2) * z[0] * z[1];
  case 15:
    return (-4.5 * z[1] + 5 * (1 - 1.5 * z[0])) * z[0] * z[1];
  case 16:
    return ((8 * z[1] + 3 * (7 * z[0] - 3)) * z[1] + ((15 * z[0] - 14) * z[0] + 2)) * z[1];
  case 17:
    return (1 - 1.5 * (z[0] + z[1])) * z[0] * z[1];
  case 18:
    return 32 * ((-2 * z[0] - 6 * z[1] + 3) * z[0] + ((4 - 3 * z[1]) * z[1] - 1)) * z[0];
  case 19:
    return -16 * sqrt(T(2.0)) * (2 - 3 * (z[0] + z[1])) * z[0] * z[1];
  case 20:
    return 32 * ((-2 * z[1] - 6 * z[0] + 3) * z[1] + ((4 - 3 * z[0]) * z[0] - 1)) * z[1];
  default:
    THROW("Not implemented")
  }
}

template<typename T, int sDim, int tDim>
T ArgyrisTriT<T, sDim, tDim>::D_yy(const PointT<T, sDim, tDim> &z, int i) const {
  switch (i) {
  case 0:
    return 60
           * (z[1] * (z[1] * (-2 * z[1] + 3) + (3 * z[0] * z[0] - 1)) + z[0] * z[0] * (z[0] - 1));
  case 1:
    return (6 * (-8 * z[1] + (-5 * z[0] + 9)) * z[1] + (z[0] + 10) * z[0] - 11) * 2 * z[0];
  case 2:
    return (12 * (-5 * z[1] + 8) * z[1] + 6 * (z[0] * z[0] - 6)) * z[1]
           + 20 * (1 - z[0]) * z[0] * z[0];
  case 3:
    return 3 * (z[0] + 2 * z[1] - 1) * z[0] * z[0];
  case 4:
    return 6 * (-4 * z[1] - 6 * z[0] + 5) * z[0] * z[1] + (4 * (-3 * z[0] + 5) * z[0] - 8) * z[0];
  case 5:
    return (2 * (-5 * z[1] + 9) * z[1] + 9 * (z[0] * z[0] - 1)) * z[1]
           + (2 * z[0] - 3) * z[0] * z[0] + 1;
  case 6:
    return 30 * (1 - z[0] - 3 * z[1]) * z[0] * z[0];
  case 7:
    return 7 * (z[0] + 3 * z[1] - 1) * z[0] * z[0];
  case 8:
    return (37 * (1 - z[0]) - 81 * z[1]) * z[0] * z[0];
  case 9:
    return (0.5 * (1 - z[0]) - 1.5 * z[1]) * z[0] * z[0];
  case 10:
    return (15 * z[1] - 7 * (1 - z[0])) * z[0] * z[0];
  case 11:
    return (2.5 - 1.5 * z[0] - 7.5 * z[1]) * z[0] * z[0];
  case 12:
    return 30
           * ((2 * (2 * z[1] - 3) * z[1] + (-3 * z[0] * z[0] + 2)) * z[1]
              + (1 - z[0]) * z[0] * z[0]);
  case 13:
    return ((-96 * z[1] + (-111 * z[0] + 84)) * z[1] + ((-27 * z[0] + 37) * z[0] - 10)) * z[0];
  case 14:
    return (12 * (-5 * z[1] + 7) * z[1] + 3 * (7 * z[0] * z[0] - 8)) * z[1]
           - 7 * (1 - z[0]) * z[0] * z[0];
  case 15:
    return (2.5 * (1 - z[0]) - 4.5 * z[1]) * z[0] * z[0];
  case 16:
    return (3 * (8 * z[1] + 7 * z[0] - 6) * z[1] + ((5 * z[0] - 7) * z[0] + 2)) * z[0];
  case 17:
    return (2 * (5 * z[1] - 6) * z[1] + (3 - 1.5 * z[0] * z[0])) * z[1]
           + 0.5 * (1 - z[0]) * z[0] * z[0];
  case 18:
    return 32 * (2 * (1 - z[0]) - 3 * z[1]) * z[0] * z[0];
  case 19:
    return -16 * sqrt(T(2.0)) * (1 - 3 * z[1] - z[0]) * z[0] * z[0];
  case 20:
    return 32 * (-6 * z[1] * z[1] + (6 * z[1] - 1 + z[0]) * (1 - z[0])) * z[0];
  default:
    THROW("Not implemented")
  }
}

#endif // ARGYRISSHAPES_HPP
