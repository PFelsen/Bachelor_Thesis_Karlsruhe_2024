#ifndef SPACETIMEVISCOACOUSTICFACEELEMENT_HPP
#define SPACETIMEVISCOACOUSTICFACEELEMENT_HPP

#include "SpaceTimeDiscretization.hpp"
#include "VectorMatrixBase.hpp"
#include "Vector.hpp"

class SpaceTimeViscoAcousticFaceElement {
public:
  const DegreePair deg;
private:
  const Quadrature &SpaceQ;
  const Quadrature &TimeQ;
  vector<double> qWeight;
  vector<Point> qLocal;
  vector<Point> qPoint;
  vector<Point> qNormal;
  Transformation T;
  vector<bool> i_zero;
  vector<bool> j_zero;

  int dim;
  const Shape &SS;
  const Shape &TS;
  const Shape &TestSpace;

  const cell &TC;
  face f;
  int fid;
  const vector<vector<Point>> &faceqpoints;
  double ww;

  const row r;

  vector<double> pressureST;
  vector<double> pressureST_TestSpace;
  vector<VectorField> velocityST;
  vector<VectorField> velocityST_TestSpace;
  int numL;
public:
  SpaceTimeViscoAcousticFaceElement(const STDiscretization &, const VectorMatrixBase &,
                                    const cell &, int f_id, int nL);

  const Point &LocalFaceCorner(int i) const {
    return TC.LocalCorner(TC.facecorner(fid, i));
  }

  const Point &QLocal(int q) const { return qLocal[q]; }

  const Point &QPoint(int q) const { return qPoint[q]; }

  const Point &QNormal(int q) const { return qNormal[q]; }

  double QWeight(int q) const { return qWeight[q]; }

  int nSpaceQ() const { return SpaceQ.size(); }

  int nTimeQ() const { return TimeQ.size(); }

  int nQ() const { return SpaceQ.size() * TimeQ.size(); }

  bool is_zero_test(int nq, int i) const { return i_zero[nq * i_dimension() + i]; }

  bool is_zero(int nq, int j) const { return j_zero[nq * j_dimension() + j]; }

  int i_dimension() const { return SS.size() * TestSpace.size() * (1 + dim + numL); }

  int j_dimension() const { return SS.size() * TS.size() * (1 + dim + numL); }

  VectorField Velocity(int nq, int j) {
    return velocityST[nq * j_dimension() + j];
  }

  VectorField Velocity(int nq, const Vector &U) {
    VectorField V = zero;
    for (int j = 0; j < j_dimension(); ++j) {
      V += U(r, j) * Velocity(nq, j);
    }
    return V;
  }

  VectorField VelocityTestSpace(int nq, int i) {
    return velocityST_TestSpace[nq * i_dimension() + i];
  }

  double Pressure(int nq, int j) {
    return pressureST[nq * j_dimension() + j];
  }

  double Pressure(int nq, const Vector &U) {
    double P = 0.0;
    for (int j = 0; j < j_dimension(); ++j) {
      P += U(r, j) * Pressure(nq, j);
    }
    return P;
  }
  
  double PressureTestSpace(int nq, int i) {
    return pressureST_TestSpace[nq * i_dimension() + i];
  }
};

#endif //SPACETIMEVISCOACOUSTICFACEELEMENT_HPP
