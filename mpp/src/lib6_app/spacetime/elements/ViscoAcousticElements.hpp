#ifndef SPACETIME_VISCOACOUSTICELEMENTS_HPP
#define SPACETIME_VISCOACOUSTICELEMENTS_HPP

#include "SpaceTimeDiscretization.hpp"
#include "VectorMatrixBase.hpp"


class SpaceTimeViscoAcousticElement {
public:
  const DegreePair deg;
  const Quadrature &SpaceQ;
  const Quadrature &TimeQ;
  vector<double> qWeight;
  vector<double> space_qWeight;
  vector<double> time_qWeight;
  vector<Point> qPoint;

  Transformation T;
  int dim;
  const Shape &SS;
  const Shape &TS;
  const Shape &TestSpace;
  const ShapeValues<Scalar> &space_value;
  const ShapeValues<Scalar> &time_value;
  const ShapeValues<Scalar> &testspace_value;
  //const cell &TC;
  double timedet;

  vector<vector<VectorField> > time_gradient;

  vector<VectorField> velocityST;
  vector<VectorField> dtVelocityST;
  vector<double> divVelocityST;
  vector<VectorField> velocityST_TestSpace;

  vector<double> pressureST;
  vector<double> dtPressureST;
  vector<VectorField> gradPressureST;
  vector<double> pressureST_TestSpace;

  int numL;

  SpaceTimeViscoAcousticElement(const STDiscretization &,
                                const VectorMatrixBase &, const cell &, int nL);

  int nQ() const { return SpaceQ.size() * TimeQ.size(); }

  int nSpaceQ() const { return SpaceQ.size(); }

  int nTimeQ() const { return TimeQ.size(); }

  double QWeight(int q) const { return qWeight[q]; }

  double SpaceQWeight(int q) const { return space_qWeight[q]; }

  double TimeQWeight(int q) const { return time_qWeight[q]; }

  const Point &QPoint(int q) const { return qPoint[q]; }

  Point LocalToGlobal(const Point &local) const;

  Point GlobalToLocal(const Point &) const;

  double Area() const {
    double a = 0;
    for (int q = 0; q < nQ(); ++q)
      a += qWeight[q];
    return a;
  }

  int dimension() const {Exit("Not implemented"); }

  int i_dimension() const { return SS.size() * TestSpace.size() * (1 + dim + numL); }

  int j_dimension() const { return SS.size() * TS.size() * (1 + dim + numL); }

  int variable(int j) {
    int tmp = j % (SS.size() * (dim + 1 + numL));
    if (tmp < dim * SS.size()) return 0;
    for (int nL = 0; nL <= numL; nL++)
      if (tmp - SS.size() * dim < (nL + 1) * SS.size()) return 1 + nL;
    Exit("BLUPP");
  }

  VectorField Velocity(int nq, int j) {
    return velocityST[nq * j_dimension() + j];
  }

  VectorField DtVelocity(int nq, int j) {
    return dtVelocityST[nq * j_dimension() + j];
  }

  double DivVelocity(int nq, int j) {
    return divVelocityST[nq * j_dimension() + j];
  }

  VectorField VelocityTestSpace(int nq, int i) {
    return velocityST_TestSpace[nq * i_dimension() + i];
  }

  double Pressure(int nq, int j) {
    return pressureST[nq * j_dimension() + j];
  }

  double DtPressure(int nq, int j) {
    return dtPressureST[nq * j_dimension() + j];
  }

  VectorField GradPressure(int nq, int j) {
    return gradPressureST[nq * j_dimension() + j];
  }

  double PressureTestSpace(int nq, int i) {
    return pressureST_TestSpace[nq * i_dimension() + i];
  }

  double Evaluate(const Point &localPT, int i) const {
    int j = i % SS.size();
    int k = i / (SS.size() * (dim + 1 + numL));
    return SS(localPT, j) * TestSpace(Point(localPT.t()), k);
  }
};

#endif //SPACETIME_VISCOACOUSTICELEMENTS_HPP