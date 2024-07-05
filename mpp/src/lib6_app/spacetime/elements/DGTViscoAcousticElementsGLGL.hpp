#ifndef SPACETIME_DGTVISCOACOUSTICELEMENTSGLGL_HPP
#define SPACETIME_DGTVISCOACOUSTICELEMENTSGLGL_HPP

#include "VectorMatrixBase.hpp"


class SpaceTimeViscoAcousticDGTElementGLGL {
public:
  const DegreePair deg;
  const Quadrature &SpaceQ;
  const Quadrature &TimeQ;
  vector<double> qWeight;
  vector<double> space_qWeight;
  vector<double> time_qWeight;
  vector<Point> qPoint;

  //vector<Transformation> T;
  int dim;
  const Shape &SS;
  const Shape &TS;
  const ShapeValues<Scalar> &spaceValue;
  const ShapeValues<Scalar> &timeValue;
  const cell &TC;
  double timedet;

  vector<vector<VectorField> > time_gradient;
  //vector<vector<VectorField> > testspace_time_gradient;

  vector<VectorField> velocityST;
  vector<VectorField> dtVelocityST;
  vector<double> divVelocityST;

  vector<double> pressureST;
  vector<double> dtPressureST;
  vector<VectorField> gradPressureST;

  int numL;

  SpaceTimeViscoAcousticDGTElementGLGL(const STDiscretization &,
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
    for (int q = 0; q < nQ(); ++q) a += qWeight[q];
    return a;
  }

  int dimension() const {Exit("Not implemented"); }

  int i_dimension() const { return SS.size() * TS.size() * (1 + dim + numL); }

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
    return velocityST[nq * i_dimension() + i];
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
    return pressureST[nq * i_dimension() + i];
  }

  double Evaluate(const Point &localPT, int i) const {
    // return TestSpace(localPT.t(),i/(SS.size()*(dim+1+numL)))*SS(localPT,i%SS.size());
    int j = i % SS.size();
    int k = i / (SS.size() * (dim + 1 + numL));
    return SS(localPT, j) * TS(localPT.t(), k);
  }
};

class SpaceTimeViscoAcousticDGTFaceElementGLGL {
  const DegreePair deg;
  const Quadrature &SpaceQ;
  const Quadrature &TimeQ;
  vector<double> qWeight;
  vector<Point> qLocal;
  vector<Point> qPoint;
  vector<Point> qNormal;
  vector<Transformation> T;
  vector<bool> i_zero;
  vector<bool> j_zero;

  int dim;
  const Shape &SS;
  const Shape &TS;

  //const cell &TC;
  //double timedet;
  //face f;
  int fid;
  double ww;

  vector<vector<int>> non_zero_indices_for_q;

  vector<vector<Point>> faceqpoints;

  //ViscoAcoustic
  vector<double> pressureST;
  vector<VectorField> velocityST;
  int numL;
public:
  SpaceTimeViscoAcousticDGTFaceElementGLGL(const STDiscretization &,
                                           const VectorMatrixBase &,
                                           const cell &,
                                           int f_id,
                                           int nL);

  SpaceTimeViscoAcousticDGTFaceElementGLGL(const STDiscretization &,
                                           const VectorMatrixBase &,
                                           const cell &,
                                           int f_id,
                                           int nL,
                                           string dgdg);

  SpaceTimeViscoAcousticDGTFaceElementGLGL(const STDiscretization &,
                                           const VectorMatrixBase &,
                                           const cell &,
                                           int f_id,
                                           int nL,
                                           string dgdg,
                                           const cell &);

  void updateCell(cell tc);

  const Point &QLocal(int q) const { return qLocal[q]; }

  const Point &QPoint(int q) const { return qPoint[q]; }

  const Point &QNormal(int q) const { return qNormal[q]; }

  double QWeight(int q) const { return qWeight[q]; }

  int nSpaceQ() const { return SpaceQ.size(); }

  int nTimeQ() const { return TimeQ.size(); }

  int nQ() const { return SpaceQ.size() * TimeQ.size(); }

  bool is_zero_test(int nq, int i) const { return i_zero[nq * i_dimension() + i]; }

  bool is_zero(int nq, int j) const { return j_zero[nq * j_dimension() + j]; }

  int i_dimension() const { return SS.size() * TS.size() * (1 + dim + numL); }

  int j_dimension() const { return SS.size() * TS.size() * (1 + dim + numL); }

  std::vector<int> get_non_zero_indices(int q) const {
    return non_zero_indices_for_q[q];
  }

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

  VectorField VelocityTestSpace(int nq, int i) {
    return velocityST[nq * i_dimension() + i];
  }

  double Pressure(int nq, int j) {
    return pressureST[nq * j_dimension() + j];
  }

  double PressureTestSpace(int nq, int i) {
    return pressureST[nq * i_dimension() + i];
  }

  VectorField Velocity(Point QP, int i) {
    int var = variable(i);
    if (var != 0) return zero;

    int ti = i / (SS.size() * (1 + dim + numL));
    int si = i % (SS.size() * (1 + dim + numL));
    int ssi = si / SS.size();
    int ssj = si % SS.size();

    VectorField V = zero;
    V[ssi] = 1.0;
    V *= SS(QP, ssj) * TS(Point(QP.t()), ti);
    return V;
  }

  double Pressure(Point QP, int i) {
    int var = variable(i);
    if (var == 0) return 0.0;

    int ti = i / (SS.size() * (1 + dim + numL));
    int si = i % SS.size();

    Scalar P = SS(QP, si) * TS(Point(QP.t()), ti);
    return P;
  }
};


#endif //SPACETIME_DGTVISCOACOUSTICELEMENTSGLGL_HPP
