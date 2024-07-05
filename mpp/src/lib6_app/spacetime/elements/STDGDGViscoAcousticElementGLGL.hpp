#ifndef SPACETIME_STDGDGVISCOACOUSTICELEMENTGLGL_HPP
#define SPACETIME_STDGDGVISCOACOUSTICELEMENTGLGL_HPP

#include "VectorMatrixBase.hpp"
#include "Vector.hpp"
#include "STDGDGViscoAcousticElement.hpp"

class STDGDGViscoAcousticElementGLGL {
public:

  double det_space;
  double det_time;

  std::vector<double> space_qWeight;
  std::vector<double> time_qWeight;
  std::vector<double> qWeight;
  std::vector<Point> qPoint;

  const int dim;
  const SpaceTimeShape &shape;
  const SpaceTimeQuadrature &quad;
  const Quadrature &space_quad;
  const cell &TC;
  const row r;

  ShapeValues<VectorFieldT<double, SpaceDimension, 1>> gradients;
  std::vector<VectorField> velocityST;
  std::vector<VectorField> dtVelocityST;
  std::vector<double> divVelocityST;

  std::vector<double> pressureST;
  std::vector<double> dtPressureST;
  std::vector<VectorField> gradPressureST;

  int numL;

  std::vector<std::vector<int>> component_to_index{3, std::vector<int>{}};

  STDGDGViscoAcousticElementGLGL(const VectorMatrixBase &,
                                 const cell &,
                                 int nL = 0,
                                 bool max_quad_order = false);

  int nQ() const;

  int nSpaceQ() const;

  int nTimeQ() const;

  double QWeight(int q) const;

  double SpaceQWeight(int qq);

  Point SpaceQPoint(int q) const;

  double TimeQWeight(int qq);

  const Point &QPoint(int q) const;

  double Area() const;

  int i_dimension() const;

  int j_dimension() const;

  int variable(int j) const;

  int GetComponent(int j) const;

  std::vector<int> GetIndicesForComponent(int component) const{
    return component_to_index[component];
  }

  double EvaluateComponentLocal(const Point &p, const Vector &u, int component) const;

  double EvaluateComponentGlobal(const Point &p, const Vector &u, int component) const;

  double EvaluateComponent(int nq, const Vector &u, int component) const;

  VectorField Velocity(int nq, int j);

  VectorField Velocity(int nq, const Vector &U);

  VectorField VelocityLocal(const Point &localPoint, int j);

  VectorField VelocityLocal(const Point &localPoint, const Vector &u);

  VectorField VelocityGlobal(const Point &globalPoint, const Vector &u);

  VectorField DtVelocity(int nq, int j);

  VectorField DtVelocity(int nq, const Vector &u);

  double DivVelocity(int nq, int j);

  double DivVelocity(int nq, const Vector &U);

  double Pressure(int nq, int j);

  double Pressure(int nq, const Vector &U);

  double PressureLocal(const Point &localPoint, int j);

  double PressureLocal(const Point &localPoint, const Vector &u);

  double PressureGlobal(const Point &globalPoint, int i);

  double PressureGlobal(const Point &globalPoint, const Vector &u);

  double DtPressure(int nq, int j);

  double DtPressure(int nq, const Vector &u);

  VectorField GradPressure(int nq, int j);

  VectorField GradPressure(int nq, const Vector &u);
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


#endif //SPACETIME_STDGDGVISCOACOUSTICELEMENTGLGL_HPP
