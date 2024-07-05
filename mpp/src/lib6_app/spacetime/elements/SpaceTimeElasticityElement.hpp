#ifndef SPACETIMEELASTICITYELEMENT_HPP
#define SPACETIMEELASTICITYELEMENT_HPP

#include "SpaceTimeDiscretization.hpp"
#include "VectorMatrixBase.hpp"


class SpaceTimeElasticityElement {
  const DegreePair deg;
  const Quadrature &SpaceQ;
  const Quadrature &TimeQ;
  vector<double> qWeight;
  vector<double> space_qWeight;
  vector<double> time_qWeight;
  vector<Point> qPoint;

  vector<Transformation> T;
  int dim;
  const Shape &SS;
  const Shape &TS;
  const Shape &TestSpace;
  const ShapeValues<Scalar> &spaceValue;
  const ShapeValues<Scalar> &timeValue;
  const ShapeValues<Scalar> &testspaceValue;
  const cell &TC;
  double timedet;

  vector<vector<VectorField> > time_gradient;
  vector<VectorField> velocityST;
  vector<VectorField> dtvelocityST;
  vector<Tensor> strainST;
  vector<VectorField> velocityST_TestSpace;
  vector<Tensor> stressST;
  vector<Tensor> dtstressST;
  vector<VectorField> divstressST;
  vector<Tensor> stressST_TestSpace;

public:
  SpaceTimeElasticityElement(const STDiscretization &,
                             const VectorMatrixBase &,
                             const cell &);

  int nQ() const { return SpaceQ.size() * TimeQ.size(); }

  int nSpaceQ() const { return SpaceQ.size(); }

  int nTimeQ() const { return TimeQ.size(); }

  double QWeight(int q) const { return qWeight[q]; }

  double SpaceQWeight(int q) const { return space_qWeight[q]; }

  double TimeQWeight(int q) const { return time_qWeight[q]; }

  const Point &QPoint(int q) const { return qPoint[q]; }

  double Area() const {
    double a = 0;
    for (int q = 0; q < nQ(); ++q)
      a += qWeight[q];
    return a;
  }

  int shape_size() const { return SS.size() * TS.size(); }

  int shape_size_dG() const { return SS.size() * (TS.size() - 1); }

  int space_shape_size_dG() const { return SS.size(); }

  int time_shape_size_dG() const { return TS.size() - 1; }

  int dimension() const {Exit("Not implemented"); }

  int mod() const {
    if (dim == 2) return 5;
    if (dim == 3) return 9;
    Exit("Not implemented");
  }

  int i_dimension() const { return SS.size() * TestSpace.size() * mod(); }

  int j_dimension() const { return SS.size() * TS.size() * mod(); }

  VectorField Velocity(int nq, int j) {
    return velocityST[nq * j_dimension() + j];
    /*int sj = j%(mod() * SS.size() );
    if ( sj < SS.size()*(mod()-dim) )
        return zero;
    sj -= SS.size()*(mod()-dim);
    int tj = j/(mod()*SS.size());
    int tq = nq/SpaceQ.size();
    int sq = nq%SpaceQ.size();
    return velocity[sq][sj] * double(timeValue[tq][tj]);*/
  }

  VectorField DtVelocity(int nq, int j) { return dtvelocityST[nq * j_dimension() + j]; }

  VectorField VelocityTestSpace(int nq, int i) {
    return velocityST_TestSpace[nq * i_dimension() + i];
  }

  Tensor Strain(int nq, int j) { return strainST[nq * j_dimension() + j]; }

  Tensor Stress(int nq, int j) { return stressST[nq * j_dimension() + j]; }

  Tensor DtStress(int nq, int j) { return dtstressST[nq * j_dimension() + j]; }

  Tensor StressTestSpace(int nq, int i) {
    return stressST_TestSpace[nq * i_dimension() + i];
  }

  VectorField DivStress(int nq, int j) { return divstressST[nq * j_dimension() + j]; }
};

class SpaceTimeElasticityFaceElement {
  const DegreePair deg;
  vector<double> qWeight;
  vector<Point> qLocal;
  vector<Point> qPoint;
  vector<Point> qNormal;
  vector<Point> qTangent;
  vector<Transformation> T;

  vector<VectorField> velocityST;
  vector<VectorField> velocityST_TestSpace;
  vector<Tensor> stressST;
  vector<Tensor> stressST_TestSpace;

  int d;
  int dim;
  const Shape &SS;
  const Shape &TS;
  const Shape &TestSpace;
  const cell &TC;
  face f;
  int fid;
  const Quadrature &SpaceQ;
  const Quadrature &TimeQ;
  const ShapeValues<Scalar> &timeValue;
  const ShapeValues<Scalar> &testspaceValue;
  double ww;
  double timedet;
public:
  SpaceTimeElasticityFaceElement(const STDiscretization &,
                                 const VectorMatrixBase &,
                                 const cell &,
                                 int);

  SpaceTimeElasticityFaceElement(const SpaceTimeElasticityFaceElement &);

  const Point &LocalFaceCorner(int i) const {
    return TC.LocalCorner(TC.facecorner(fid, i));
  }

  Point Local(const Point &) const;

  const Point &QLocal(int q) const { return qLocal[q]; }

  const Point &QPoint(int q) const { return qPoint[q]; }

  const Point &QNormal(int q) const { return qNormal[q]; }

  const Point &QTangent(int q) const { return qTangent[q]; }

  double QWeight(int q) const { return qWeight[q]; }

  int nSpaceQ() const { return SpaceQ.size(); }

  int nTimeQ() const { return TimeQ.size(); }

  int nQ() const { return SpaceQ.size() * TimeQ.size(); }

  int shape_size() const { return SS.size(); }

  int mod() const {
    if (dim == 2) return 5;
    if (dim == 3) return 9;
    Exit("Not implemented");
  }

  int i_dimension() const { return SS.size() * TestSpace.size() * mod(); }

  int j_dimension() const { return SS.size() * TS.size() * mod(); }

  VectorField Velocity(int nq, int j) {
    return velocityST[nq * j_dimension() + j];
  }

  VectorField VelocityTestSpace(int nq, int i) {
    return velocityST_TestSpace[nq * i_dimension() + i];
  }

  Tensor Stress(int nq, int j) {
    return stressST[nq * j_dimension() + j];
  }

  Tensor StressTestSpace(int nq, int i) {
    return stressST_TestSpace[nq * i_dimension() + i];
  }
};

#endif //SPACETIMEELASTICITYELEMENT_HPP
