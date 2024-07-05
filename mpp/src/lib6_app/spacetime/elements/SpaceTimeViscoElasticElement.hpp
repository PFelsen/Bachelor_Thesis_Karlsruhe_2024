#ifndef SPACETIMEVISCOELASTICELEMENT_HPP
#define SPACETIMEVISCOELASTICELEMENT_HPP

#include "SpaceTimeDiscretization.hpp"
#include "VectorMatrixBase.hpp"

class SpaceTimeViscoElasticElement {
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
  const Shape &SS2;
  const Shape &TS;
  const Shape &TestSpace;
  const ShapeValues<Scalar> &spaceValue;
  const ShapeValues<Scalar> &space_value2;
  const ShapeValues<Scalar> &timeValue;
  const ShapeValues<Scalar> &testspaceValue;

  const cell &TC;
  double timedet;

  vector<vector<VectorField> > time_gradient;

  vector<VectorField> velocityST;
  vector<VectorField> dtVelocityST;
  vector<Tensor> strainST;
  vector<VectorField> velocityST_TestSpace;

  vector<Tensor> stressST;
  vector<Tensor> dtStressST;
  vector<VectorField> divStressST;
  vector<Tensor> stressST_TestSpace;
  /*vector<vector<Tensor>> stressST;
  vector<vector<Tensor>> dtStressST;
  vector<vector<VectorField>> divStressST;
  vector<vector<Tensor>> stressST_TestSpace;*/

  //ViscoElastic
  int numL;
public:

  SpaceTimeViscoElasticElement(const STDiscretization &,
                               const VectorMatrixBase &, const cell &, int nL);

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

  int dimension() const {Exit("Not implemented")}

  int td() const { //TensorDimension
    return std::max(1, dim + (dim - 1) + (dim - 2));
  }

  int i_dimension() const {
    return (SS.size() * (dim + td()) + SS2.size() * td() * numL) * TestSpace.size();
  }

  int j_dimension() const {
    return (SS.size() * (dim + td()) + SS2.size() * td() * numL) * TS.size();
  }

  int variable(int j) {
    int tmp = j % (SS.size() * (dim + td()) + SS2.size() * td() * numL);
    if (tmp < dim * SS.size()) return 0;
    else if (tmp - dim * SS.size() < td() * SS.size()) return 1;
    for (int nL = 0; nL < numL; nL++)
      if (tmp - SS.size() * (dim + td()) < (nL + 1) * td() * SS2.size())
        return 2 + nL;
    Exit("BLUPP");
  }

  VectorField Velocity(int nq, int j) {
    return velocityST[nq * j_dimension() + j];
  }

  VectorField DtVelocity(int nq, int j) {
    return dtVelocityST[nq * j_dimension() + j];
  }

  Tensor Strain(int nq, int j) {
    return strainST[nq * j_dimension() + j];
  }

  VectorField VelocityTestSpace(int nq, int i) {
    return velocityST_TestSpace[nq * i_dimension() + i];
  }

  Tensor Stress(int nq, int j) {
    return stressST[nq * j_dimension() + j];
  }

  Tensor DtStress(int nq, int j) {
    return dtStressST[nq * j_dimension() + j];
  }

  Tensor StressTestSpace(int nq, int i) {
    return stressST_TestSpace[nq * i_dimension() + i];
  }

  VectorField DivStress(int nq, int j) {
    return divStressST[nq * j_dimension() + j];
  }
};

class SpaceTimeViscoElasticFaceElement {
  const DegreePair deg;
  const Quadrature &SpaceQ;
  const Quadrature &TimeQ;
  vector<double> qWeight;
  vector<Point> qLocal;
  vector<Point> qPoint;
  vector<Point> qNormal;
  vector<Point> qTangent;
  vector<Transformation> T;
  vector<bool> i_zero;
  vector<bool> j_zero;

  int dim;
  const Shape &SS;
  const Shape &SS2;
  const Shape &TS;
  const Shape &TestSpace;
  //const ScalarShapeValues &timeValue;
  //const ScalarShapeValues &testspaceValue;

  const cell &TC;
  double timedet;
  face f;
  int fid;
  const vector<vector<Point>> &faceqpoints;
  double ww;

  //ViscoElastic
  vector<VectorField> velocityST;
  vector<VectorField> velocityST_TestSpace;
  vector<Tensor> stressST;
  vector<Tensor> stressST_TestSpace;
  int numL;
public:
  SpaceTimeViscoElasticFaceElement(const STDiscretization &,
                                   const VectorMatrixBase &,
                                   const cell &,
                                   int f_id,
                                   int nL);

  const Point &LocalFaceCorner(int i) const {
    return TC.LocalCorner(TC.facecorner(fid, i));
  }

  const Point &QLocal(int q) const { return qLocal[q]; }

  const Point &QPoint(int q) const { return qPoint[q]; }

  const Point &QNormal(int q) const { return qNormal[q]; }

  const Point &QTangent(int q) const { return qTangent[q]; }

  double QWeight(int q) const { return qWeight[q]; }

  int nSpaceQ() const { return SpaceQ.size(); }

  int nTimeQ() const { return TimeQ.size(); }

  int nQ() const { return SpaceQ.size() * TimeQ.size(); }

  int td() const { //TensorDimension
    return std::max(1, dim + (dim - 1) + (dim - 2));
  }

  bool is_zero_test(int nq, int i) const { return i_zero[nq * i_dimension() + i]; }

  bool is_zero(int nq, int j) const { return j_zero[nq * j_dimension() + j]; }

  int i_dimension() const {
    return (SS.size() * (dim + td()) + SS2.size() * td() * numL) * TestSpace.size();
  }

  int j_dimension() const {
    return (SS.size() * (dim + td()) + SS2.size() * td() * numL) * TS.size();
  }

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
#endif //SPACETIMEVISCOELASTICELEMENT_HPP
