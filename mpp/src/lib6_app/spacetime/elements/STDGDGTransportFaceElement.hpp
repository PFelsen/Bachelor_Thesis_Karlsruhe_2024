#ifndef SPACETIME_STDGDGTRANSPORTFACEELEMENT_HPP
#define SPACETIME_STDGDGTRANSPORTFACEELEMENT_HPP

#include "Quadrature.hpp"
#include "Shapes.hpp"
#include "SpaceTimeDiscretization.hpp"
#include "VectorMatrixBase.hpp"
#include "Vector.hpp"

class SpaceTimeTransportDGTFaceElement {
  const DegreePair deg;
  const Quadrature &SpaceQ;
  const Quadrature &TimeQ;
  vector<double> qWeight;
  vector<Point> qPoint;
  vector<Point> qLocal;
  vector<Point> qNormal;

  int dim;
  const SpaceTimeShape &shape;

  const cell &TC;
  const row r;
  face f;
  int fid;
  const vector<vector<Point>> &faceqpoints;

  vector<double> densityST;

public:
  SpaceTimeTransportDGTFaceElement(const STDiscretization &,
                                   const VectorMatrixBase &,
                                   const cell &,
                                   int f_id);

  SpaceTimeTransportDGTFaceElement(const STDiscretization &,
                                   const VectorMatrixBase &,
                                   const cell &,
                                   int f_id,
                                   string);

  SpaceTimeTransportDGTFaceElement(const STDiscretization &,
                                   const VectorMatrixBase &,
                                   const cell &,
                                   int f_id,
                                   string,
                                   const cell &);

  const Point &QPoint(int q) const { return qPoint[q]; }

  const Point &QLocal(int q) const { return qLocal[q]; }

  const Point &QNormal(int q) const { return qNormal[q]; }

  double QWeight(int q) const { return qWeight[q]; }

  int nQ() const { return SpaceQ.size() * TimeQ.size(); }

  int i_dimension() const { return shape.size(); }

  int j_dimension() const { return shape.size(); }

  double Density(int nq, int j) {
    return densityST[nq * j_dimension() + j];
  }

  double Density(int nq, const Vector &U) {
    double P = 0.0;
    for (int j = 0; j < j_dimension(); ++j) {
      P += U(r, j) * Density(nq, j);
    }
    return P;
  }
};


#endif //SPACETIME_STDGDGVISCOACOUSTICFACEELEMENT_HPP
