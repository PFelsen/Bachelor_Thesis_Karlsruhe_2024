#ifndef SPACETIME_STDGDGVISCOACOUSTICFACEELEMENTGLGL_HPP
#define SPACETIME_STDGDGVISCOACOUSTICFACEELEMENTGLGL_HPP

#include "Quadrature.hpp"
#include "Shapes.hpp"
#include "SpaceTimeDiscretization.hpp"
#include "VectorMatrixBase.hpp"
#include "Vector.hpp"
#include "STDGDGViscoAcousticElementGLGL.hpp"

class STDGDGViscoAcousticFaceElementGLGL {

  const DegreePair deg;

  double det_time;

  vector<double> space_qWeight;
  vector<double> time_qWeight;
  vector<double> qWeight;
  vector<Point> qPoint;
  vector<Point> qLocal;

  int dim;
  const SpaceTimeShape &shape;

  const Quadrature &SpaceQ;
  const Quadrature &TimeQ;
  vector<Point> qNormal;
  vector<bool> i_zero;
  vector<bool> j_zero;

  const cell &TC;
  const row r;
  face f;
  int fid;
  const vector<vector<Point>> &faceqpoints;

  vector<double> pressureST;
  vector<VectorField> velocityST;
  int numL;
public:
  STDGDGViscoAcousticFaceElementGLGL(const STDiscretization &,
                                     const VectorMatrixBase &,
                                     const cell &,
                                     int f_id,
                                     int nL);

  STDGDGViscoAcousticFaceElementGLGL(const STDiscretization &,
                                     const VectorMatrixBase &,
                                     const cell &,
                                     int f_id,
                                     int nL,
                                     string dgdg);

  STDGDGViscoAcousticFaceElementGLGL(const STDiscretization &,
                                     const VectorMatrixBase &,
                                     const cell &,
                                     int f_id,
                                     int nL,
                                     string dgdg,
                                     const cell &);

  const Point &QLocal(int q) const;

  const Point &QPoint(int q) const;

  const Point &QNormal(int q) const;

  double QWeight(int q) const;

  int nSpaceQ() const;

  int nTimeQ() const;

  int nQ() const;

  bool is_zero_test(int nq, int i) const;

  bool is_zero(int nq, int j) const;

  int i_dimension() const;

  int j_dimension() const;

  int variable(int j) const;

  VectorField Velocity(int nq, int j);

  double Pressure(int nq, int j);

  double Pressure(int nq, const double *u);

  double Pressure(int nq, const Vector &U) {
    double P = 0.0;
    for (int j = 0; j < j_dimension(); ++j) {
      P += U(r, j) * Pressure(nq, j);
    }
    return P;
  }

  VectorField Velocity(Point QP, int i);

  VectorField Velocity(int nq, const double *u);

  VectorField Velocity(int nq, const Vector &U) {
    VectorField V = zero;
    for (int j = 0; j < j_dimension(); ++j) {
      V += U(r, j) * Velocity(nq, j);
    }
    return V;
  }

  double Pressure(Point QP, int i);
};


#endif //SPACETIME_STDGDGVISCOACOUSTICFACEELEMENT_HPP
