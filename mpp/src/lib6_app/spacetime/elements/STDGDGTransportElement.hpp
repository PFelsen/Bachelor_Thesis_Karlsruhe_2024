#ifndef SPACETIME_STDGDGTRANSPORTELEMENT_HPP
#define SPACETIME_STDGDGTRANSPORTELEMENT_HPP

#include "SpaceTimeDiscretization.hpp"
#include "SpaceTimeQuadrature.hpp"
#include "SpaceTimeShape.hpp"
#include "Vector.hpp"

class STDGDGTransportElement {
public:
  vector<double> qWeight;
  vector<Point> qPoint;
  vector<Point> qLocal;
  vector<Point> qNormal;

  int dim;
  const SpaceTimeShape &shape;
  const SpaceTimeQuadrature &quad;
  ShapeValues<VectorField> gradients;

  const cell &TC;
  const row r;

  std::vector<double> densityST;
  std::vector<double> dtDensityST;
  std::vector<VectorField> gradDensityST;

  STDGDGTransportElement(const VectorMatrixBase &,
                         const cell &,
                         bool max_quad_order = false);


  int nQ() const { return quad.size(); }

  double QWeight(int q) const { return qWeight[q]; }

  const Point &QPoint(int q) const { return qPoint[q]; }

  int i_dimension() const { return shape.size(); }

  int j_dimension() const { return shape.size(); }

  double Density(int nq, int j);

  double Density(const Point &globalPoint, const Vector &U);

  double Density(int nq, const Vector &U);

  double DtDensity(int nq, int j);

  double DtDensity(int nq, const Vector &u);

  VectorField GradDensity(int nq, int j);

  VectorField GradDensity(int nq, const Vector &u);
};

using SpaceTimeTransportDGTElement = STDGDGTransportElement;

#endif 
