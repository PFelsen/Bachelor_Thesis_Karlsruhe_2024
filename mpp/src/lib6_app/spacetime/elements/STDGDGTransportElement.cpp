#include "STDGDGTransportElement.hpp"

const SpaceTimeShape &Get_SpaceTimeShape(const VectorMatrixBase &base, const Cell &c) {
  auto &disc = dynamic_cast<const STDiscretization &>(base.GetDisc());
  DegreePair degree = base.GetDoF().GetDegree(c);
  return disc.GetSpaceTimeShape(degree);
}

const SpaceTimeShape &Get_SpaceTimeShapeHighOrder(const VectorMatrixBase &base, const Cell &c) {
  auto &disc = dynamic_cast<const STDiscretization &>(base.GetDisc());
  return disc.GetSpaceTimeShape({6, 6});
}

const SpaceTimeQuadrature &Get_SpaceTimeQuadrature(const VectorMatrixBase &base, const Cell &c) {
  auto &disc = dynamic_cast<const STDiscretization &>(base.GetDisc());
  DegreePair degree = base.GetDoF().GetDegree(c);
  return disc.GetSpaceTimeQuad(degree);
}

const SpaceTimeQuadrature &Get_SpaceTimeQuadratureHighOrder(const VectorMatrixBase &base) {
  auto &disc = dynamic_cast<const STDiscretization &>(base.GetDisc());
  return disc.GetSpaceTimeQuad({6, 6});
}

STDGDGTransportElement::STDGDGTransportElement(const VectorMatrixBase &g,
                                               const cell &tc, bool max_quad_order) :
    shape(Get_SpaceTimeShape(g, *tc)),
    quad(Get_SpaceTimeQuadrature(g, *tc)),
    TC(tc), r(g.find_row(tc())), dim(tc.dim()),
    gradients(quad.size(), shape.size()) {
  qWeight.resize(nQ());
  qPoint.resize(nQ());
  densityST.resize(nQ() * j_dimension());
  dtDensityST.resize(nQ() * j_dimension());
  gradDensityST.resize(nQ() * j_dimension());

  const Transformation &T = TC.GetTransformation();
  for (int q = 0; q < quad.size(); q++) {
    qWeight[q] = T.Det() * quad.Weight(q);
    qPoint[q] = TC.LocalToGlobal(quad.QPoint(q));
    for (int i = 0; i < shape.size(); i++) {
      VectorField STGrad = shape.LocalSpaceTimeGradient(q, i);
      gradients[q][i] = T * STGrad;
    }
  }
  for (int tq = 0; tq < quad.GetTimeQuad().size(); ++tq) {
    for (int sq = 0; sq < quad.GetSpaceQuad().size(); ++sq) {
      int quad_i = tq * quad.GetSpaceQuad().size() + sq;
      for (int ti = 0; ti < shape.GetTimeShape().size(); ++ti) {
        for (int si = 0; si < shape.GetSpaceShape().size(); ++si) {
          int shape_i = ti * shape.GetSpaceShape().size() + si;
          int index = (tq * quad.GetSpaceQuad().size() + sq) * j_dimension()
                      + ti * shape.GetSpaceShape().size() + si;
          densityST[index] = shape.GetSpaceShape()(sq, si) * shape.GetTimeShape()(tq, ti);
          dtDensityST[index] = gradients[quad_i][shape_i].t();
          gradDensityST[index] = gradients[quad_i][shape_i].sliceTime();
        }
      }
    }
  }
}

double STDGDGTransportElement::Density(int nq, int j) {
  return densityST[nq * j_dimension() + j];
}

double STDGDGTransportElement::Density(const Point &globalPoint, const Vector &U){
  const Point localPoint = TC.GlobalToLocal(globalPoint);
  double P = 0.0;
  for (int j = 0; j < j_dimension(); ++j) {
    P += U(r, j) * shape(localPoint, j);
  }
  return P;
}

double STDGDGTransportElement::Density(int nq, const Vector &U) {
  double P = 0.0;
  for (int j = 0; j < j_dimension(); ++j) {
    P += U(r, j) * Density(nq, j);
  }
  return P;
}

double STDGDGTransportElement::DtDensity(int nq, int j) {
  return dtDensityST[nq * j_dimension() + j];
}

double STDGDGTransportElement::DtDensity(int nq, const Vector &u) {
  double dtP = 0;
  for (int j = 0; j < j_dimension(); ++j)
    dtP += u(r, j) * DtDensity(nq, j);
  return dtP;
}

VectorField STDGDGTransportElement::GradDensity(int nq, int j) {
  return gradDensityST[nq * j_dimension() + j];
}

VectorField STDGDGTransportElement::GradDensity(int nq, const Vector &u) {
  VectorField gradP = zero;
  for (int j = 0; j < j_dimension(); ++j)
    gradP += u(r, j) * GradDensity(nq, j);
  return gradP;
}
