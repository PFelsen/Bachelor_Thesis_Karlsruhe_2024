#include "STDGDGViscoAcousticElement.hpp"

const SpaceTimeShape &GetSpaceTimeShape(const VectorMatrixBase &base, const Cell &c) {
  auto &disc = dynamic_cast<const STDiscretization &>(base.GetDisc());
  DegreePair degree = base.GetDoF().GetDegree(c);
  return disc.GetSpaceTimeShape(degree);
}

const SpaceTimeShape &GetSpaceTimeShapeHighOrder(const VectorMatrixBase &base, const Cell &c) {
  auto &disc = dynamic_cast<const STDiscretization &>(base.GetDisc());
  return disc.GetSpaceTimeShape({6, 6});
}

const SpaceTimeQuadrature &GetSpaceTimeQuadrature(const VectorMatrixBase &base, const Cell &c) {
  auto &disc = dynamic_cast<const STDiscretization &>(base.GetDisc());
  DegreePair degree = base.GetDoF().GetDegree(c);
  return disc.GetSpaceTimeQuad(degree);
}

const SpaceTimeQuadrature &GetSpaceTimeQuadratureHighOrder(const VectorMatrixBase &base) {
  auto &disc = dynamic_cast<const STDiscretization &>(base.GetDisc());
  return disc.GetSpaceTimeQuad({6, 6});
}

STDGDGViscoAcousticElement::STDGDGViscoAcousticElement
    (const VectorMatrixBase &g, const cell &tc, int nL, bool max_quad_order)
    :
    shape(GetSpaceTimeShape(g, *tc)),
    quad(GetSpaceTimeQuadrature(g, *tc)),
    TC(tc), r(g.find_row(tc())), dimension(tc.dim()),
    gradients(quad.size(), shape.size()),
    dampingComponentCount(nL),
    space_quad(
        max_quad_order ? GetSpaceTimeQuadratureHighOrder(g).GetSpaceQuad() : quad.GetSpaceQuad()),
    component_to_index(tc.SpaceCell().dim() + 1 + nL, std::vector<int>{}) {
  qWeight.resize(nQ());
  qPoint.resize(nQ());
  velocityST.resize(nQ() * j_dimension());
  dtVelocityST.resize(nQ() * j_dimension());
  divVelocityST.resize(nQ() * j_dimension());

  pressureST.resize(nQ() * j_dimension());
  dtPressureST.resize(nQ() * j_dimension());
  gradPressureST.resize(nQ() * j_dimension());

  for (int i = 0; i < i_dimension(); i++) {
    component_to_index[GetComponent(i)].push_back(i);
  }


  const Transformation &T = TC.GetTransformation();


  det_time = T.DetTime();
  det_space = T.DetSpace();
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
          for (int k = 0; k < dimension; ++k) {
            int index = (tq * quad.GetSpaceQuad().size() + sq) * j_dimension()
                        + ti * shape.GetSpaceShape().size() * (1 + dampingComponentCount + dimension)
                        + k * shape.GetSpaceShape().size() + si;
            velocityST[index][k] = shape(quad_i, shape_i);
            dtVelocityST[index][k] = gradients[quad_i][shape_i].t();
            divVelocityST[index] = gradients[quad_i][shape_i][k];
          }
          for (int k = 0; k <= dampingComponentCount; ++k) {
            int index = (tq * quad.GetSpaceQuad().size() + sq) * j_dimension()
                        + ti * shape.GetSpaceShape().size() * (1 + dampingComponentCount + dimension)
                        + (dimension + k) * shape.GetSpaceShape().size()
                        + si;
            pressureST[index] = shape.GetSpaceShape()(sq, si) * shape.GetTimeShape()(tq, ti);
            dtPressureST[index] = gradients[quad_i][shape_i].t();
            gradPressureST[index] = gradients[quad_i][shape_i].sliceTime();
          }
        }
      }
    }
  }
}

int STDGDGViscoAcousticElement::nQ() const { return quad.size(); }

int STDGDGViscoAcousticElement::nSpaceQ() const { return space_quad.size(); }

int STDGDGViscoAcousticElement::nTimeQ() const { return quad.GetTimeQuad().size(); }

double STDGDGViscoAcousticElement::QWeight(int q) const { return qWeight[q]; }

double STDGDGViscoAcousticElement::SpaceQWeight(int qq) {
  if (space_qWeight.empty()) {
    space_qWeight.resize(space_quad.size());
    for (int q = 0; q < space_qWeight.size(); q++) {
      space_qWeight[q] = det_space * space_quad.Weight(q);
    }
  }
  return space_qWeight[qq];
}

Point STDGDGViscoAcousticElement::SpaceQPoint(int q) const {
  return space_quad.QPoint(q);
}
Point STDGDGViscoAcousticElement::TimeQPoint(int q) const {
  return quad.GetTimeQuad().QPoint(q);
}
// TimeQPoint
double STDGDGViscoAcousticElement::TimeQWeight(int qq) {
  if (time_qWeight.empty()) {
    time_qWeight.resize(quad.GetTimeQuad().size());
    for (int q = 0; q < time_qWeight.size(); q++) {
      time_qWeight[q] = det_time * quad.GetTimeQuad().Weight(q);
    }
  }
  return time_qWeight[qq];
}

Point STDGDGViscoAcousticElement::LocalTimeQPoint(int q) const {
  return quad.GetTimeQuad().QPoint(q);
}

const Point &STDGDGViscoAcousticElement::QPoint(int q) const { return qPoint[q]; }

int STDGDGViscoAcousticElement::i_dimension() const { return shape.size() * (1 + dimension + dampingComponentCount); }

double STDGDGViscoAcousticElement::Area() const {
  double a = 0;
  for (int q = 0; q < nQ(); ++q) a += qWeight[q];
  return a;
}

int STDGDGViscoAcousticElement::j_dimension() const { return shape.size() * (1 + dimension + dampingComponentCount); }

int STDGDGViscoAcousticElement::variable(int j) const {
  int ss_size = shape.GetSpaceShape().size();
  int tmp = j % (ss_size * (dimension + 1 + dampingComponentCount));
  if (tmp < dimension * ss_size) return 0;
  for (int nL = 0; nL <= dampingComponentCount; nL++)
    if (tmp - ss_size * dimension < (nL + 1) * ss_size) return 1 + nL;
  Exit("BLUPP");
}

int STDGDGViscoAcousticElement::GetComponent(int j) const {
  int ss_size = shape.GetSpaceShape().size();
  return (j % (ss_size * (dimension + 1 + dampingComponentCount))) / ss_size;
}

double STDGDGViscoAcousticElement::EvaluateComponentGlobal(const Point &p,
                                                           const Vector &u,
                                                           int component) const {
  Point localPoint = TC.GlobalToLocal(p);
  return EvaluateComponentLocal(localPoint, u, component);
}

double STDGDGViscoAcousticElement::EvaluateComponentLocal(const Point &p,
                                                          const Vector &u,
                                                          int component) const {
  double value = 0.0;

  for (int i : GetIndicesForComponent(component)) {
    int ti = i / (shape.GetSpaceShape().size() * (1 + dimension + dampingComponentCount));
    int si = i % (shape.GetSpaceShape().size() * (1 + dimension + dampingComponentCount));
    int ssj = si % shape.GetSpaceShape().size();
    int shape_iddx = ti * shape.GetSpaceShape().size() + ssj;
    value += u(r, i) * shape(p, shape_iddx);
  }

  return value;
}

double STDGDGViscoAcousticElement::EvaluateComponent(int nq, const Vector &u,
                                                     int component) const {
  double value = 0.0;

  for (int i : GetIndicesForComponent(component)) {
    int ti = i / (shape.GetSpaceShape().size() * GetComponentCount());
    int si = i % (shape.GetSpaceShape().size() * GetComponentCount());
    int ssj = si % shape.GetSpaceShape().size();
    int shape_iddx = ti * shape.GetSpaceShape().size() + ssj;
    value += u(r, i) * shape(nq, shape_iddx);
  }

  return value;
}

VectorField STDGDGViscoAcousticElement::Velocity(int nq, int j) const {
  return velocityST[nq * j_dimension() + j];
}

VectorField STDGDGViscoAcousticElement::Velocity(int nq, const Vector &U) const {
  VectorField V = zero;
  for (int j = 0; j < j_dimension(); ++j) {
    V += U(r, j) * Velocity(nq, j);
  }
  return V;
}

VectorField STDGDGViscoAcousticElement::VelocityLocal(const Point &localPoint, int j) const {
  int var = variable(j);
  if (var != 0) return zero;
  int ti = j / (shape.GetSpaceShape().size() * (1 + dimension + dampingComponentCount));
  int si = j % (shape.GetSpaceShape().size() * (1 + dimension + dampingComponentCount));
  int vComponent = si / shape.GetSpaceShape().size();
  int ssj = si % shape.GetSpaceShape().size();
  int shape_iddx = ti * shape.GetSpaceShape().size() + ssj;

  VectorField V = zero;
  V[vComponent] = shape(localPoint, shape_iddx);
  return V;
}

VectorField
STDGDGViscoAcousticElement::VelocityLocal(const Point &localPoint, const Vector &u) const {
  VectorField V = zero;
  for (int j = 0; j < j_dimension(); ++j)
    V += u(r, j) * VelocityLocal(localPoint, j);
  return V;
}

VectorField
STDGDGViscoAcousticElement::VelocityGlobal(const Point &globalPoint, const Vector &u) const {
  Point localPoint = TC.GlobalToLocal(globalPoint);
  return VelocityLocal(localPoint, u);
}

VectorField STDGDGViscoAcousticElement::DtVelocity(int nq, int j) const {
  return dtVelocityST[nq * j_dimension() + j];
}

VectorField STDGDGViscoAcousticElement::DtVelocity(int nq, const Vector &u) const {
  VectorField V = zero;
  for (int j = 0; j < j_dimension(); ++j)
    V += u(r, j) * DtVelocity(nq, j);
  return V;
}

double STDGDGViscoAcousticElement::DivVelocity(int nq, int j) const {
  return divVelocityST[nq * j_dimension() + j];
}

double STDGDGViscoAcousticElement::DivVelocity(int nq, const Vector &u) const {
  double divV = 0.0;
  for (int j = 0; j < j_dimension(); ++j)
    divV += u(r, j) * DivVelocity(nq, j);
  return divV;
}

double STDGDGViscoAcousticElement::Pressure(int nq, int j) const {
  return pressureST[nq * j_dimension() + j];
}

double STDGDGViscoAcousticElement::Pressure(int nq, const Vector &U) const {
  double P = 0.0;
  for (int j : component_to_index[TC.SpaceCell().dim()]) {
    P += U(r, j) * Pressure(nq, j);
  }
  return P;
}

DampingVector STDGDGViscoAcousticElement::DampingPressure(int nq, const Vector &U) const {
  DampingVector DP(dampingComponentCount);
  const auto spaceDimension = TC.SpaceCell().dim();
  for (int i = 0; i < dampingComponentCount; ++i) {
    for (int j : component_to_index[spaceDimension + 1 + i]) {
      DP[i] += U(r, j) * Pressure(nq, j);
    }
  }
  return DP;
}
double STDGDGViscoAcousticElement::EvaluateTimeLocal(int q, int i) const {
  return shape.GetTimeShape()(q,i);
}
double STDGDGViscoAcousticElement::PressureLocal(const Point &localPoint, int j) const {
  int var = variable(j);
  if (var == 0) return 0.0;
  int ti = j / (shape.GetSpaceShape().size() * (1 + dimension + dampingComponentCount));
  int si = j % (shape.GetSpaceShape().size() * (1 + dimension + dampingComponentCount));
  int ssj = si % shape.GetSpaceShape().size();
  int shape_iddx = ti * shape.GetSpaceShape().size() + ssj;

  return shape(localPoint, shape_iddx);
}

double STDGDGViscoAcousticElement::PressureLocal(const Point &localPoint, const Vector &u) const {
  double P = 0;
  for (int j = 0; j < j_dimension(); ++j)
    P += u(r, j) * PressureLocal(localPoint, j);
  return P;
}

double STDGDGViscoAcousticElement::PressureGlobal(const Point &globalPoint, int i) const {
  Point localPoint = TC.GlobalToLocal(globalPoint);
  return PressureLocal(localPoint, i);
}

double STDGDGViscoAcousticElement::PressureGlobal(const Point &globalPoint, const Vector &u) const {
  Point localPoint = TC.GlobalToLocal(globalPoint);
  return PressureLocal(localPoint, u);
}

double STDGDGViscoAcousticElement::DtPressure(int nq, int j) const {
  return dtPressureST[nq * j_dimension() + j];
}

double STDGDGViscoAcousticElement::DtPressure(int nq, const Vector &u) const {
  double dtP = 0;
  for (int j = 0; j < j_dimension(); ++j)
    dtP += u(r, j) * DtPressure(nq, j);
  return dtP;
}

VectorField STDGDGViscoAcousticElement::GradPressure(int nq, int j) const {
  return gradPressureST[nq * j_dimension() + j];
}

VectorField STDGDGViscoAcousticElement::GradPressure(int nq, const Vector &u) const {
  VectorField gradP = zero;
  for (int j = 0; j < j_dimension(); ++j)
    gradP += u(r, j) * GradPressure(nq, j);
  return gradP;
}
int STDGDGViscoAcousticElement::GetTimeDimension() const {return shape.GetTimeShape().size();}
