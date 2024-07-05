#include "STDGDGViscoAcousticElementGLGL.hpp"


STDGDGViscoAcousticElementGLGL::STDGDGViscoAcousticElementGLGL
    (const VectorMatrixBase &g, const cell &tc, int nL, bool max_quad_order)
    :
    shape(GetSpaceTimeShape(g, *tc)),
    quad(GetSpaceTimeQuadrature(g, *tc)),
    TC(tc), r(g.find_row(tc())), dim(tc.dim()),
    gradients(quad.size(), shape.size()),
    numL(nL),
    space_quad(
        max_quad_order ? GetSpaceTimeQuadratureHighOrder(g).GetSpaceQuad() : quad.GetSpaceQuad()) {
  qWeight.resize(nQ());
  qPoint.resize(nQ());
  component_to_index.resize(1 + tc.SpaceCell().dim() + numL);
  //mout << quad.Name() << quad.size() << DOUT(quad.GetSpaceQuad().Name()) << DOUT(quad.GetTimeQuad().Name()) << endl;
  velocityST.resize(nQ() * j_dimension());
  dtVelocityST.resize(nQ() * j_dimension());
  divVelocityST.resize(nQ() * j_dimension());

  pressureST.resize(nQ() * j_dimension());
  dtPressureST.resize(nQ() * j_dimension());
  gradPressureST.resize(nQ() * j_dimension());

  for(int i = 0; i < i_dimension(); i++){
    component_to_index[GetComponent(i)].push_back(i);
  }


  const auto &T = TC.GetTransformation();

  det_time = T.DetTime();
  det_space = T.DetSpace();
  for (int q = 0; q < quad.size(); q++) {
    qWeight[q] = T.Det() * quad.Weight(q);
    qPoint[q] = TC.LocalToGlobal(quad.QPoint(q));
    for (int i = 0; i < shape.size(); i++) {
      VectorFieldT<double, SpaceDimension, 1> STGrad = shape.LocalSpaceTimeGradient(q, i);
      gradients[q][i] = T * STGrad;
    }
  }
  for (int tq = 0; tq < quad.GetTimeQuad().size(); ++tq) {
    for (int sq = 0; sq < quad.GetSpaceQuad().size(); ++sq) {
      int quad_i = tq * quad.GetSpaceQuad().size() + sq;
      for (int ti = 0; ti < shape.GetTimeShape().size(); ++ti) {
        for (int si = 0; si < shape.GetSpaceShape().size(); ++si) {
          int shape_i = ti * shape.GetSpaceShape().size() + si;
          for (int k = 0; k < dim; ++k) {
            int index = (tq * quad.GetSpaceQuad().size() + sq) * j_dimension()
                        + ti * shape.GetSpaceShape().size() * (1 + numL + dim)
                        + k * shape.GetSpaceShape().size() + si;
            velocityST[index][k] = shape(quad_i, shape_i);
            dtVelocityST[index][k] = gradients[quad_i][shape_i].t();
            divVelocityST[index] = gradients[quad_i][shape_i][k];
          }
          for (int k = 0; k <= numL; ++k) {
            int index = (tq * quad.GetSpaceQuad().size() + sq) * j_dimension()
                        + ti * shape.GetSpaceShape().size() * (1 + numL + dim)
                        + (dim + k) * shape.GetSpaceShape().size()
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

int STDGDGViscoAcousticElementGLGL::nQ() const { return quad.size(); }

int STDGDGViscoAcousticElementGLGL::nSpaceQ() const { return space_quad.size(); }

int STDGDGViscoAcousticElementGLGL::nTimeQ() const { return quad.GetTimeQuad().size(); }

double STDGDGViscoAcousticElementGLGL::QWeight(int q) const { return qWeight[q]; }

double STDGDGViscoAcousticElementGLGL::SpaceQWeight(int qq) {
  if (space_qWeight.empty()) {
    space_qWeight.resize(space_quad.size());
    for (int q = 0; q < space_qWeight.size(); q++) {
      space_qWeight[q] = det_space * space_quad.Weight(q);
    }
  }
  return space_qWeight[qq];
}

Point STDGDGViscoAcousticElementGLGL::SpaceQPoint(int q) const {
  return space_quad.QPoint(q);
}

double STDGDGViscoAcousticElementGLGL::TimeQWeight(int qq) {
  if (time_qWeight.empty()) {
    time_qWeight.resize(quad.GetTimeQuad().size());
    for (int q = 0; q < time_qWeight.size(); q++) {
      time_qWeight[q] = det_time * quad.GetTimeQuad().Weight(q);
    }
  }
  return time_qWeight[qq];
}

const Point &STDGDGViscoAcousticElementGLGL::QPoint(int q) const { return qPoint[q]; }

int STDGDGViscoAcousticElementGLGL::i_dimension() const { return shape.size() * (1 + dim + numL); }

double STDGDGViscoAcousticElementGLGL::Area() const {
  double a = 0;
  for (int q = 0; q < nQ(); ++q) a += qWeight[q];
  return a;
}

int STDGDGViscoAcousticElementGLGL::j_dimension() const { return shape.size() * (1 + dim + numL); }

int STDGDGViscoAcousticElementGLGL::variable(int j) const {
  int ss_size = shape.GetSpaceShape().size();
  int tmp = j % (ss_size * (dim + 1 + numL));
  if (tmp < dim * ss_size) return 0;
  for (int nL = 0; nL <= numL; nL++)
    if (tmp - ss_size * dim < (nL + 1) * ss_size) return 1 + nL;
  Exit("BLUPP");
}

int STDGDGViscoAcousticElementGLGL::GetComponent(int j) const {
  int ss_size = shape.GetSpaceShape().size();
  return (j % (ss_size * (dim + 1 + numL))) / ss_size;
}

double STDGDGViscoAcousticElementGLGL::EvaluateComponentGlobal(const Point &p,
                                                           const Vector &u,
                                                           int component) const {
  Point localPoint = TC.GlobalToLocal(p);
  return EvaluateComponentLocal(localPoint, u, component);
}

double STDGDGViscoAcousticElementGLGL::EvaluateComponentLocal(const Point &p,
                                                          const Vector &u,
                                                          int component) const {
  double value = 0.0;

  for (int i : GetIndicesForComponent(component)) {
    int ti = i / (shape.GetSpaceShape().size() * (1 + dim + numL));
    int si = i % (shape.GetSpaceShape().size() * (1 + dim + numL));
    int ssj = si % shape.GetSpaceShape().size();
    int shape_iddx = ti * shape.GetSpaceShape().size() + ssj;
    value += u(r, i) * shape(p, shape_iddx);
  }

  return value;
}

double STDGDGViscoAcousticElementGLGL::EvaluateComponent(int nq, const Vector &u,
                                                     int component) const {
  double value = 0.0;

  for (int i : GetIndicesForComponent(component)) {
    int ti = i / (shape.GetSpaceShape().size() * (1 + dim + numL));
    int si = i % (shape.GetSpaceShape().size() * (1 + dim + numL));
    int ssj = si % shape.GetSpaceShape().size();
    int shape_iddx = ti * shape.GetSpaceShape().size() + ssj;
    value += u(r, i) * shape(nq, shape_iddx);
  }

  return value;
}

VectorField STDGDGViscoAcousticElementGLGL::Velocity(int nq, int j) {
  return velocityST[nq * j_dimension() + j];
}

VectorField STDGDGViscoAcousticElementGLGL::Velocity(int nq, const Vector &U) {
  VectorField V = zero;
  for (int j = 0; j < j_dimension(); ++j) {
    V += U(r, j) * Velocity(nq, j);
  }
  return V;
}

VectorField STDGDGViscoAcousticElementGLGL::VelocityLocal(const Point &localPoint, int j) {
  int var = variable(j);
  if (var != 0) return zero;
  int ti = j / (shape.GetSpaceShape().size() * (1 + dim + numL));
  int si = j % (shape.GetSpaceShape().size() * (1 + dim + numL));
  int vComponent = si / shape.GetSpaceShape().size();
  int ssj = si % shape.GetSpaceShape().size();
  int shape_iddx = ti * shape.GetSpaceShape().size() + ssj;

  VectorField V = zero;
  V[vComponent] = shape(localPoint, shape_iddx);
  return V;
}

VectorField STDGDGViscoAcousticElementGLGL::VelocityLocal(const Point &localPoint, const Vector &u) {
  VectorField V = zero;
  for (int j = 0; j < j_dimension(); ++j)
    V += u(r, j) * VelocityLocal(localPoint, j);
  return V;
}

VectorField STDGDGViscoAcousticElementGLGL::VelocityGlobal(const Point &globalPoint, const Vector &u) {
  Point localPoint = TC.GlobalToLocal(globalPoint);
  return VelocityLocal(localPoint, u);
}

VectorField STDGDGViscoAcousticElementGLGL::DtVelocity(int nq, int j) {
  return dtVelocityST[nq * j_dimension() + j];
}

VectorField STDGDGViscoAcousticElementGLGL::DtVelocity(int nq, const Vector &u) {
  VectorField V = zero;
  for (int j = 0; j < j_dimension(); ++j)
    V += u(r, j) * DtVelocity(nq, j);
  return V;
}

double STDGDGViscoAcousticElementGLGL::DivVelocity(int nq, int j) {
  return divVelocityST[nq * j_dimension() + j];
}

double STDGDGViscoAcousticElementGLGL::DivVelocity(int nq, const Vector &u) {
  double divV = 0.0;
  for (int j = 0; j < j_dimension(); ++j)
    divV += u(r, j) * DivVelocity(nq, j);
  return divV;
}

double STDGDGViscoAcousticElementGLGL::Pressure(int nq, int j) {
  return pressureST[nq * j_dimension() + j];
}

double STDGDGViscoAcousticElementGLGL::Pressure(int nq, const Vector &U) {
  double P = 0.0;
  for (int j = 0; j < j_dimension(); ++j) {
    P += U(r, j) * Pressure(nq, j);
  }
  return P;
}

double STDGDGViscoAcousticElementGLGL::PressureLocal(const Point &localPoint, int j) {
  int var = variable(j);
  if (var == 0) return 0.0;
  int ti = j / (shape.GetSpaceShape().size() * (1 + dim + numL));
  int si = j % (shape.GetSpaceShape().size() * (1 + dim + numL));
  int ssj = si % shape.GetSpaceShape().size();
  int shape_iddx = ti * shape.GetSpaceShape().size() + ssj;

  return shape(localPoint, shape_iddx);
}

double STDGDGViscoAcousticElementGLGL::PressureLocal(const Point &localPoint, const Vector &u) {
  double P = 0;
  for (int j = 0; j < j_dimension(); ++j)
    P += u(r, j) * PressureLocal(localPoint, j);
  return P;
}

double STDGDGViscoAcousticElementGLGL::PressureGlobal(const Point &globalPoint, int i) {
  Point localPoint = TC.GlobalToLocal(globalPoint);
  return PressureLocal(localPoint, i);
}

double STDGDGViscoAcousticElementGLGL::PressureGlobal(const Point &globalPoint, const Vector &u) {
  Point localPoint = TC.GlobalToLocal(globalPoint);
  return PressureLocal(localPoint, u);
}

double STDGDGViscoAcousticElementGLGL::DtPressure(int nq, int j) {
  return dtPressureST[nq * j_dimension() + j];
}

double STDGDGViscoAcousticElementGLGL::DtPressure(int nq, const Vector &u) {
  double dtP = 0;
  for (int j = 0; j < j_dimension(); ++j)
    dtP += u(r, j) * DtPressure(nq, j);
  return dtP;
}

VectorField STDGDGViscoAcousticElementGLGL::GradPressure(int nq, int j) {
  return gradPressureST[nq * j_dimension() + j];
}

VectorField STDGDGViscoAcousticElementGLGL::GradPressure(int nq, const Vector &u) {
  VectorField gradP = zero;
  for (int j = 0; j < j_dimension(); ++j)
    gradP += u(r, j) * GradPressure(nq, j);
  return gradP;
}