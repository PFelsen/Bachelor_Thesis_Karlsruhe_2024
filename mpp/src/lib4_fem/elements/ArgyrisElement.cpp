#include "ArgyrisElement.hpp"

template<typename TT, int sDim, int tDim>
void ArgyrisBasicElementT<TT, sDim, tDim>::H20BC(const std::list<Point> &corners, Vector &u) {
  u.ClearDirichletFlags();
  for (cell c = u.cells(); c != u.cells_end(); ++c) {
    H20BC(corners, u, *c);
  }
  u.DirichletConsistent();
}

template<typename TT, int sDim, int tDim>
void ArgyrisBasicElementT<TT, sDim, tDim>::H20BC(const std::list<Point> &corners, Vector &u,
                                                 const Cell &C) {
  if (!u.OnBoundary(C)) return;
  for (int face = 0; face < C.Faces(); ++face) {
    if (!u.OnBoundary(C, face)) continue;
    rows R(u.GetMatrixGraph(), C, face);
    for (const row &r : R) {
      if (r.NumberOfDofs() == 1) {
        // Nodal point is on face
        u(r, 0) = 0.0;
        u.D(r, 0) = true;
      } else {
        // Nodal point is a corner
        if (std::find(corners.begin(), corners.end(), r()) != corners.end()) {
          // nodal point coincides with a corner of domain
          for (int k = 0; k < r.NumberOfDofs(); ++k) {
            u(r, k) = 0.0;
            u.D(r, k) = true;
          }
        } else {
          // nodal point is on the boundary but not a corner of domain
          for (int k = 0; k < 3; ++k) {
            u(r, k) = 0.0;
            u.D(r, k) = true;
          }
          Point tangent = C.FaceCorner(face, 1) - C.FaceCorner(face, 0);
          if (abs(tangent * Point(1.0, 0.0)) < 1.0e-8) {
            // face parallel to y axis
            u(r, 4) = 0.0;
            u.D(r, 4) = true;
            u(r, 5) = 0.0;
            u.D(r, 5) = true;
          } else if (abs(tangent * Point(0.0, 1.0)) < 1.0e-8) {
            // face parallel to x axis
            u(r, 3) = 0.0;
            u.D(r, 3) = true;
            u(r, 4) = 0.0;
            u.D(r, 4) = true;
          } else {
            // face not parallel to any axis
            double lambda = tangent[1] / tangent[0];
            u.BCadd(r, 3, 4, lambda);
            u.BCadd(r, 5, 4, 1.0 / lambda);
          }
        }
      }
    }
  }
}

template<typename TT, int sDim, int tDim>
ArgyrisBasicElementT<TT, sDim, tDim>::ArgyrisBasicElementT(const VectorMatrixBase &base,
                                                           const Cell &c) :
    ElementT<TT, sDim, tDim>(base, c) {
  init();
}

template<typename TT, int sDim, int tDim>
ArgyrisBasicElementT<TT, sDim, tDim>::ArgyrisBasicElementT(
    const VectorMatrixBase &base, const BasicElementT<TT, sDim, tDim> &baseElement) :
    ElementT<TT, sDim, tDim>(base, baseElement) {
  init();
}

template<typename TT, int sDim, int tDim>
std::vector<TT>
ArgyrisBasicElementT<TT, sDim, tDim>::values(const PointT<TT, sDim, tDim> &z) const {
  std::vector<TT> v(21);
  for (int j = 0; j < 21; ++j)
    v[j] = this->shape(z, j);
  std::vector<TT> V(21);
  for (int j = 0; j < 21; ++j)
    applyC(this->C0, this->C1, this->C2, this->C3, this->C4, this->C5, V[j], v, j);
  return V;
}

template<typename TT, int sDim, int tDim>
std::vector<VectorFieldT<TT, sDim>>
ArgyrisBasicElementT<TT, sDim, tDim>::derivatives(const PointT<TT, sDim, tDim> &z) const {
  vector<VectorFieldT<TT, sDim>> G(21);
  for (int j = 0; j < 21; ++j)
    G[j] = this->GetTransformation(0) * this->shape.LocalGradient(z, j);
  std::vector<VectorFieldT<TT, sDim>> D(21);
  for (int j = 0; j < 21; ++j)
    applyC(this->C0, this->C1, this->C2, this->C3, this->C4, this->C5, D[j], G, j);
  return D;
}

template<typename TT, int sDim, int tDim>
std::vector<SymTensorT<TT, sDim>>
ArgyrisBasicElementT<TT, sDim, tDim>::hessians(const PointT<TT, sDim, tDim> &z) const {
  vector<SymTensorT<TT, sDim>> H(21);
  for (int j = 0; j < 21; ++j)
    H[j] = this->multiplyTheta(this->shape.LocalHessian(z, j));
  vector<SymTensorT<TT, sDim>> D(21);
  for (int j = 0; j < 21; ++j)
    applyC(this->C0, this->C1, this->C2, this->C3, this->C4, this->C5, D[j], H, j);
  return D;
}

template<typename TT, int sDim, int tDim>
int ArgyrisBasicElementT<TT, sDim, tDim>::get_maxk(int i) const {
  if (i < 3) return 6;
  return 1;
}

template<typename TT, int sDim, int tDim>
int ArgyrisBasicElementT<TT, sDim, tDim>::indexing_k(int i, int k, int m) const {
  if (i < 3) return 6 * m + k;
  return m;
}

template<typename TT, int sDim, int tDim>
TT ArgyrisBasicElementT<TT, sDim, tDim>::getC(int i, int j) const {
  if ((i == 0 && j == 0) || (i == 6 && j == 6) || (i == 12 && j == 12)) return TT(1.0);
  if ((i == 1 && j == 1) || (i == 7 && j == 7) || (i == 13 && j == 13)) return C0[0][0];
  if ((i == 2 && j == 1) || (i == 8 && j == 7) || (i == 14 && j == 13)) return C0[1][0];
  if ((i == 1 && j == 2) || (i == 7 && j == 8) || (i == 13 && j == 14)) return C0[0][1];
  if ((i == 2 && j == 2) || (i == 8 && j == 8) || (i == 14 && j == 14)) return C0[1][1];
  if ((i == 3 && j == 3) || (i == 9 && j == 9) || (i == 15 && j == 15)) return C1[0][0];
  if ((i == 4 && j == 3) || (i == 10 && j == 9) || (i == 16 && j == 15)) return C1[1][0];
  if ((i == 5 && j == 3) || (i == 11 && j == 9) || (i == 17 && j == 15)) return C1[2][0];
  if ((i == 3 && j == 4) || (i == 9 && j == 10) || (i == 15 && j == 16)) return C1[0][1];
  if ((i == 4 && j == 4) || (i == 10 && j == 10) || (i == 16 && j == 16)) return C1[1][1];
  if ((i == 5 && j == 4) || (i == 11 && j == 10) || (i == 17 && j == 16)) return C1[2][1];
  if ((i == 3 && j == 5) || (i == 9 && j == 11) || (i == 15 && j == 17)) return C1[0][2];
  if ((i == 4 && j == 5) || (i == 10 && j == 11) || (i == 16 && j == 17)) return C1[1][2];
  if ((i == 5 && j == 5) || (i == 11 && j == 11) || (i == 17 && j == 17)) return C1[2][2];
  if (i == 18) {
    if (0 <= j && j < 6) return C2[0][j];
    if (6 <= j && j < 12) return C3[0][j - 6];
  }
  if (i == 19) {
    if (6 <= j && j < 12) return C3[1][j - 6];
    if (12 <= j && j < 18) return C4[0][j - 12];
  }
  if (i == 20) {
    if (0 <= j && j < 6) return C2[1][j];
    if (12 <= j && j < 18) return C4[1][j - 12];
  }
  if ((i == 18 && j == 18) || (i == 19 && j == 19) || (i == 20 && j == 20)) return C5[i - 18];
  return TT(0.0);
}

template<typename TT, int sDim, int tDim>
void ArgyrisBasicElementT<TT, sDim, tDim>::init() {
  if (this->c.dim() != 2 || this->c.ReferenceType() != TRIANGLE) Exit("Error in Argyris element!");

  const TransformationT<TT, sDim> &trafo = this->GetTransformation(0);

  theta[0][0] = trafo[0][0] * trafo[0][0];
  theta[0][1] = 2 * trafo[0][0] * trafo[0][1];
  theta[0][2] = trafo[0][1] * trafo[0][1];
  theta[1][0] = trafo[0][0] * trafo[1][0];
  theta[1][1] = trafo[0][0] * trafo[1][1] + trafo[0][1] * trafo[1][0];
  theta[1][2] = trafo[0][1] * trafo[1][1];
  theta[2][0] = trafo[1][0] * trafo[1][0];
  theta[2][1] = 2 * trafo[1][0] * trafo[1][1];
  theta[2][2] = trafo[1][1] * trafo[1][1];

  const PointT<TT, sDim, tDim> &z0 = this->c[0];
  const PointT<TT, sDim, tDim> &z1 = this->c[1];
  const PointT<TT, sDim, tDim> &z2 = this->c[2];
  PointT<TT, sDim, tDim> v[3];
  TT l[3];

  v[0] = z1 - z0;
  v[1] = z2 - z1;
  v[2] = z0 - z2;
  l[0] = v[0] * v[0];
  l[1] = v[1] * v[1];
  l[2] = v[2] * v[2];
  for (int i = 0; i < 3; ++i) {
    normal[i][0] = v[i][1];
    normal[i][1] = -v[i][0];
    normal[i] /= sqrt(l[i]);
    if ((this->c.Face(i) + mid(normal[i])) > this->c.Face(i)) sign[i] = TT(1.0);
    else sign[i] = TT(-1.0);
  }

  C0[0][0] = trafo(0)[0];
  C0[0][1] = trafo(1)[0];
  C0[1][0] = trafo(0)[1];
  C0[1][1] = trafo(1)[1];

  C1[0][0] = trafo(0)[0] * trafo(0)[0];
  C1[0][1] = 2 * trafo(0)[0] * trafo(1)[0];
  C1[0][2] = trafo(1)[0] * trafo(1)[0];
  C1[1][0] = trafo(0)[1] * trafo(0)[0];
  C1[1][1] = trafo(0)[1] * trafo(1)[0] + trafo(0)[0] * trafo(1)[1];
  C1[1][2] = trafo(1)[0] * trafo(1)[1];
  C1[2][0] = trafo(0)[1] * trafo(0)[1];
  C1[2][1] = 2 * trafo(1)[1] * trafo(0)[1];
  C1[2][2] = trafo(1)[1] * trafo(1)[1];

  TT c18 = -(trafo(0)[1] * v[0][0] + trafo(1)[1] * v[0][1]) / l[0];
  TT c19 = ((trafo(0)[0] + trafo(0)[1]) * v[1][0] + (trafo(1)[0] + trafo(1)[1]) * v[1][1]) / l[1]
           / sqrt(TT(2.0));
  TT c20 = -(trafo(0)[0] * v[2][0] + trafo(1)[0] * v[2][1]) / l[2];

  C3[0][0] = (15.0 / TT(8.0)) * c18;
  C2[0][0] = -C3[0][0];
  C4[0][0] = (15.0 / TT(8.0)) * c19;
  C3[1][0] = -C4[0][0];
  C2[1][0] = (15.0 / TT(8.0)) * c20;
  C4[1][0] = -C2[1][0];
  C2[0][1] = C3[0][1] = (-7.0 / TT(16.0)) * v[0][0] * c18;
  C2[0][2] = C3[0][2] = (-7.0 / TT(16.0)) * v[0][1] * c18;
  C3[1][1] = C4[0][1] = (-7.0 / TT(16.0)) * v[1][0] * c19;
  C3[1][2] = C4[0][2] = (-7.0 / TT(16.0)) * v[1][1] * c19;
  C2[1][1] = C4[1][1] = (-7.0 / TT(16.0)) * v[2][0] * c20;
  C2[1][2] = C4[1][2] = (-7.0 / TT(16.0)) * v[2][1] * c20;
  C3[0][3] = v[0][0] * v[0][0] * c18 / 32;
  C2[0][3] = -C3[0][3];
  C3[0][4] = v[0][0] * v[0][1] * c18 / 16;
  C2[0][4] = -C3[0][4];
  C3[0][5] = v[0][1] * v[0][1] * c18 / 32;
  C2[0][5] = -C3[0][5];
  C4[0][3] = v[1][0] * v[1][0] * c19 / 32;
  C3[1][3] = -C4[0][3];
  C4[0][4] = v[1][0] * v[1][1] * c19 / 16;
  C3[1][4] = -C4[0][4];
  C4[0][5] = v[1][1] * v[1][1] * c19 / 32;
  C3[1][5] = -C4[0][5];
  C2[1][3] = v[2][0] * v[2][0] * c20 / 32;
  C4[1][3] = -C2[1][3];
  C2[1][4] = v[2][0] * v[2][1] * c20 / 16;
  C4[1][4] = -C2[1][4];
  C2[1][5] = v[2][1] * v[2][1] * c20 / 32;
  C4[1][5] = -C2[1][5];

  C5[0] = -sign[0] * (trafo(0)[1] * normal[0][0] + trafo(1)[1] * normal[0][1]);
  C5[1] =
      sign[1]
      * ((trafo(0)[0] + trafo(0)[1]) * normal[1][0] + (trafo(1)[0] + trafo(1)[1]) * normal[1][1])
      / sqrt(TT(2.0));
  C5[2] = -sign[2] * (trafo(0)[0] * normal[2][0] + trafo(1)[0] * normal[2][1]);
}

template<typename TT, int sDim, int tDim>
SymTensorT<TT, sDim>
ArgyrisBasicElementT<TT, sDim, tDim>::multiplyTheta(const SymTensorT<TT, sDim> &S) const {
  return SymTensorT<TT, sDim>(theta[0][0] * S(0, 0) + theta[0][1] * S(0, 1) + theta[0][2] * S(1, 1),
                              theta[1][0] * S(0, 0) + theta[1][1] * S(0, 1) + theta[1][2] * S(1, 1),
                              theta[2][0] * S(0, 0) + theta[2][1] * S(0, 1)
                                  + theta[2][2] * S(1, 1));
}

template<typename TT, int sDim, int tDim>
int ArgyrisBasicElementT<TT, sDim, tDim>::indexingShape(int i, int k) const {
  if (i < 3) return i * 6 + k;
  return 15 + i;
}

template<typename TT, int sDim, int tDim>
ArgyrisElementT<TT, sDim, tDim>::ArgyrisElementT(const VectorMatrixBase &base, const Cell &c) :
    ArgyrisBasicElementT<TT, sDim, tDim>(base, c), value(this->nQ(), this->shape.size()),
    gradient(this->nQ(), this->shape.size()), hessian(this->nQ(), this->shape.size()) {
  const TransformationT<TT, sDim> &trafo = this->GetTransformation(0);
  for (int q = 0; q < this->nQ(); ++q) {
    vector<VectorFieldT<TT, sDim>> G(21);
    vector<SymTensorT<TT, sDim>> H(21);
    for (int i = 0; i < 21; ++i) {
      G[i] = trafo * this->shape.LocalGradient(q, i);
      H[i] = this->multiplyTheta(this->shape.LocalHessian(q, i));
    }
    for (int i = 0; i < 21; ++i) {
      applyC(this->C0, this->C1, this->C2, this->C3, this->C4, this->C5, value[q][i],
             this->shape.values()[q], i);
      applyC(this->C0, this->C1, this->C2, this->C3, this->C4, this->C5, gradient[q][i], G, i);
      applyC(this->C0, this->C1, this->C2, this->C3, this->C4, this->C5, hessian[q][i], H, i);
    }
  }
}

template<typename TT, int sDim, int tDim>
ArgyrisElementT<TT, sDim, tDim>::ArgyrisElementT(const VectorMatrixBase &base,
                                                 const BasicElementT<TT, sDim, tDim> &baseElement) :
    ArgyrisBasicElementT<TT, sDim, tDim>(base, baseElement), value(this->nQ(), this->shape.size()),
    gradient(this->nQ(), this->shape.size()), hessian(this->nQ(), this->shape.size()) {
  const TransformationT<TT, sDim> &trafo = this->GetTransformation(0);
  for (int q = 0; q < this->nQ(); ++q) {
    vector<VectorFieldT<TT, sDim>> G(21);
    vector<SymTensorT<TT, sDim>> H(21);
    for (int i = 0; i < 21; ++i) {
      G[i] = trafo * this->shape.LocalGradient(q, i);
      H[i] = this->multiplyTheta(this->shape.LocalHessian(q, i));
    }
    for (int i = 0; i < 21; ++i) {
      applyC(this->C0, this->C1, this->C2, this->C3, this->C4, this->C5, value[q][i],
             this->shape.values()[q], i);
      applyC(this->C0, this->C1, this->C2, this->C3, this->C4, this->C5, gradient[q][i], G, i);
      applyC(this->C0, this->C1, this->C2, this->C3, this->C4, this->C5, hessian[q][i], H, i);
    }
  }
}

template<typename TT, int sDim, int tDim>
TT ArgyrisElementT<TT, sDim, tDim>::Value(const PointT<TT, sDim, tDim> &z, int i, int k) const {
  return Value(z, ShapeId{this->indexingShape(i, k)});
}

template<typename TT, int sDim, int tDim>
TT ArgyrisElementT<TT, sDim, tDim>::Value(const PointT<TT, sDim, tDim> &z,
                                          const ShapeId &iter) const {
  vector<TT> V(21);
  for (int j = 0; j < 21; ++j)
    V[j] = this->shape(z, j);
  TT U = TT(0.0);
  applyC(this->C0, this->C1, this->C2, this->C3, this->C4, this->C5, U, V, iter.id);
  return U;
}

template<typename TT, int sDim, int tDim>
TT ArgyrisElementT<TT, sDim, tDim>::Value(int q, const Vector &u, int m) const {
  TT U = TT(0.0);
  for (int i = 0, shapeId = 0; i < this->size(); ++i)
    for (int k = 0; k < this->get_maxk(i); ++k, ++shapeId)
      U += u(this->r(i), this->indexing_k(i, k, m)) * value[q][shapeId];
  return U;
}

template<typename TT, int sDim, int tDim>
TT ArgyrisElementT<TT, sDim, tDim>::Value(const PointT<TT, sDim, tDim> &z, const Vector &u,
                                          int m) const {
  vector<TT> v_z = this->values(z);
  TT U{};
  for (int i = 0, shapeId = 0; i < this->size(); ++i)
    for (int k = 0; k < this->get_maxk(i); ++k, ++shapeId)
      U += u(this->r(i), this->indexing_k(i, k, m)) * v_z[shapeId];
  return U;
}

template<typename TT, int sDim, int tDim>
VectorFieldT<TT, sDim> ArgyrisElementT<TT, sDim, tDim>::Derivative(const PointT<TT, sDim, tDim> &z,
                                                                   int i, int k) const {
  return Derivative(z, ShapeId{this->indexingShape(i, k)});
}

template<typename TT, int sDim, int tDim>
VectorFieldT<TT, sDim> ArgyrisElementT<TT, sDim, tDim>::Derivative(const PointT<TT, sDim, tDim> &z,
                                                                   const ShapeId &iter) const {
  vector<VectorFieldT<TT, sDim>> G(21);
  for (int j = 0; j < 21; ++j)
    G[j] = this->GetTransformation(0) * this->shape.LocalGradient(z, j);
  VectorFieldT<TT, sDim> D;
  applyC(this->C0, this->C1, this->C2, this->C3, this->C4, this->C5, D, G, iter.id);
  return D;
}

template<typename TT, int sDim, int tDim>
VectorFieldT<TT, sDim> ArgyrisElementT<TT, sDim, tDim>::Derivative(int q, const Vector &u,
                                                                   int m) const {
  VectorFieldT<TT, sDim> D;
  for (int i = 0, shapeId = 0; i < this->size(); ++i)
    for (int k = 0; k < this->get_maxk(i); ++k, ++shapeId)
      D += u(this->r(i), this->indexing_k(i, k, m)) * gradient[q][shapeId];
  return D;
}

template<typename TT, int sDim, int tDim>
VectorFieldT<TT, sDim> ArgyrisElementT<TT, sDim, tDim>::Derivative(const PointT<TT, sDim, tDim> &z,
                                                                   const Vector &u, int m) const {

  vector<VectorFieldT<TT, sDim>> d_z = this->derivatives(z);
  VectorFieldT<TT, sDim> D;
  for (int i = 0, shapeId = 0; i < this->size(); ++i)
    for (int k = 0; k < this->get_maxk(i); ++k, ++shapeId)
      D += u(this->r(i), this->indexing_k(i, k, m)) * d_z[shapeId];
  return D;
}

template<typename TT, int sDim, int tDim>
SymTensorT<TT, sDim> ArgyrisElementT<TT, sDim, tDim>::Hessian(const PointT<TT, sDim, tDim> &z,
                                                              int i, int k) const {
  return Hessian(z, ShapeId{this->indexingShape(i, k)});
}

template<typename TT, int sDim, int tDim>
SymTensorT<TT, sDim> ArgyrisElementT<TT, sDim, tDim>::Hessian(const PointT<TT, sDim, tDim> &z,
                                                              const ShapeId &iter) const {
  vector<SymTensorT<TT, sDim>> H(21);
  for (int j = 0; j < 21; ++j)
    H[j] = this->multiplyTheta(this->shape.LocalHessian(z, j));
  SymTensorT<TT, sDim> D;
  applyC(this->C0, this->C1, this->C2, this->C3, this->C4, this->C5, D, H, iter.id);
  return D;
}

template<typename TT, int sDim, int tDim>
SymTensorT<TT, sDim> ArgyrisElementT<TT, sDim, tDim>::Hessian(int q, const Vector &u, int m) const {
  SymTensorT<TT, sDim> H;
  for (int i = 0, shapeId = 0; i < this->size(); ++i)
    for (int k = 0; k < this->get_maxk(i); ++k, ++shapeId)
      H += TT(u(this->r(i), this->indexing_k(i, k, m))) * hessian[q][shapeId];
  return H;
}

template<typename TT, int sDim, int tDim>
SymTensorT<TT, sDim> ArgyrisElementT<TT, sDim, tDim>::Hessian(const PointT<TT, sDim, tDim> &z,
                                                              const Vector &u, int m) const {
  SymTensorT<TT, sDim> H;
  vector<SymTensorT<TT, sDim>> h_z = this->hessians(z);
  for (int i = 0, shapeId = 0; i < this->size(); ++i)
    for (int k = 0; k < this->get_maxk(i); ++k, ++shapeId)
      H += u(this->r(i), this->indexing_k(i, k, m)) * h_z[shapeId];
  return H;
}

template<typename TT, int sDim, int tDim>
TT ArgyrisElementT<TT, sDim, tDim>::Laplace(int q, int i, int k) const {
  return Laplace(q, ShapeId{this->indexingShape(i, k)});
}

template<typename TT, int sDim, int tDim>
TT ArgyrisElementT<TT, sDim, tDim>::Laplace(int q, const ShapeId &iter) const {
  return hessian[q][iter.id](0, 0) + hessian[q][iter.id](1, 1);
}

template<typename TT, int sDim, int tDim>
TT ArgyrisElementT<TT, sDim, tDim>::Laplace(const PointT<TT, sDim, tDim> &z, int i, int k) const {
  SymTensorT<TT, sDim> H = Hessian(z, i, k);
  return H(0, 0) + H(1, 1);
}

template<typename TT, int sDim, int tDim>
TT ArgyrisElementT<TT, sDim, tDim>::Laplace(const PointT<TT, sDim, tDim> &z,
                                            const ShapeId &iter) const {
  SymTensorT<TT, sDim> H = Hessian(z, iter);
  return H(0, 0) + H(1, 1);
}

template<typename TT, int sDim, int tDim>
TT ArgyrisElementT<TT, sDim, tDim>::Laplace(int q, const Vector &u, int m) const {
  TT L = TT(0.0);
  for (int i = 0, shapeId = 0; i < this->size(); ++i)
    for (int k = 0; k < this->get_maxk(i); ++k, ++shapeId)
      L += u(this->r(i), this->indexing_k(i, k, m))
           * (hessian[q][shapeId](0, 0) + hessian[q][shapeId](1, 1));
  return L;
}

template<typename TT, int sDim, int tDim>
TT ArgyrisElementT<TT, sDim, tDim>::Laplace(const PointT<TT, sDim, tDim> &z, const Vector &u,
                                            int m) const {
  SymTensorT<TT, sDim> H = Hessian(z, u, m);
  return H(0, 0) + H(1, 1);
}


template class ArgyrisBasicElementT<>;

template class ArgyrisElementT<>;

#ifdef BUILD_IA

template class ArgyrisBasicElementT<IAInterval, SpaceDimension, TimeDimension>;

template class ArgyrisElementT<IAInterval, SpaceDimension, TimeDimension>;

template<typename T>
std::vector<T> B5_TMP(const BernsteinCoeff &c, const T &sqrt2) {
  return std::vector<T>{c[0],
                        c[0] + c[1] / 5,
                        c[0] + c[2] / 5,
                        c[0] + (2 * c[1] + c[3] / 4) / 5,
                        c[0] + (c[1] + c[2] + c[4] / 4) / 5,
                        c[0] + (2 * c[2] + c[5] / 4) / 5,
                        c[6] - (2 * c[7] - c[9] / 4) / 5,
                        c[0] + 2 * c[1] / 5 - (c[2] + c[8]) / 6 + c[3] / 20 - (c[4] - c[10]) / 30
                            - 8 * c[18] / 15,
                        c[0] - (c[1] + c[13]) / 6 + 2 * c[2] / 5 - (c[4] - c[16]) / 30 + c[5] / 20
                            - 8 * c[20] / 15,
                        c[12] - (2 * c[14] - c[17] / 4) / 5,
                        c[6] - c[7] / 5,
                        c[6] - (2 * c[7] - c[8]) / 5 + (c[9] - c[10]) / 20,
                        (c[6] + c[12]) / 2 + (-7 * (c[7] + c[14]) + 17 * (c[8] + c[13])) / 60
                            + (c[9] + c[17]) / 120 - (c[10] + c[16]) / 20 + (c[11] + c[15]) / 24
                            - 4 * sqrt2 * c[19] / 15,
                        c[12] + (c[13] - 2 * c[14]) / 5 - (c[16] - c[17]) / 20,
                        c[12] - c[14] / 5,
                        c[6],
                        c[6] - (c[7] - c[8]) / 5,
                        c[6] - 2 * (c[7] - c[8]) / 5 + ((c[9] + c[11]) / 2 - c[10]) / 10,
                        c[12] + 2 * (c[13] - c[14]) / 5 + ((c[15] + c[17]) / 2 - c[16]) / 10,
                        c[12] + (c[13] - c[14]) / 5,
                        c[12]};
}

template<typename T>
std::vector<T> B4_TMP_x(const BernsteinCoeff &c, const T &sqrt2) {
  return std::vector<T>{c[1],
                        c[1] + c[3] / 4,
                        c[1] + c[4] / 4,
                        5 * (c[6] - c[0]) - 2 * (c[1] + c[7]) - (c[3] - c[9]) / 4,
                        c[1] + c[3] / 4
                            + ((-5 * (c[4] / 2 + c[8]) - 11 * c[2] + c[10]) / 2 - 8 * c[18]) / 3,
                        ((-5 * (c[1] + c[13]) - c[4] + c[16]) / 2 - 8 * c[20]) / 3,
                        c[7] - c[9] / 4,
                        5 * (c[6] - c[0]) - 2 * (c[1] + c[7])
                            + ((5 * (c[2] - c[10] / 2) + c[4] + 11 * c[8]) / 2 + 8 * c[18]) / 3
                            - (c[3] - c[9]) / 4,
                        -5 * c[0] - 2 * c[2]
                            + ((9 * c[13] - c[5] - c[10]) / 2 + 5 * (c[6] + c[12])) / 2
                            + ((((c[9] + 5 * c[11] + 5 * c[15] + c[17]) / 2
                                 + (-7 * (c[7] + c[14]) + 17 * c[8] - 5 * c[16]))
                                    / 2
                                + (5 * c[1] + c[4]))
                                   / 2
                               - 4 * (c[19] * sqrt2 - 2 * c[20]))
                                  / 3,
                        c[13] - c[16] / 4,
                        c[7],
                        c[7] - (c[9] - c[10]) / 4,
                        ((c[16] - c[10]) / 2 + 5 * (c[6] - c[12])) / 2
                            + (((5 * (c[9] - c[15]) + c[11] - c[17]) / 2
                                + (-17 * (c[7] + c[13]) + 7 * (c[8] + c[14])))
                                   / 4
                               + 4 * c[19] * sqrt2)
                                  / 3,
                        c[13] + (c[15] - c[16]) / 4,
                        c[13]};
}

template<typename T>
std::vector<T> B4_TMP_y(const BernsteinCoeff &c, const T &sqrt2) {
  return std::vector<T>{c[2],
                        c[2] + c[4] / 4,
                        c[2] + c[5] / 4,
                        ((-5 * (c[2] + c[8]) - c[4] + c[10]) / 2 - 8 * c[18]) / 3,
                        ((-5 * (c[4] / 2 + c[13]) - 11 * c[1] + c[16]) / 2 - 8 * c[20]) / 3 + c[2]
                            + c[5] / 4,
                        5 * (c[12] - c[0]) - 2 * (c[2] + c[14]) - (c[5] - c[17]) / 4,
                        c[8] - c[10] / 4,
                        -5 * c[0] - 2 * c[1]
                            + ((9 * c[8] - c[3] - c[16]) / 2 + 5 * (c[6] + c[12])) / 2
                            + ((((c[9] + 5 * (c[11] + c[15]) + c[17]) / 2
                                 + (17 * c[13] - 7 * (c[7] + c[14]) - 5 * c[10]))
                                    / 2
                                + (5 * c[2] + c[4]))
                                   / 2
                               + 4 * (2 * c[18] - c[19] * sqrt2))
                                  / 3,
                        5 * (c[12] - c[0]) - 2 * (c[2] + c[14]) + (c[17] - c[5]) / 4
                            + ((-5 * c[16] / 2 + (5 * c[1] + c[4] + 11 * c[13])) / 2 + 8 * c[20])
                                  / 3,
                        c[14] - c[17] / 4,
                        c[8],
                        c[8] - (c[10] - c[11]) / 4,
                        ((c[10] - c[16]) / 2 + 5 * (c[12] - c[6])) / 2
                            + (((c[15] - c[9] + 5 * (c[17] - c[11])) / 2
                                + (7 * (c[7] + c[13]) - 17 * (c[8] + c[14])))
                                   / 4
                               + 4 * c[19] * sqrt2)
                                  / 3,
                        c[14] + (c[16] - c[17]) / 4,
                        c[14]};
}

template<typename T>
std::vector<T> B3_TMP_xx(const BernsteinCoeff &c, const T &sqrt2) {
  return std::vector<T>{c[3],
                        2 * (10 * (c[6] - c[0]) - 6 * c[1] - c[3] - 4 * c[7]) + c[9],
                        (-22 * c[2] - 8 * c[4] - 10 * c[8] + 2 * c[10] - 32 * c[18]) / 3 + c[3],
                        2 * (10 * (c[0] - c[6]) + 4 * c[1] + 6 * c[7] - c[9]) + c[3],
                        2 * (10 * (c[6] - c[0]) - 6 * c[1] - c[3] - 4 * c[7]) + c[9]
                            + (32 * (c[2] + c[8] + 2 * c[18]) + 7 * (c[4] - c[10])) / 3,
                        10 * (c[6] + c[12]) - 20 * c[0] - 8 * c[2] - c[5] - c[10]
                            + ((c[9] + 5 * (c[11] + c[15]) + c[17]) / 2
                               + (20 * c[1] + 4 * c[4] + 17 * c[8] + 37 * c[13]
                                  - 7 * (c[7] + c[14] + c[16]) - 16 * c[19] * sqrt2 + 64 * c[20]))
                                  / 3,
                        c[9],
                        2 * (10 * (c[0] - c[6]) + 4 * c[1] + 6 * c[7] - c[9]) + c[3]
                            + 2 * (4 * c[10] - 5 * c[2] - c[4] - 11 * c[8] - 16 * c[18]) / 3,
                        2 * (10 * (c[0] - c[12]) + 4 * c[2]) + c[5]
                            + (2
                                   * (c[9] - c[4] - c[11] - 5 * (c[1] + c[7] + c[8]) - 22 * c[13]
                                      + 7 * c[14] + 4 * c[16])
                               - 5 * c[15] - c[17] + 32 * (c[19] * sqrt2 - c[20]))
                                  / 3,
                        c[15]};
}

template<typename T>
std::vector<T> B3_TMP_xy(const BernsteinCoeff &c, const T &sqrt2) {
  return std::vector<T>{c[4],
                        (-22 * c[2] - 5 * c[4] - 10 * c[8] + 2 * c[10] - 32 * c[18]) / 3,
                        (-22 * c[1] - 5 * c[4] - 10 * c[13] + 2 * c[16] - 32 * c[20]) / 3,
                        (10 * c[2] + 2 * c[4] + 22 * c[8] - 5 * c[10] + 32 * c[18]) / 3,
                        10 * (c[12] - 2 * c[0] + c[6]) + 9 * (c[8] + c[13]) - c[3] - c[5]
                            + ((c[9] + 5 * (c[11] + c[15]) + c[17]) / 2
                               + (-2 * (c[1] + c[2]) + 7 * (c[4] - c[7] - c[14])
                                  - 5 * (c[10] + c[16])
                                  + 16 * (2 * (c[18] + c[20]) - c[19] * sqrt2)))
                                  / 3,
                        (10 * c[1] + 2 * c[4] + 22 * c[13] - 5 * c[16] + 32 * c[20]) / 3,
                        c[10],
                        10 * (2 * c[0] - c[6] - c[12]) + 8 * c[1] - 5 * c[8] + c[3] + c[16]
                            + ((c[11] - c[9] - 5 * c[15] - c[17]) / 2
                               + (-10 * c[2] + 2 * (c[10] - c[4]) + 7 * (c[7] + c[14]) - 17 * c[13]
                                  + 16 * (c[19] * sqrt2 - 2 * c[18])))
                                  / 3,
                        10 * (2 * c[0] - c[6] - c[12]) + 8 * c[2] - 5 * c[13] + c[5] + c[10]
                            + ((c[15] - c[9] - 5 * c[11] - c[17]) / 2
                               + (-10 * c[1] + 2 * (c[16] - c[4]) + 7 * (c[7] + c[14]) - 17 * c[8]
                                  + 16 * (c[19] * sqrt2 - 2 * c[20])))
                                  / 3,
                        c[16]};
}

template<typename T>
std::vector<T> B3_TMP_yy(const BernsteinCoeff &c, const T &sqrt2) {
  return std::vector<T>{c[5],
                        c[5] + 2 * (-11 * c[1] - 4 * c[4] - 5 * c[13] + c[16] - 16 * c[20]) / 3,
                        2 * (10 * (c[12] - c[0]) - 6 * c[2] - 4 * c[14] - c[5]) + c[17],
                        2 * (-10 * c[0] - 4 * c[1] + 5 * (c[6] + c[12]))
                            + (20 * c[2] + 4 * c[4] - 7 * (c[7] + c[10] + c[14]) + 37 * c[8]
                               + 17 * c[13] + 64 * c[18] - 16 * c[19] * sqrt2
                               + (c[9] + 5 * (c[11] + c[15]) + c[17]) / 2)
                                  / 3
                            - c[3] - c[16],
                        2 * (10 * (c[12] - c[0]) - 6 * c[2] - c[5] - 4 * c[14])
                            + (7 * (c[4] - c[16]) + 32 * (c[1] + c[13] + 2 * c[20])) / 3 + c[17],
                        2 * (10 * (c[0] - c[12]) + 4 * c[2] + 6 * c[14] - c[17]) + c[5],
                        c[11],
                        4 * (5 * (c[0] - c[6]) + 2 * c[1])
                            + (-10 * (c[2] + c[13] + c[14]) + 14 * c[7] - 44 * c[8] - c[9]
                               + 8 * c[10] - 5 * c[11] - 2 * (c[4] + c[15] - c[17])
                               - 32 * (c[18] - c[19] * sqrt2))
                                  / 3
                            + c[3],
                        2
                                * (10 * (c[0] - c[12]) + 4 * c[2] + 6 * c[14] - c[17]
                                   + (-5 * c[1] - c[4] - 11 * c[13] + 4 * c[16] - 16 * c[20]) / 3)
                            + c[5],
                        c[17]};
}

template<typename TT, int sDim, int tDim>
inline TT vectorValue(const ArgyrisBasicElementT<TT, sDim, tDim> &element, const Vector &vec, int i,
                      int k) {
  return TT(vec(element[i], k));
}

template<typename TT, int sDim, int tDim>
BernsteinRangesT<TT, sDim, tDim>::BernsteinRangesT(
    const ArgyrisBasicElementT<TT, sDim, tDim> &element, std::vector<BernsteinDataT<TT>> data) {
  IAInterval sqrt2 = IAInterval::Sqrt2();
  BernsteinCoeff c(21);
  for (int n = 0; n < data.size(); ++n) {
    for (int i = 0; i < 3; ++i) {
      c[6 * i] += data[n].lambda * vectorValue(element, data[n].vec, i, 0);
      c[6 * i + 1] += data[n].lambda
                      * (element.C0[0][0] * vectorValue(element, data[n].vec, i, 1)
                         + element.C0[0][1] * vectorValue(element, data[n].vec, i, 2));
      c[6 * i + 2] += data[n].lambda
                      * (element.C0[1][0] * vectorValue(element, data[n].vec, i, 1)
                         + element.C0[1][1] * vectorValue(element, data[n].vec, i, 2));
      c[6 * i + 3] += data[n].lambda
                      * (element.C1[0][0] * vectorValue(element, data[n].vec, i, 3)
                         + element.C1[0][1] * vectorValue(element, data[n].vec, i, 4)
                         + element.C1[0][2] * vectorValue(element, data[n].vec, i, 5));
      c[6 * i + 4] += data[n].lambda
                      * (element.C1[1][0] * vectorValue(element, data[n].vec, i, 3)
                         + element.C1[1][1] * vectorValue(element, data[n].vec, i, 4)
                         + element.C1[1][2] * vectorValue(element, data[n].vec, i, 5));
      c[6 * i + 5] += data[n].lambda
                      * (element.C1[2][0] * vectorValue(element, data[n].vec, i, 3)
                         + element.C1[2][1] * vectorValue(element, data[n].vec, i, 4)
                         + element.C1[2][2] * vectorValue(element, data[n].vec, i, 5));
    }
    std::vector<TT> tmp(3);
    for (int k = 0; k < 6; ++k) {
      tmp[0] += element.C2[0][k] * vectorValue(element, data[n].vec, 0, k)
                + element.C3[0][k] * vectorValue(element, data[n].vec, 1, k);
      tmp[1] += element.C3[1][k] * vectorValue(element, data[n].vec, 1, k)
                + element.C4[0][k] * vectorValue(element, data[n].vec, 2, k);
      tmp[2] += element.C2[1][k] * vectorValue(element, data[n].vec, 0, k)
                + element.C4[1][k] * vectorValue(element, data[n].vec, 2, k);
    }
    for (int i = 0; i < 3; ++i) {
      tmp[i] += element.C5[i] * vectorValue(element, data[n].vec, i + 3, 0);
      c[18 + i] += data[n].lambda * tmp[i];
    }
  }

  coeff = B5_TMP(c, sqrt2);
  const BernsteinCoeff &tmp_x = B4_TMP_x(c, sqrt2);
  const BernsteinCoeff &tmp_y = B4_TMP_y(c, sqrt2);
  const BernsteinCoeff &tmp_xx = B3_TMP_xx(c, sqrt2);
  const BernsteinCoeff &tmp_xy = B3_TMP_xy(c, sqrt2);
  const BernsteinCoeff &tmp_yy = B3_TMP_yy(c, sqrt2);

  const TransformationT<TT, sDim> &trafo = element.GetTransformation(0);
  for (int i = 0; i < tmp_x.size(); ++i) {
    coeffDx.push_back(trafo[0][0] * tmp_x[i] + trafo[0][1] * tmp_y[i]);
    coeffDy.push_back(trafo[1][0] * tmp_x[i] + trafo[1][1] * tmp_y[i]);
  }
  for (int i = 0; i < tmp_xx.size(); ++i) {
    coeffDxx.push_back(element.theta[0][0] * tmp_xx[i] + element.theta[0][1] * tmp_xy[i]
                       + element.theta[0][2] * tmp_yy[i]);
    coeffDxy.push_back(element.theta[1][0] * tmp_xx[i] + element.theta[1][1] * tmp_xy[i]
                       + element.theta[1][2] * tmp_yy[i]);
    coeffDyy.push_back(element.theta[2][0] * tmp_xx[i] + element.theta[2][1] * tmp_xy[i]
                       + element.theta[2][2] * tmp_yy[i]);
  }
}

template<typename TT, int sDim, int tDim>
BernsteinRangesT<TT, sDim, tDim>::BernsteinRangesT(
    const ArgyrisBasicElementT<TT, sDim, tDim> &element, const Vector &vec) :
    BernsteinRangesT(element, {{TT(1.0), vec}}) {}

template<typename TT, int sDim, int tDim>
void BernsteinRangesT<TT, sDim, tDim>::Range(IAInterval &range) const {
  range = coeff[0];
  for (int j = 1; j < coeff.size(); ++j) {
    range = range | coeff[j];
  }
}

template<typename TT, int sDim, int tDim>
void BernsteinRangesT<TT, sDim, tDim>::RangeGradient(IAInterval &rangeDx,
                                                     IAInterval &rangeDy) const {
  rangeDx = coeffDx[0];
  rangeDy = coeffDy[0];
  for (int j = 1; j < coeffDx.size(); ++j) {
    rangeDx = rangeDx | coeffDx[j];
    rangeDy = rangeDy | coeffDy[j];
  }
}

template<typename TT, int sDim, int tDim>
void BernsteinRangesT<TT, sDim, tDim>::RangeHessian(IAInterval &rangeDxx, IAInterval &rangeDxy,
                                                    IAInterval &rangeDyy) const {
  rangeDxx = coeffDxx[0];
  rangeDxy = coeffDxy[0];
  rangeDyy = coeffDyy[0];
  for (int j = 1; j < coeffDx.size(); ++j) {
    rangeDxx = rangeDxx | coeffDxx[j];
    rangeDxy = rangeDxy | coeffDxy[j];
    rangeDyy = rangeDyy | coeffDyy[j];
  }
}

template<typename TT, int sDim, int tDim>
void BernsteinRangesT<TT, sDim, tDim>::AbsMax(double &absMax) const {
  IAInterval range;
  Range(range);
  absMax = sup(abs(range));
}

template<typename TT, int sDim, int tDim>
void BernsteinRangesT<TT, sDim, tDim>::AbsMaxGradient(double &absMaxDx, double &absMaxDy) const {
  IAInterval rangeDx, rangeDy;
  RangeGradient(rangeDx, rangeDy);
  absMaxDx = sup(abs(rangeDx));
  absMaxDy = sup(abs(rangeDy));
}

template<typename TT, int sDim, int tDim>
void BernsteinRangesT<TT, sDim, tDim>::AbsMaxHessian(double &absMaxDxx, double &absMaxDxy,
                                                     double &absMaxDyy) const {
  IAInterval rangeDxx, rangeDxy, rangeDyy;
  RangeHessian(rangeDxx, rangeDxy, rangeDyy);
  absMaxDxx = sup(abs(rangeDxx));
  absMaxDxy = sup(abs(rangeDxy));
  absMaxDyy = sup(abs(rangeDyy));
}

template class BernsteinRangesT<>;

template class BernsteinRangesT<IAInterval, SpaceDimension, TimeDimension>;

#endif