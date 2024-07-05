#include "DivergenceFreeElement.hpp"

template<typename TT, int sDim, int tDim>
DivergenceFreeElementT<TT, sDim, tDim>::DivergenceFreeElementT(const VectorMatrixBase &base,
                                                               const Cell &c) :
    ArgyrisBasicElementT<TT, sDim, tDim>(base, c), gradient(this->nQ(), this->shape.size()),
    hessian(this->nQ(), this->shape.size()) {
  const TransformationT<TT, sDim> &trafo = this->GetTransformation(0);
  for (int q = 0; q < this->nQ(); ++q) {
    vector<VectorFieldT<TT, sDim>> G(21);
    vector<SymTensorT<TT, sDim>> H(21);
    for (int i = 0; i < 21; ++i) {
      G[i] = trafo * this->shape.LocalGradient(q, i);
      H[i] = this->multiplyTheta(this->shape.LocalHessian(q, i));
    }
    for (int i = 0; i < 21; ++i) {
      applyC(this->C0, this->C1, this->C2, this->C3, this->C4, this->C5, gradient[q][i], G, i);
      applyC(this->C0, this->C1, this->C2, this->C3, this->C4, this->C5, hessian[q][i], H, i);
    }
  }
}

template<typename TT, int sDim, int tDim>
DivergenceFreeElementT<TT, sDim, tDim>::DivergenceFreeElementT(
    const VectorMatrixBase &base, const BasicElementT<TT, sDim, tDim> &baseElement) :
    ArgyrisBasicElementT<TT, sDim, tDim>(base, baseElement),
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
      applyC(this->C0, this->C1, this->C2, this->C3, this->C4, this->C5, gradient[q][i], G, i);
      applyC(this->C0, this->C1, this->C2, this->C3, this->C4, this->C5, hessian[q][i], H, i);
    }
  }
}

template<typename TT, int sDim, int tDim>
VelocityT<TT, sDim> DivergenceFreeElementT<TT, sDim, tDim>::VelocityField(int q, int i,
                                                                          int k) const {
  return VelocityField(q, ShapeId{this->indexingShape(i, k)});
}

template<typename TT, int sDim, int tDim>
VelocityT<TT, sDim>
DivergenceFreeElementT<TT, sDim, tDim>::VelocityField(int q, const ShapeId &iter) const {
  return VelocityT<TT, sDim>(-gradient[q][iter.id][1], gradient[q][iter.id][0]);
}

template<typename TT, int sDim, int tDim>
VelocityT<TT, sDim>
DivergenceFreeElementT<TT, sDim, tDim>::VelocityField(const PointT<TT, sDim, tDim> &z, int i,
                                                      int k) const {
  return VelocityField(z, ShapeId(this->indexingShape(i, k)));
}

template<typename TT, int sDim, int tDim>
VelocityT<TT, sDim>
DivergenceFreeElementT<TT, sDim, tDim>::VelocityField(const PointT<TT, sDim, tDim> &z,
                                                      const ShapeId &iter) const {
  vector<VectorFieldT<TT, sDim>> G(21);
  for (int j = 0; j < 21; ++j)
    G[j] = this->GetTransformation(0) * this->shape.LocalGradient(z, j);
  VectorFieldT<TT, sDim> D;
  applyC(this->C0, this->C1, this->C2, this->C3, this->C4, this->C5, D, G, iter.id);
  return VelocityT<TT, sDim>(-D[1], D[0]);
}

template<typename TT, int sDim, int tDim>
VelocityT<TT, sDim> DivergenceFreeElementT<TT, sDim, tDim>::VelocityField(int q, const Vector &u,
                                                                          int m) const {
  VectorFieldT<TT, sDim> D;
  for (int i = 0, shapeId = 0; i < this->size(); ++i)
    for (int k = 0; k < this->get_maxk(i); ++k, ++shapeId)
      D += u(this->r(i), this->indexing_k(i, k, m)) * VelocityField(q, ShapeId{shapeId});
  return D;
}

template<typename TT, int sDim, int tDim>
VelocityT<TT, sDim>
DivergenceFreeElementT<TT, sDim, tDim>::VelocityField(const PointT<TT, sDim, tDim> &z,
                                                      const Vector &u, int m) const {
  std::vector<VectorFieldT<TT, sDim>> d_z = this->derivatives(z);
  VectorFieldT<TT, sDim> D;
  for (int i = 0, shapeId = 0; i < this->size(); ++i)
    for (int k = 0; k < this->get_maxk(i); ++k, ++shapeId)
      D += u(this->r(i), this->indexing_k(i, k, m))
           * VelocityT<TT, sDim>(-d_z[shapeId][1], d_z[shapeId][0]);
  return D;
}

template<typename TT, int sDim, int tDim>
VelocityGradientT<TT, sDim>
DivergenceFreeElementT<TT, sDim, tDim>::VelocityFieldGradient(int q, int i, int k) const {
  return VelocityFieldGradient(q, this->indexingShape(i, k));
}

template<typename TT, int sDim, int tDim>
VelocityGradientT<TT, sDim>
DivergenceFreeElementT<TT, sDim, tDim>::VelocityFieldGradient(int q, const ShapeId &iter) const {
  return VelocityGradientT<TT, sDim>(-hessian[q][iter.id](1, 0), -hessian[q][iter.id](1, 1),
                                     hessian[q][iter.id](0, 0), hessian[q][iter.id](0, 1));
}

template<typename TT, int sDim, int tDim>
VelocityGradientT<TT, sDim>
DivergenceFreeElementT<TT, sDim, tDim>::VelocityFieldGradient(const PointT<TT, sDim, tDim> &z,
                                                              int i, int k) const {
  return VelocityFieldGradient(z, ShapeId{this->indexingShape(i, k)});
}

template<typename TT, int sDim, int tDim>
VelocityGradientT<TT, sDim>
DivergenceFreeElementT<TT, sDim, tDim>::VelocityFieldGradient(const PointT<TT, sDim, tDim> &z,
                                                              const ShapeId &iter) const {
  vector<SymTensorT<TT, sDim>> H(21);
  for (int j = 0; j < 21; ++j)
    H[j] = this->multiplyTheta(this->shape.LocalHessian(z, j));
  SymTensorT<TT, sDim> D;
  applyC(this->C0, this->C1, this->C2, this->C3, this->C4, this->C5, D, H, iter.id);

  return VelocityGradientT<TT, sDim>(-D(1, 0), -D(1, 1), D(0, 0), D(0, 1));
}

template<typename TT, int sDim, int tDim>
VelocityGradientT<TT, sDim>
DivergenceFreeElementT<TT, sDim, tDim>::VelocityFieldGradient(int q, const Vector &u, int m) const {
  VelocityGradientT<TT, sDim> H;
  for (int i = 0, shapeId = 0; i < this->size(); ++i)
    for (int k = 0; k < this->get_maxk(i); ++k, ++shapeId)
      H += u(this->r(i), this->indexing_k(i, k, m)) * VelocityFieldGradient(q, ShapeId{shapeId});
  return H;
}

template<typename TT, int sDim, int tDim>
VelocityGradientT<TT, sDim>
DivergenceFreeElementT<TT, sDim, tDim>::VelocityFieldGradient(const PointT<TT, sDim, tDim> &z,
                                                              const Vector &u, int m) const {
  vector<SymTensorT<TT, sDim>> h_z = this->hessians(z);
  VelocityGradientT<TT, sDim> H;
  for (int i = 0, shapeId = 0; i < this->size(); ++i)
    for (int k = 0; k < this->get_maxk(i); ++k, ++shapeId) {
      const SymTensorT<TT, sDim> &h_z_ik = h_z[shapeId];
      H += u(this->r(i), this->indexing_k(i, k, m))
           * VelocityGradientT<TT, sDim>(-h_z_ik(1, 0), -h_z_ik(1, 1), h_z_ik(0, 0), h_z_ik(0, 1));
    }
  return H;
}

template class DivergenceFreeElementT<>;

#ifdef BUILD_IA

template class DivergenceFreeElementT<IAInterval, SpaceDimension, TimeDimension>;

#endif