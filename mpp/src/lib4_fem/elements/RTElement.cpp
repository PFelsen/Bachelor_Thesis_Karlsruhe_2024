#include "RTElement.hpp"

template<typename TT, int sDim, int tDim>
void RTElementT<TT, sDim, tDim>::init(const VectorMatrixBase &base) {
  nodalPointsOnFacesTotal = 0;
  for (int face = 0; face < this->c.Faces(); ++face) {
    TransformationT<TT, sDim> T = this->GetTransformation((*this)[face]());
    // TODO: interval arithmetic normal
    normal[face] = T * PointT<TT, sDim, tDim>(this->c.LocalFaceNormal(face));
    normal[face] /= norm(normal[face]);
    if (mid(normal[face]) > Origin) sign[face] = 1;
    else sign[face] = -1;
    // TODO: Check if this is correct for Hybrid
    nodalPointsOnFacesTotal += base.GetDoF().NumberOfStoragePointsOnFace(this->c, face);
  }
}

template<typename TT, int sDim, int tDim>
void RTElementT<TT, sDim, tDim>::initValues(const VectorMatrixBase &base) {
  vShapeValues.resize(this->nQ(), 2 * this->size() - nodalPointsOnFacesTotal);
  divShapeValues.resize(this->nQ(), 2 * this->size() - nodalPointsOnFacesTotal);
  for (int q = 0; q < this->nQ(); ++q) {
    for (int i = 0; i < nodalPointsOnFacesTotal; ++i) {
      vShapeValues[q][i] = sign[i % this->c.Faces()] * (1.0 / this->GetTransformationDet(q))
                           * (this->GetTransformation(q).ApplyJ(this->shape.LocalVector(q, i)));
      divShapeValues[q][i] = sign[i % this->c.Faces()] * (1.0 / this->GetTransformationDet(q))
                             * this->shape.LocalDiv(q, i);
    }
    for (int i = nodalPointsOnFacesTotal; i < vShapeValues[q].size(); ++i) {
      vShapeValues[q][i] = (1.0 / this->GetTransformationDet(q))
                           * (this->GetTransformation(q).ApplyJ(this->shape.LocalVector(q, i)));
      divShapeValues[q][i] = (1.0 / this->GetTransformationDet(q)) * this->shape.LocalDiv(q, i);
    }
  }
}

template<typename TT, int sDim, int tDim>
RTElementT<TT, sDim, tDim>::RTElementT(const VectorMatrixBase &base, const Cell &c, int face,
                                       int n) :
    ElementT<TT, sDim, tDim>(base, c, face, n), dim(c.dim()), normal(c.Faces()), sign(c.Faces()) {
  init(base);
}

template<typename TT, int sDim, int tDim>
RTElementT<TT, sDim, tDim>::RTElementT(const VectorMatrixBase &base, const Cell &c) :
    ElementT<TT, sDim, tDim>(base, c), dim(c.dim()), normal(c.Faces()), sign(c.Faces()) {
  init(base);
  initValues(base);
}

template<typename TT, int sDim, int tDim>
RTElementT<TT, sDim, tDim>::RTElementT(const VectorMatrixBase &base,
                                       const BasicElementT<TT, sDim, tDim> &baseElement) :
    ElementT<TT, sDim, tDim>(base, baseElement), dim(this->c.dim()), normal(this->c.Faces()),
    sign(this->c.Faces()) {
  init(base);
  initValues(base);
}

template<typename TT, int sDim, int tDim>
VectorFieldT<TT, sDim> RTElementT<TT, sDim, tDim>::VelocityField(const PointT<TT, sDim, tDim> &z,
                                                                 int i, int k) const {
  return VelocityField(z, ShapeId{indexingShape(i, k)});
}

template<typename TT, int sDim, int tDim>
VectorFieldT<TT, sDim> RTElementT<TT, sDim, tDim>::VelocityField(const PointT<TT, sDim, tDim> &z,
                                                                 const ShapeId &iter) const {
  TransformationT<TT, sDim> trafo = this->GetTransformation(z);
  if (iter.id < nodalPointsOnFacesTotal) {
    return sign[iter.id % this->c.Faces()] * (1.0 / trafo.Det())
           * (trafo.ApplyJ(this->shape.LocalVector(z, iter.id)));
  }
  return (1.0 / trafo.Det()) * (trafo.ApplyJ(this->shape.LocalVector(z, iter.id)));
}

template<typename TT, int sDim, int tDim>
VectorFieldT<TT, sDim> RTElementT<TT, sDim, tDim>::VelocityField(int q, const Vector &u) const {
  VectorFieldT<TT, sDim> V;
  for (int i = 0, shapeId = 0; i < this->size(); ++i)
    for (int k = 0; k < get_maxk(i); ++k, ++shapeId)
      V += u(this->r(i), k) * vShapeValues[q][shapeId];
  return V;
}

template<typename TT, int sDim, int tDim>
VectorFieldT<TT, sDim> RTElementT<TT, sDim, tDim>::VelocityField(int q, const Vector &u,
                                                                 int m) const {
  VectorFieldT<TT, sDim> V;
  for (int i = 0, shapeId = 0; i < this->size(); ++i)
    for (int k = 0; k < get_maxk(i); ++k, ++shapeId)
      V += u(this->r(i), indexing_k(i, k, m)) * vShapeValues[q][shapeId];
  return V;
}

template<typename TT, int sDim, int tDim>
VectorFieldT<TT, sDim> RTElementT<TT, sDim, tDim>::VelocityField(const PointT<TT, sDim, tDim> &z,
                                                                 const Vector &u) const {
  std::vector<VectorFieldT<TT, sDim>> values(2 * this->size() - nodalPointsOnFacesTotal);
  TransformationT<TT, sDim> trafo = this->GetTransformation(z);
  for (int i = 0; i < nodalPointsOnFacesTotal; ++i) {
    values[i] = sign[i % this->c.Faces()] * (1.0 / trafo.Det())
                * (trafo.ApplyJ(this->shape.LocalVector(z, i)));
  }
  for (int i = nodalPointsOnFacesTotal; i < values.size(); ++i) {
    values[i] = (1.0 / trafo.Det()) * (trafo.ApplyJ(this->shape.LocalVector(z, i)));
  }
  VectorFieldT<TT, sDim> V;
  for (int i = 0, shapeId = 0; i < this->size(); ++i)
    for (int k = 0; k < get_maxk(i); ++k, ++shapeId)
      V += u(this->r(i), k) * values[shapeId];
  return V;
}

template<typename TT, int sDim, int tDim>
VectorFieldT<TT, sDim> RTElementT<TT, sDim, tDim>::VelocityField(const PointT<TT, sDim, tDim> &z,
                                                                 const Vector &u, int m) const {
  std::vector<VectorFieldT<TT, sDim>> values(2 * this->size() - nodalPointsOnFacesTotal);
  TransformationT<TT, sDim> trafo = this->GetTransformation(z);
  for (int i = 0; i < nodalPointsOnFacesTotal; ++i) {
    values[i] = sign[i % this->c.Faces()] * (1.0 / trafo.Det())
                * (trafo.ApplyJ(this->shape.LocalVector(z, i)));
  }
  for (int i = nodalPointsOnFacesTotal; i < values.size(); ++i) {
    values[i] = (1.0 / trafo.Det()) * (trafo.ApplyJ(this->shape.LocalVector(z, i)));
  }
  VectorFieldT<TT, sDim> V;
  for (int i = 0, shapeId = 0; i < this->size(); ++i)
    for (int k = 0; k < get_maxk(i); ++k, ++shapeId)
      V += u(this->r(i), indexing_k(i, k, m)) * values[shapeId];
  return V;
}

template<typename TT, int sDim, int tDim>
TT RTElementT<TT, sDim, tDim>::VelocityFieldDivergence(const PointT<TT, sDim, tDim> &z, int i,
                                                       int k) const {
  return VelocityFieldDivergence(z, ShapeId{indexingShape(i, k)});
}

template<typename TT, int sDim, int tDim>
TT RTElementT<TT, sDim, tDim>::VelocityFieldDivergence(const PointT<TT, sDim, tDim> &z,
                                                       const ShapeId &iter) const {
  TransformationT<TT, sDim> trafo = this->GetTransformation(z);
  if (iter.id < nodalPointsOnFacesTotal) {
    return sign[iter.id % this->c.Faces()] * (1.0 / trafo.Det()) * this->shape.LocalDiv(z, iter.id);
  }
  return (1.0 / trafo.Det()) * this->shape.LocalDiv(z, iter.id);
}

template<typename TT, int sDim, int tDim>
TT RTElementT<TT, sDim, tDim>::VelocityFieldDivergence(int q, const Vector &u) const {
  TT d = 0.0;
  for (int i = 0, shapeId = 0; i < this->size(); ++i)
    for (int k = 0; k < get_maxk(i); ++k, ++shapeId)
      d += u(this->r(i), k) * divShapeValues[q][shapeId];
  return d;
}

template<typename TT, int sDim, int tDim>
TT RTElementT<TT, sDim, tDim>::VelocityFieldDivergence(int q, const Vector &u, int m) const {
  TT d = 0.0;
  for (int i = 0, shapeId = 0; i < this->size(); ++i)
    for (int k = 0; k < get_maxk(i); ++k, ++shapeId)
      d += u(this->r(i), indexing_k(i, k, m)) * divShapeValues[q][shapeId];
  return d;
}

template<typename TT, int sDim, int tDim>
TT RTElementT<TT, sDim, tDim>::VelocityFieldDivergence(const PointT<TT, sDim, tDim> &z,
                                                       const Vector &u) const {
  std::vector<TT> values(2 * this->size() - nodalPointsOnFacesTotal);
  TransformationT<TT, sDim> trafo = this->GetTransformation(z);
  for (int i = 0; i < nodalPointsOnFacesTotal; ++i) {
    values[i] = sign[i % this->c.Faces()] * (1.0 / trafo.Det()) * this->shape.LocalDiv(z, i);
  }
  for (int i = nodalPointsOnFacesTotal; i < values.size(); ++i) {
    values[i] = (1.0 / trafo.Det()) * this->shape.LocalDiv(z, i);
  }
  TT d = 0.0;
  for (int i = 0, shapeId = 0; i < this->size(); ++i)
    for (int k = 0; k < get_maxk(i); ++k, ++shapeId)
      d += u(this->r(i), k) * values[shapeId];
  return d;
}

template<typename TT, int sDim, int tDim>
TT RTElementT<TT, sDim, tDim>::VelocityFieldDivergence(const PointT<TT, sDim, tDim> &z,
                                                       const Vector &u, int m) const {
  std::vector<TT> values(2 * this->size() - nodalPointsOnFacesTotal);
  TransformationT<TT, sDim> trafo = this->GetTransformation(z);
  for (int i = 0; i < nodalPointsOnFacesTotal; ++i) {
    values[i] = sign[i % this->c.Faces()] * (1.0 / trafo.Det()) * this->shape.LocalDiv(z, i);
  }
  for (int i = nodalPointsOnFacesTotal; i < values.size(); ++i) {
    values[i] = (1.0 / trafo.Det()) * this->shape.LocalDiv(z, i);
  }
  TT d = 0.0;
  for (int i = 0, shapeId = 0; i < this->size(); ++i)
    for (int k = 0; k < get_maxk(i); ++k, ++shapeId)
      d += u(this->r(i), indexing_k(i, k, m)) * values[shapeId];
  return d;
}

template class RTElementT<>;

#ifdef BUILD_IA

template class RTElementT<IAInterval, SpaceDimension, TimeDimension>;

#endif

template<typename TT, int sDim, int tDim>
RTFaceElementT<TT, sDim, tDim>::RTFaceElementT(const VectorMatrixBase &base, const Cell &c,
                                               int face, int n) :
    RTElementT<TT, sDim, tDim>(base, c, face, n) {

  this->vShapeValues.resize(this->nQ(), 2 * this->size() - this->nodalPointsOnFacesTotal);
  this->divShapeValues.resize(this->nQ(), 2 * this->size() - this->nodalPointsOnFacesTotal);
  for (int q = 0; q < this->nQ(); ++q) {
    for (int i = 0; i < this->nodalPointsOnFacesTotal; ++i) {
      this->vShapeValues[q][i] =
          this->sign[i % this->c.Faces()] * (1.0 / this->GetTransformationDet(q))
          * (this->GetTransformation(q).ApplyJ(this->shape.LocalVector(this->QLocal(q), i)));
      this->divShapeValues[q][i] = this->sign[i % this->c.Faces()]
                                   * (1.0 / this->GetTransformationDet(q))
                                   * this->shape.LocalDiv(this->QLocal(q), i);
    }
    for (int i = this->nodalPointsOnFacesTotal; i < this->vShapeValues[q].size(); ++i) {
      this->vShapeValues[q][i] =
          (1.0 / this->GetTransformationDet(q))
          * (this->GetTransformation(q).ApplyJ(this->shape.LocalVector(this->QLocal(q), i)));
      this->divShapeValues[q][i] =
          (1.0 / this->GetTransformationDet(q)) * this->shape.LocalDiv(this->QLocal(q), i);
    }
  }
}

template class RTFaceElementT<>;

#ifdef BUILD_IA

template class RTFaceElementT<IAInterval, SpaceDimension, TimeDimension>;

#endif