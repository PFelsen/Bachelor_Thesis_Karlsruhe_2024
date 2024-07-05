#include "RTLagrangeElement.hpp"
#include "RTLagrangeDiscretization.hpp"

template<typename TT, int sDim, int tDim>
void RTLagrangeElementT<TT, sDim, tDim>::init(const VectorMatrixBase &base) {
  velocitySize = this->mixedDoF.NumberOfStoragePoints(VelocityIndex(), this->c);
  lagrangeDegree = dynamic_cast<const RTLagrangeDiscretization &>(base.GetDisc()).LagrangeDegree();
  nodalPointsOnFacesTotalRT = 0;
  for (int face = 0; face < this->c.Faces(); ++face) {
    TransformationT<TT, sDim> T = this->GetTransformation((*this)[face]());
    // TODO: interval arithmetic normal
    normal[face] = T * PointT<TT, sDim, tDim>(this->c.LocalFaceNormal(face));
    normal[face] /= norm(normal[face]);
    if (mid(normal[face]) > Origin) sign[face] = 1;
    else sign[face] = -1;
    nodalPointsOnFacesTotalRT +=
        this->mixedDoF.NumberOfStoragePointsOnFace(VelocityIndex(), this->c, face);
  }
}

template<typename TT, int sDim, int tDim>
void RTLagrangeElementT<TT, sDim, tDim>::initValues(const VectorMatrixBase &base) {
  vShapeValues.resize(this->nQ(), 2 * VelocitySize() - nodalPointsOnFacesTotalRT);
  divShapeValues.resize(this->nQ(), 2 * VelocitySize() - nodalPointsOnFacesTotalRT);
  for (int q = 0; q < this->nQ(); ++q) {
    for (int i = 0; i < nodalPointsOnFacesTotalRT; ++i) {
      vShapeValues[q][i] =
          sign[i % this->c.Faces()] * (1.0 / this->GetTransformationDet(q))
          * (this->GetTransformation(q).ApplyJ(this->shapes[0]->LocalVector(q, i)));
      divShapeValues[q][i] = sign[i % this->c.Faces()] * (1.0 / this->GetTransformationDet(q))
                             * this->shapes[0]->LocalDiv(q, i);
    }
    for (int i = nodalPointsOnFacesTotalRT; i < vShapeValues[q].size(); ++i) {
      vShapeValues[q][i] =
          (1.0 / this->GetTransformationDet(q))
          * (this->GetTransformation(q).ApplyJ(this->shapes[0]->LocalVector(q, i)));
      divShapeValues[q][i] =
          (1.0 / this->GetTransformationDet(q)) * this->shapes[0]->LocalDiv(q, i);
    }
  }

  if (lagrangeDegree > 0) {
    pShapeGradients.resize(this->nQ(), this->size(1));
    for (int q = 0; q < this->nQ(); ++q)
      for (int i = 0; i < this->size(1); ++i)
        pShapeGradients[q][i] = this->GetTransformation(q) * this->shapes[1]->LocalGradient(q, i);
  }
}

template<typename TT, int sDim, int tDim>
RTLagrangeElementT<TT, sDim, tDim>::RTLagrangeElementT(const VectorMatrixBase &base, const Cell &c,
                                                       const ShapeValues<TT> &pShapeValuesFace,
                                                       int face) :
    MixedElementT<TT, sDim, tDim>(base, c, face), dim(c.dim()), normal(c.Faces()), sign(c.Faces()),
    vectorfield_face(c.Faces()), pShapeValues(pShapeValuesFace),
    vShapeValues(this->nQ(), this->size(0)), divShapeValues(this->nQ(), this->size(0)) {
  init(base);
}

template<typename TT, int sDim, int tDim>
RTLagrangeElementT<TT, sDim, tDim>::RTLagrangeElementT(const VectorMatrixBase &base,
                                                       const Cell &c) :
    MixedElementT<TT, sDim, tDim>(base, c), dim(c.dim()), normal(c.Faces()), sign(c.Faces()),
    vectorfield_face(c.Faces()), pShapeValues(this->shapes[1]->values()),
    vShapeValues(this->nQ(), this->size(0)), divShapeValues(this->nQ(), this->size(0)) {
  init(base);
  initValues(base);
}

template<typename TT, int sDim, int tDim>
TT RTLagrangeElementT<TT, sDim, tDim>::PressureValue(int q, const Vector &u, int m) const {
  TT p = 0.0;
  for (int i = 0; i < this->PressureSize(); ++i)
    p += u(this->r(this->nodalPointData[PressureIndex()][i].index),
           this->nodalPointData[PressureIndex()][i].shift + m)
         * pShapeValues[q][i];
  return p;
}

template<typename TT, int sDim, int tDim>
TT RTLagrangeElementT<TT, sDim, tDim>::PressureValue(const PointT<TT, sDim, tDim> &z,
                                                     const Vector &u, int m) const {
  TT p = 0.0;
  for (int i = 0; i < this->PressureSize(); ++i)
    p += u(this->r(this->nodalPointData[PressureIndex()][i].index),
           this->nodalPointData[PressureIndex()][i].shift + m)
         * (*this->shapes[PressureIndex()])(z, i);
  return p;
}

template<typename TT, int sDim, int tDim>
const VectorFieldT<TT, sDim> &RTLagrangeElementT<TT, sDim, tDim>::PressureGradient(int q,
                                                                                   int i) const {
  if (lagrangeDegree == 0) {
    THROW("PressureGradient does not exist for degree 0")
  } else {
    return pShapeGradients[q][i];
  }
}

template<typename TT, int sDim, int tDim>
const VectorFieldT<TT, sDim> &
RTLagrangeElementT<TT, sDim, tDim>::PressureGradient(int q, const ShapeId &iter) const {
  if (lagrangeDegree == 0) {
    THROW("PressureGradient does not exist for degree 0")
  } else {
    return pShapeGradients[q][iter.id];
  }
}

template<typename TT, int sDim, int tDim>
VectorFieldT<TT, sDim> RTLagrangeElementT<TT, sDim, tDim>::PressureGradient(int q, const Vector &u,
                                                                            int m) const {
  if (lagrangeDegree == 0) {
    THROW("PressureGradient does not exist for degree 0")
  } else {
    VectorFieldT<TT, sDim> D;
    for (int i = 0; i < this->PressureSize(); ++i)
      D += u(this->r(this->nodalPointData[PressureIndex()][i].index),
             this->nodalPointData[PressureIndex()][i].shift + m)
           * pShapeGradients[q][i];
    return D;
  }
}

template<typename TT, int sDim, int tDim>
VectorFieldT<TT, sDim>
RTLagrangeElementT<TT, sDim, tDim>::PressureGradient(const PointT<TT, sDim, tDim> &z,
                                                     const Vector &u, int m) const {
  if (lagrangeDegree == 0) {
    THROW("PressureGradient does not exist for degree 0")
  } else {
    TransformationT<TT, sDim> trafo = this->GetTransformation(z);
    VectorFieldT<TT, sDim> D;
    for (int i = 0; i < this->PressureSize(); ++i)
      D += u(this->r(this->nodalPointData[PressureIndex()][i].index),
             this->nodalPointData[PressureIndex()][i].shift + m)
           * (trafo * this->shapes[1]->LocalGradient(z, i));
    return D;
  }
}

template<typename TT, int sDim, int tDim>
VectorFieldT<TT, sDim> RTLagrangeElementT<TT, sDim, tDim>::VelocityField(int q,
                                                                         const Vector &u) const {
  VectorFieldT<TT, sDim> V;
  for (int i = 0; i < this->VelocitySize(); ++i)
    for (int k = 0; k < VelocityMaxk(i); ++k)
      V += u(this->r(this->nodalPointData[VelocityIndex()][i].index),
             this->nodalPointData[VelocityIndex()][i].shift + k)
           * vShapeValues[q][IndexingVelocityShape(i, k)];
  return V;
}

template<typename TT, int sDim, int tDim>
VectorFieldT<TT, sDim> RTLagrangeElementT<TT, sDim, tDim>::VelocityField(int q, const Vector &u,
                                                                         int m) const {
  VectorFieldT<TT, sDim> V;
  for (int i = 0; i < this->VelocitySize(); ++i)
    for (int k = 0; k < VelocityMaxk(i); ++k)
      V += u(this->r(this->nodalPointData[VelocityIndex()][i].index),
             this->nodalPointData[VelocityIndex()][i].shift + IndexingVelocity_k(i, k, m))
           * vShapeValues[q][IndexingVelocityShape(i, k)];
  return V;
}

template<typename TT, int sDim, int tDim>
VectorFieldT<TT, sDim>
RTLagrangeElementT<TT, sDim, tDim>::VelocityField(const PointT<TT, sDim, tDim> &z,
                                                  const Vector &u) const {
  std::vector<VectorFieldT<TT, sDim>> values(2 * VelocitySize() - nodalPointsOnFacesTotalRT);
  TransformationT<TT, sDim> trafo = this->GetTransformation(z);

  for (int i = 0; i < nodalPointsOnFacesTotalRT; ++i) {
    values[i] = sign[i % this->c.Faces()] * (1.0 / trafo.Det())
                * (trafo.ApplyJ(this->shapes[0]->LocalVector(z, i)));
  }
  for (int i = nodalPointsOnFacesTotalRT; i < values.size(); ++i) {
    values[i] = (1.0 / trafo.Det()) * (trafo.ApplyJ(this->shapes[0]->LocalVector(z, i)));
  }
  VectorFieldT<TT, sDim> V;
  for (int i = 0; i < this->VelocitySize(); ++i)
    for (int k = 0; k < VelocityMaxk(i); ++k)
      V += u(this->r(this->nodalPointData[VelocityIndex()][i].index),
             this->nodalPointData[VelocityIndex()][i].shift + k)
           * values[IndexingVelocityShape(i, k)];
  return V;
}

template<typename TT, int sDim, int tDim>
VectorFieldT<TT, sDim>
RTLagrangeElementT<TT, sDim, tDim>::VelocityField(const PointT<TT, sDim, tDim> &z, const Vector &u,
                                                  int m) const {
  std::vector<VectorFieldT<TT, sDim>> values(2 * VelocitySize() - nodalPointsOnFacesTotalRT);
  TransformationT<TT, sDim> trafo = this->GetTransformation(z);

  for (int i = 0; i < nodalPointsOnFacesTotalRT; ++i) {
    values[i] = sign[i % this->c.Faces()] * (1.0 / trafo.Det())
                * (trafo.ApplyJ(this->shapes[0]->LocalVector(z, i)));
  }
  for (int i = nodalPointsOnFacesTotalRT; i < values.size(); ++i) {
    values[i] = (1.0 / trafo.Det()) * (trafo.ApplyJ(this->shapes[0]->LocalVector(z, i)));
  }
  VectorFieldT<TT, sDim> V;
  for (int i = 0; i < this->VelocitySize(); ++i)
    for (int k = 0; k < VelocityMaxk(i); ++k)
      V += u(this->r(this->nodalPointData[VelocityIndex()][i].index),
             this->nodalPointData[VelocityIndex()][i].shift + IndexingVelocity_k(i, k, m))
           * values[IndexingVelocityShape(i, k)];
  return V;
}

template<typename TT, int sDim, int tDim>
TT RTLagrangeElementT<TT, sDim, tDim>::VelocityFieldDivergence(int q, const Vector &u) const {
  TT d = 0.0;
  for (int i = 0; i < this->VelocitySize(); ++i)
    for (int k = 0; k < VelocityMaxk(i); ++k)
      d += u(this->r(this->nodalPointData[VelocityIndex()][i].index),
             this->nodalPointData[VelocityIndex()][i].shift + k)
           * divShapeValues[q][IndexingVelocityShape(i, k)];
  return d;
}

template<typename TT, int sDim, int tDim>
TT RTLagrangeElementT<TT, sDim, tDim>::VelocityFieldDivergence(int q, const Vector &u,
                                                               int m) const {
  TT d = 0.0;
  for (int i = 0; i < this->VelocitySize(); ++i)
    for (int k = 0; k < VelocityMaxk(i); ++k)
      d += u(this->r(this->nodalPointData[VelocityIndex()][i].index),
             this->nodalPointData[VelocityIndex()][i].shift + IndexingVelocity_k(i, k, m))
           * divShapeValues[q][IndexingVelocityShape(i, k)];
  return d;
}

template<typename TT, int sDim, int tDim>
TT RTLagrangeElementT<TT, sDim, tDim>::VelocityFieldDivergence(const PointT<TT, sDim, tDim> &z,
                                                               const Vector &u) const {
  std::vector<TT> values(2 * VelocitySize() - nodalPointsOnFacesTotalRT);
  TransformationT<TT, sDim> trafo = this->GetTransformation(z);

  for (int i = 0; i < nodalPointsOnFacesTotalRT; ++i) {
    values[i] = sign[i % this->c.Faces()] * (1.0 / trafo.Det()) * this->shapes[0]->LocalDiv(z, i);
  }
  for (int i = nodalPointsOnFacesTotalRT; i < values.size(); ++i) {
    values[i] = (1.0 / trafo.Det()) * this->shapes[0]->LocalDiv(z, i);
  }

  TT d = 0.0;
  for (int i = 0; i < this->VelocitySize(); ++i)
    for (int k = 0; k < VelocityMaxk(i); ++k)
      d += u(this->r(this->nodalPointData[VelocityIndex()][i].index),
             this->nodalPointData[VelocityIndex()][i].shift + k)
           * values[IndexingVelocityShape(i, k)];
  return d;
}

template<typename TT, int sDim, int tDim>
TT RTLagrangeElementT<TT, sDim, tDim>::VelocityFieldDivergence(const PointT<TT, sDim, tDim> &z,
                                                               const Vector &u, int m) const {
  std::vector<TT> values(2 * VelocitySize() - nodalPointsOnFacesTotalRT);
  TransformationT<TT, sDim> trafo = this->GetTransformation(z);

  for (int i = 0; i < nodalPointsOnFacesTotalRT; ++i) {
    values[i] = sign[i % this->c.Faces()] * (1.0 / trafo.Det()) * this->shapes[0]->LocalDiv(z, i);
  }
  for (int i = nodalPointsOnFacesTotalRT; i < values.size(); ++i) {
    values[i] = (1.0 / trafo.Det()) * this->shapes[0]->LocalDiv(z, i);
  }

  TT d = 0.0;
  for (int i = 0; i < this->VelocitySize(); ++i)
    for (int k = 0; k < VelocityMaxk(i); ++k)
      d += u(this->r(this->nodalPointData[VelocityIndex()][i].index),
             this->nodalPointData[VelocityIndex()][i].shift + IndexingVelocity_k(i, k, m))
           * values[IndexingVelocityShape(i, k)];
  return d;
}

template<typename TT, int sDim, int tDim>
VectorFieldT<TT, sDim> RTLagrangeElementT<TT, sDim, tDim>::CellFlux(const Vector &u) const {
  VectorFieldT<TT, sDim> F = zero;
  TT area = 0;
  for (int q = 0; q < this->nQ(); ++q) {
    TT w = this->QWeight(q);
    area += w;
    F += w * this->VelocityField(q, u);
  }
  F *= (1 / area);
  return F;
}

template class RTLagrangeElementT<>;

#ifdef BUILD_IA

template class RTLagrangeElementT<IAInterval, SpaceDimension, TimeDimension>;

#endif


template<typename TT, int sDim, int tDim>
RTLagrangeFaceElementT<TT, sDim, tDim>::RTLagrangeFaceElementT(const VectorMatrixBase &base,
                                                               const Cell &c, int face) :
    IMixedFaceElementT<TT, sDim, tDim>(1),
    RTLagrangeElementT<TT, sDim, tDim>(base, c, this->values()[0], face) {
  this->vShapeValues.resize(this->nQ(), 2 * this->VelocitySize() - this->nodalPointsOnFacesTotalRT);
  this->divShapeValues.resize(this->nQ(),
                              2 * this->VelocitySize() - this->nodalPointsOnFacesTotalRT);
  for (int q = 0; q < this->nQ(); ++q) {
    for (int i = 0; i < this->nodalPointsOnFacesTotalRT; ++i) {
      this->vShapeValues[q][i] =
          this->sign[i % this->c.Faces()] * (1.0 / this->GetTransformationDet(q))
          * (this->GetTransformation(q).ApplyJ(this->shapes[0]->LocalVector(this->QLocal(q), i)));
      this->divShapeValues[q][i] = this->sign[i % this->c.Faces()]
                                   * (1.0 / this->GetTransformationDet(q))
                                   * this->shapes[0]->LocalDiv(this->QLocal(q), i);
    }
    for (int i = this->nodalPointsOnFacesTotalRT; i < this->vShapeValues[q].size(); ++i) {
      this->vShapeValues[q][i] =
          (1.0 / this->GetTransformationDet(q))
          * (this->GetTransformation(q).ApplyJ(this->shapes[0]->LocalVector(this->QLocal(q), i)));
      this->divShapeValues[q][i] =
          (1.0 / this->GetTransformationDet(q)) * this->shapes[0]->LocalDiv(this->QLocal(q), i);
    }
  }
  for (int q = 0; q < this->nQ(); ++q) {
    for (int i = 0; i < this->size(1); ++i) {
      this->faceValues[0].resize(this->nQ(), this->size(1));
      this->SetFaceValue(q, 0, i, (*this->shapes[1])(this->QLocal(q), i));
    }
  }
  if (this->lagrangeDegree > 0) {
    this->pShapeGradients.resize(this->nQ(), this->size(1));
    for (int q = 0; q < this->nQ(); ++q)
      for (int i = 0; i < this->size(1); ++i)
        this->pShapeGradients[q][i] =
            this->GetTransformation(q) * this->shapes[1]->LocalGradient(this->QLocal(q), i);
  }
}

template class RTLagrangeFaceElementT<>;

#ifdef BUILD_IA

template class RTLagrangeFaceElementT<IAInterval, SpaceDimension, TimeDimension>;

#endif