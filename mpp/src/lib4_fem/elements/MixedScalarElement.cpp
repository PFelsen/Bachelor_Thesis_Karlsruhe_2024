#include "MixedScalarElement.hpp"

template<typename TT, int sDim, int tDim>
void MixedScalarElementT<TT, sDim, tDim>::init() {
  for (int n = 0; n < shapeValues.size(); ++n) {
    shapeValues[n] = &(this->shapes[n]->values());
    shapeGradients[n].resize(this->nQ(), this->size(n));
    for (int q = 0; q < this->nQ(); ++q)
      for (int i = 0; i < this->size(n); ++i)
        shapeGradients[n][q][i] = this->GetTransformation(q) * this->shapes[n]->LocalGradient(q, i);
  }
}

template<typename TT, int sDim, int tDim>
MixedScalarElementT<TT, sDim, tDim>::MixedScalarElementT(const VectorMatrixBase &base,
                                                         const Cell &c) :
    MixedElementT<TT, sDim, tDim>(base, c), shapeGradients(this->mixedDoF.NumberOfDoFs()),
    shapeValues(this->mixedDoF.NumberOfDoFs()) {
  init();
}

template<typename TT, int sDim, int tDim>
MixedScalarElementT<TT, sDim, tDim>::MixedScalarElementT(
    const VectorMatrixBase &base, const BasicElementT<TT, sDim, tDim> &baseElement) :
    MixedElementT<TT, sDim, tDim>(base, baseElement), shapeGradients(this->mixedDoF.NumberOfDoFs()),
    shapeValues(this->mixedDoF.NumberOfDoFs()) {
  init();
}

template<typename TT, int sDim, int tDim>
MixedScalarElementT<TT, sDim, tDim>::MixedScalarElementT(
    const VectorMatrixBase &base, const Cell &c, int face,
    const vector<ShapeValues<TT>> &faceValues) :
    MixedElementT<TT, sDim, tDim>(base, c, face), shapeValues(faceValues.size()),
    shapeGradients(this->mixedDoF.NumberOfDoFs()) {
  for (int n = 0; n < faceValues.size(); ++n) {
    shapeValues[n] = &(faceValues[n]);
  }
}

template<typename TT, int sDim, int tDim>
TT MixedScalarElementT<TT, sDim, tDim>::Value(int q, int n, int i) const {
  return (*shapeValues[n])[q][i];
}

template<typename TT, int sDim, int tDim>
TT MixedScalarElementT<TT, sDim, tDim>::Value(int q, int n, const ShapeId &iter) const {
  return (*shapeValues[n])[q][iter.id];
}

template<typename TT, int sDim, int tDim>
TT MixedScalarElementT<TT, sDim, tDim>::Value(const PointT<TT, sDim, tDim> &z, int n, int i) const {
  return (*this->shapes[n])(z, i);
}

template<typename TT, int sDim, int tDim>
TT MixedScalarElementT<TT, sDim, tDim>::Value(const PointT<TT, sDim, tDim> &z, int n,
                                              const ShapeId &iter) const {
  return (*this->shapes[n])(z, iter.id);
}

template<typename TT, int sDim, int tDim>
TT MixedScalarElementT<TT, sDim, tDim>::Value(int q, const Vector &u, int n, int k) const {
  TT U = 0.0;
  for (int i = 0; i < this->size(n); ++i)
    U += u(this->r(this->nodalPointData[n][i].index), this->nodalPointData[n][i].shift + k)
         * (*shapeValues[n])[q][i];
  return U;
}

template<typename TT, int sDim, int tDim>
TT MixedScalarElementT<TT, sDim, tDim>::Value(const PointT<TT, sDim, tDim> &z, const Vector &u,
                                              int n, int k) const {
  TT U = 0.0;
  for (int i = 0; i < this->size(n); ++i)
    U += u(this->r(this->nodalPointData[n][i].index), this->nodalPointData[n][i].shift + k)
         * (*this->shapes[n])(z, i);
  return U;
}

template<typename TT, int sDim, int tDim>
const VectorFieldT<TT, sDim> &MixedScalarElementT<TT, sDim, tDim>::Derivative(int q, int n,
                                                                              int i) const {
  return shapeGradients[n][q][i];
}

template<typename TT, int sDim, int tDim>
const VectorFieldT<TT, sDim> &
MixedScalarElementT<TT, sDim, tDim>::Derivative(int q, int n, const ShapeId &iter) const {
  return shapeGradients[n][q][iter.id];
}

template<typename TT, int sDim, int tDim>
VectorFieldT<TT, sDim>
MixedScalarElementT<TT, sDim, tDim>::Derivative(const PointT<TT, sDim, tDim> &z, int n,
                                                int i) const {
  return this->GetTransformation(z) * this->shapes[n]->LocalGradient(z, i);
}

template<typename TT, int sDim, int tDim>
VectorFieldT<TT, sDim>
MixedScalarElementT<TT, sDim, tDim>::Derivative(const PointT<TT, sDim, tDim> &z, int n,
                                                const ShapeId &iter) const {
  return this->GetTransformation(z) * this->shapes[n]->LocalGradient(z, iter.id);
}

template<typename TT, int sDim, int tDim>
VectorFieldT<TT, sDim> MixedScalarElementT<TT, sDim, tDim>::Derivative(int q, const Vector &u,
                                                                       int n, int k) const {
  VectorFieldT<TT, sDim> Du;
  for (int i = 0; i < this->size(n); ++i)
    Du += u(this->r(this->nodalPointData[n][i].index), this->nodalPointData[n][i].shift + k)
          * shapeGradients[n][q][i];
  return Du;
}

template<typename TT, int sDim, int tDim>
VectorFieldT<TT, sDim>
MixedScalarElementT<TT, sDim, tDim>::Derivative(const PointT<TT, sDim, tDim> &z, const Vector &u,
                                                int n, int k) const {
  VectorFieldT<TT, sDim> Du;
  const TransformationT<TT, sDim> &T = this->GetTransformation(z);
  for (int i = 0; i < this->size(n); ++i)
    Du += u(this->r(this->nodalPointData[n][i].index), this->nodalPointData[n][i].shift + k)
          * (T * this->shapes[n]->LocalGradient(z, i));
  return Du;
}


template class MixedScalarElementT<>;

#ifdef BUILD_IA

template class MixedScalarElementT<IAInterval, SpaceDimension, TimeDimension>;

#endif

template<typename TT, int sDim, int tDim>
MixedScalarFaceElementT<TT, sDim, tDim>::MixedScalarFaceElementT(const VectorMatrixBase &base,
                                                                 const Cell &c, int face) :
    IMixedFaceElementT<TT, sDim, tDim>(this->mixedDoF.NumberOfDoFs()),
    MixedScalarElementT<TT, sDim, tDim>(base, c, face, this->values()) {
  for (int n = 0; n < this->shapeValues.size(); ++n) {
    this->faceValues[n].resize(this->nQ(), this->size(n));
    this->shapeGradients[n].resize(this->nQ(), this->size(n));
    for (int q = 0; q < this->nQ(); ++q) {
      for (int i = 0; i < this->size(n); ++i) {
        this->SetFaceValue(q, n, i, (*this->shapes[n])(this->QLocal(q), i));
        this->shapeGradients[n][q][i] =
            this->GetTransformation(q) * this->shapes[n]->LocalGradient(this->QLocal(q), i);
      }
    }
  }
}


template class MixedScalarFaceElementT<>;

#ifdef BUILD_IA

template class MixedScalarFaceElementT<IAInterval, SpaceDimension, TimeDimension>;

#endif