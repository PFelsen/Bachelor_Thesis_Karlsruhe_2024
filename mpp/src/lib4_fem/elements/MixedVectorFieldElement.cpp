#include "MixedVectorFieldElement.hpp"

template<typename TT, int sDim, int tDim>
MixedVectorFieldElementT<TT, sDim, tDim>::MixedVectorFieldElementT(
    const VectorMatrixBase &base, const Cell &c, int face,
    const vector<ShapeValues<TT>> &faceValues) :
    MixedScalarElementT<TT, sDim, tDim>(base, c, face, faceValues), dim(c.dim()) {}

template<typename TT, int sDim, int tDim>
MixedVectorFieldElementT<TT, sDim, tDim>::MixedVectorFieldElementT(const VectorMatrixBase &base,
                                                                   const Cell &c) :
    MixedScalarElementT<TT, sDim, tDim>(base, c), dim(c.dim()) {}

template<typename TT, int sDim, int tDim>
MixedVectorFieldElementT<TT, sDim, tDim>::MixedVectorFieldElementT(
    const VectorMatrixBase &base, const BasicElementT<TT, sDim, tDim> &baseElement) :
    MixedScalarElementT<TT, sDim, tDim>(base, baseElement), dim(this->c.dim()) {}

template<typename TT, int sDim, int tDim>
VectorFieldComponentT<TT>
MixedVectorFieldElementT<TT, sDim, tDim>::VectorComponentValue(int q, int n, int i, int k) const {
  return VectorFieldComponentT<TT>(k, this->Value(q, n, i));
}

template<typename TT, int sDim, int tDim>
VectorFieldComponentT<TT>
MixedVectorFieldElementT<TT, sDim, tDim>::VectorComponentValue(int q, int n, const ShapeId &iter,
                                                               int k) const {
  return VectorFieldComponentT<TT>(k, this->Value(q, n, iter));
}

template<typename TT, int sDim, int tDim>
VectorFieldComponentT<TT>
MixedVectorFieldElementT<TT, sDim, tDim>::VectorComponentValue(const PointT<TT, sDim, tDim> &z,
                                                               int n, int i, int k) const {
  return VectorFieldComponentT<TT>(k, this->Value(z, n, i));
}

template<typename TT, int sDim, int tDim>
VectorFieldComponentT<TT> MixedVectorFieldElementT<TT, sDim, tDim>::VectorComponentValue(
    const PointT<TT, sDim, tDim> &z, int n, const ShapeId &iter, int k) const {
  return VectorFieldComponentT<TT>(k, this->Value(z, n, iter));
}

template<typename TT, int sDim, int tDim>
VectorFieldT<TT, sDim> MixedVectorFieldElementT<TT, sDim, tDim>::VectorValue(int q, const Vector &u,
                                                                             int n) const {
  VectorFieldT<TT, sDim> V;
  for (int i = 0; i < this->size(n); ++i)
    for (int k = 0; k < dim; ++k)
      V[k] += this->Value(q, n, i)
              * u(this->r(this->nodalPointData[n][i].index), this->nodalPointData[n][i].shift + k);
  return V;
}

template<typename TT, int sDim, int tDim>
VectorFieldT<TT, sDim>
MixedVectorFieldElementT<TT, sDim, tDim>::VectorValue(const PointT<TT, sDim, tDim> &z,
                                                      const Vector &u, int n) const {
  VectorFieldT<TT, sDim> V;
  for (int i = 0; i < this->size(n); ++i) {
    TT v = this->Value(z, n, i);
    for (int k = 0; k < dim; ++k)
      V[k] +=
          v * u(this->r(this->nodalPointData[n][i].index), this->nodalPointData[n][i].shift + k);
  }
  return V;
}

template<typename TT, int sDim, int tDim>
TensorRowT<TT> MixedVectorFieldElementT<TT, sDim, tDim>::VectorRowGradient(int q, int n, int i,
                                                                           int k) const {
  return TensorRowT<TT>(this->Derivative(q, n, i), k);
}

template<typename TT, int sDim, int tDim>
TensorRowT<TT> MixedVectorFieldElementT<TT, sDim, tDim>::VectorRowGradient(int q, int n,
                                                                           const ShapeId &iter,
                                                                           int k) const {
  return TensorRowT<TT>(this->Derivative(q, n, iter), k);
}

template<typename TT, int sDim, int tDim>
TensorRowT<TT>
MixedVectorFieldElementT<TT, sDim, tDim>::VectorRowGradient(const PointT<TT, sDim, tDim> &z, int n,
                                                            int i, int k) const {
  return TensorRowT<TT>(this->Derivative(z, n, i), k);
}

template<typename TT, int sDim, int tDim>
TensorRowT<TT>
MixedVectorFieldElementT<TT, sDim, tDim>::VectorRowGradient(const PointT<TT, sDim, tDim> &z, int n,
                                                            const ShapeId &iter, int k) const {
  return TensorRowT<TT>(this->Derivative(z, n, iter), k);
}

template<typename TT, int sDim, int tDim>
TensorT<TT, sDim> MixedVectorFieldElementT<TT, sDim, tDim>::VectorGradient(int q, const Vector &u,
                                                                           int n) const {
  TensorT<TT, sDim> T;
  for (int i = 0; i < this->size(n); ++i)
    for (int k = 0; k < dim; ++k)
      T[k] += u(this->r(this->nodalPointData[n][i].index), this->nodalPointData[n][i].shift + k)
              * this->Derivative(q, n, i);
  return T;
}

template<typename TT, int sDim, int tDim>
TensorT<TT, sDim>
MixedVectorFieldElementT<TT, sDim, tDim>::VectorGradient(const PointT<TT, sDim, tDim> &z,
                                                         const Vector &u, int n) const {
  TensorT<TT, sDim> T;
  for (int i = 0; i < this->size(n); ++i) {
    VectorFieldT<TT, sDim> D = this->Derivative(z, n, i);
    for (int k = 0; k < dim; ++k)
      T[k] +=
          u(this->r(this->nodalPointData[n][i].index), this->nodalPointData[n][i].shift + k) * D;
  }
  return T;
}

template<typename TT, int sDim, int tDim>
TT MixedVectorFieldElementT<TT, sDim, tDim>::Divergence(int q, int n, int i, int k) const {
  return this->Derivative(q, n, i)[k];
}

template<typename TT, int sDim, int tDim>
TT MixedVectorFieldElementT<TT, sDim, tDim>::Divergence(int q, int n, const ShapeId &iter,
                                                        int k) const {
  return this->Derivative(q, n, iter)[k];
}

template<typename TT, int sDim, int tDim>
TT MixedVectorFieldElementT<TT, sDim, tDim>::Divergence(const PointT<TT, sDim, tDim> &z, int n,
                                                        int i, int k) const {
  return this->Derivative(z, n, i)[k];
}

template<typename TT, int sDim, int tDim>
TT MixedVectorFieldElementT<TT, sDim, tDim>::Divergence(const PointT<TT, sDim, tDim> &z, int n,
                                                        const ShapeId &iter, int k) const {
  return this->Derivative(z, n, iter)[k];
}

template<typename TT, int sDim, int tDim>
TT MixedVectorFieldElementT<TT, sDim, tDim>::Divergence(int q, const Vector &u, int n) const {
  TT s = 0.0;
  for (int i = 0; i < this->size(n); ++i) {
    VectorFieldT<TT, sDim> G = this->Derivative(q, n, i);
    for (int k = 0; k < dim; ++k)
      s +=
          u(this->r(this->nodalPointData[n][i].index), this->nodalPointData[n][i].shift + k) * G[k];
  }
  return s;
}

template<typename TT, int sDim, int tDim>
TT MixedVectorFieldElementT<TT, sDim, tDim>::Divergence(const PointT<TT, sDim, tDim> &z,
                                                        const Vector &u, int n) const {
  TT s = 0.0;
  for (int i = 0; i < this->size(n); ++i) {
    VectorFieldT<TT, sDim> G = this->Derivative(z, n, i);
    for (int k = 0; k < dim; ++k)
      s +=
          u(this->r(this->nodalPointData[n][i].index), this->nodalPointData[n][i].shift + k) * G[k];
  }
  return s;
}


template class MixedVectorFieldElementT<>;

#ifdef BUILD_IA

template class MixedVectorFieldElementT<IAInterval, SpaceDimension, TimeDimension>;

#endif

template<typename TT, int sDim, int tDim>
MixedVectorFieldFaceElementT<TT, sDim, tDim>::MixedVectorFieldFaceElementT(
    const VectorMatrixBase &base, const Cell &c, int face) :
    IMixedFaceElementT<TT, sDim, tDim>(this->mixedDoF.NumberOfDoFs()),
    MixedVectorFieldElementT<TT, sDim, tDim>(base, c, face, this->values()) {
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
  //  }
}

template class MixedVectorFieldFaceElementT<>;

#ifdef BUILD_IA

template class MixedVectorFieldFaceElementT<IAInterval, SpaceDimension, TimeDimension>;

#endif