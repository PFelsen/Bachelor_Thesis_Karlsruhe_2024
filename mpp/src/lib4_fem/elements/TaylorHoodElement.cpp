#include "TaylorHoodElement.hpp"

template<typename TT, int sDim, int tDim>
TaylorHoodElementT<TT, sDim, tDim>::TaylorHoodElementT(const VectorMatrixBase &base, const Cell &c,
                                                       int face,
                                                       const vector<ShapeValues<TT>> &faceValues) :
    MixedVectorFieldElementT<TT, sDim, tDim>(base, c, face, faceValues) {}

template<typename TT, int sDim, int tDim>
TT TaylorHoodElementT<TT, sDim, tDim>::PressureValue(int q, int i) const {
  return this->Value(q, 1, i);
}

template<typename TT, int sDim, int tDim>
TT TaylorHoodElementT<TT, sDim, tDim>::PressureValue(int q, const ShapeId &iter) const {
  return this->Value(q, 1, iter);
}

template<typename TT, int sDim, int tDim>
TT TaylorHoodElementT<TT, sDim, tDim>::PressureValue(int q, const Vector &u) const {
  return this->Value(q, u, 1);
}

template<typename TT, int sDim, int tDim>
const VectorFieldT<TT, sDim> &TaylorHoodElementT<TT, sDim, tDim>::PressureGradient(int q,
                                                                                   int i) const {
  return this->Derivative(q, 1, i);
}

template<typename TT, int sDim, int tDim>
const VectorFieldT<TT, sDim> &
TaylorHoodElementT<TT, sDim, tDim>::PressureGradient(int q, const ShapeId &iter) const {
  return this->Derivative(q, 1, iter);
}

template<typename TT, int sDim, int tDim>
VectorFieldT<TT, sDim> TaylorHoodElementT<TT, sDim, tDim>::PressureGradient(int q,
                                                                            const Vector &u) const {
  return this->Derivative(q, u, 1);
}

template<typename TT, int sDim, int tDim>
VectorFieldComponentT<TT> TaylorHoodElementT<TT, sDim, tDim>::VelocityFieldComponent(int q, int i,
                                                                                     int k) const {
  return this->VectorComponentValue(q, 0, i, k);
}

template<typename TT, int sDim, int tDim>
VectorFieldComponentT<TT>
TaylorHoodElementT<TT, sDim, tDim>::VelocityFieldComponent(int q, const ShapeId &iter,
                                                           int k) const {
  return this->VectorComponentValue(q, 0, iter, k);
}

template<typename TT, int sDim, int tDim>
VectorFieldT<TT, sDim> TaylorHoodElementT<TT, sDim, tDim>::VelocityField(int q,
                                                                         const Vector &u) const {
  return this->VectorValue(q, u, 0);
}

template<typename TT, int sDim, int tDim>
TensorRowT<TT> TaylorHoodElementT<TT, sDim, tDim>::VelocityFieldRowGradient(int q, int i,
                                                                            int k) const {
  return this->VectorRowGradient(q, 0, i, k);
}

template<typename TT, int sDim, int tDim>
TensorRowT<TT> TaylorHoodElementT<TT, sDim, tDim>::VelocityFieldRowGradient(int q,
                                                                            const ShapeId &iter,
                                                                            int k) const {
  return this->VectorRowGradient(q, 0, iter, k);
}

template<typename TT, int sDim, int tDim>
TensorT<TT, sDim> TaylorHoodElementT<TT, sDim, tDim>::VelocityFieldGradient(int q,
                                                                            const Vector &u) const {
  return this->VectorGradient(q, u, 0);
}

template<typename TT, int sDim, int tDim>
TT TaylorHoodElementT<TT, sDim, tDim>::VelocityFieldDivergence(int q, int i, int k) const {
  return this->Divergence(q, 0, i, k);
}

template<typename TT, int sDim, int tDim>
TT TaylorHoodElementT<TT, sDim, tDim>::VelocityFieldDivergence(int q, const ShapeId &iter,
                                                               int k) const {
  return this->Divergence(q, 0, iter, k);
}

template<typename TT, int sDim, int tDim>
TT TaylorHoodElementT<TT, sDim, tDim>::VelocityFieldDivergence(int q, const Vector &u) const {
  return this->Divergence(q, u, 0);
}

template class TaylorHoodElementT<>;

#ifdef BUILD_IA

template class TaylorHoodElementT<IAInterval, SpaceDimension, TimeDimension>;

#endif

template<typename TT, int sDim, int tDim>
TaylorHoodFaceElementT<TT, sDim, tDim>::TaylorHoodFaceElementT(const VectorMatrixBase &base,
                                                               const Cell &c, int face) :
    IMixedFaceElementT<TT, sDim, tDim>(this->mixedDoF.NumberOfDoFs()),
    TaylorHoodElementT<TT, sDim, tDim>(base, c, face, this->values()) {
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

template class TaylorHoodFaceElementT<>;

#ifdef BUILD_IA

template class TaylorHoodFaceElementT<IAInterval, SpaceDimension, TimeDimension>;

#endif