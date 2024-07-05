#include "VectorFieldElement.hpp"

template<typename TT, int sDim, int tDim>
VectorFieldElementT<TT, sDim, tDim>::VectorFieldElementT(const VectorMatrixBase &base,
                                                         const Cell &c, int face,
                                                         const ShapeValues<TT> &faceValues) :
    ScalarElementT<TT, sDim, tDim>(base, c, face, faceValues), dim(c.dim()) {}

template<typename TT, int sDim, int tDim>
VectorFieldElementT<TT, sDim, tDim>::VectorFieldElementT(const VectorMatrixBase &base,
                                                         const Cell &c) :
    ScalarElementT<TT, sDim, tDim>(base, c), dim(c.dim()) {}

template<typename TT, int sDim, int tDim>
VectorFieldElementT<TT, sDim, tDim>::VectorFieldElementT(
    const VectorMatrixBase &base, const BasicElementT<TT, sDim, tDim> &baseElement) :
    ScalarElementT<TT, sDim, tDim>(base, baseElement), dim(this->c.dim()) {}

template<typename TT, int sDim, int tDim>
VectorFieldComponentT<TT> VectorFieldElementT<TT, sDim, tDim>::VectorComponentValue(int q, int i,
                                                                                    int k) const {
  return VectorFieldComponentT<TT>(k, this->Value(q, i));
}

template<typename TT, int sDim, int tDim>
VectorFieldT<TT, sDim, tDim>
VectorFieldElementT<TT, sDim, tDim>::VectorValue(int q, const Vector &u) const {
  VectorFieldT<TT, sDim> V;
  for (int i = 0; i < this->size(); ++i) {
    TT s = this->Value(q, i);
    for (int k = 0; k < dim; ++k)
      V[k] += s * u(this->r(i), k);
  }
  return V;
}

template<typename TT, int sDim, int tDim>
VectorFieldT<TT, sDim, tDim>
VectorFieldElementT<TT, sDim, tDim>::VectorValue(const PointT<TT, sDim, tDim> &z,
                                                 const Vector &u) const {
  VectorFieldT<TT, sDim> V;
  for (int k = 0; k < dim; ++k)
    V[k] = this->Value(z, u, k);
  return V;
}

template<typename TT, int sDim, int tDim>
VectorFieldT<TT, sDim, tDim>
VectorFieldElementT<TT, sDim, tDim>::VectorValue(const PointT<TT, sDim, tDim> &z, int i) const {
  VectorFieldT<TT, sDim> V;
  return VectorFieldT<TT, sDim, tDim>(this->Value(z, i));
}

template<typename TT, int sDim, int tDim>
TensorRowT<TT, sDim> VectorFieldElementT<TT, sDim, tDim>::VectorRowGradient(int q, int i,
                                                                            int k) const {
  return TensorRowT<TT, sDim>(k, this->Derivative(q, i));
}

template<typename TT, int sDim, int tDim>
TensorT<TT, sDim> VectorFieldElementT<TT, sDim, tDim>::VectorGradient(int q,
                                                                      const Vector &u) const {
  TensorT<TT, sDim> T;
  for (int i = 0; i < this->size(); ++i)
    for (int k = 0; k < dim; ++k) {
      VectorFieldT<TT, sDim> G = u(this->r(i), k) * this->Derivative(q, i);
      for (int l = 0; l < dim; ++l)
        T[k][l] += G[l];
    }
  return T;
}

template<typename TT, int sDim, int tDim>
TensorT<TT, sDim>
VectorFieldElementT<TT, sDim, tDim>::VectorGradient(const PointT<TT, sDim, tDim> &z,
                                                    const Vector &u) const {
  TensorT<TT, sDim> T;
  for (int k = 0; k < dim; ++k) {
    VectorFieldT<TT, sDim> G = this->Derivative(z, u, k);
    for (int l = 0; l < dim; ++l)
      T[k][l] += G[l];
  }
  return T;
}

template<typename TT, int sDim, int tDim>
TT VectorFieldElementT<TT, sDim, tDim>::Divergence(int q, int i, int k) const {
  VectorFieldT<TT, sDim> G = this->Derivative(q, i);
  return G[k];
}

template<typename TT, int sDim, int tDim>
TT VectorFieldElementT<TT, sDim, tDim>::Divergence(int q, const Vector &u) const {
  TT s = TT(0.0);
  for (int i = 0; i < this->size(); ++i) {
    VectorFieldT<TT, sDim> G = this->Derivative(q, i);
    for (int k = 0; k < dim; ++k)
      s += u(this->r(i), k) * G[k];
  }
  return s;
}


template class VectorFieldElementT<>;

#ifdef BUILD_IA

template class VectorFieldElementT<IAInterval, SpaceDimension, TimeDimension>;

#endif


template<typename TT, int sDim, int tDim>
VectorFieldFaceElementT<TT, sDim, tDim>::VectorFieldFaceElementT(const VectorMatrixBase &base,
                                                                 const Cell &c, int face) :
    IFaceElementT<TT, sDim, tDim>(),
    VectorFieldElementT<TT, sDim, tDim>(base, c, face, this->values()) {
  this->ResizeFaceValues(this->nQ(), this->size());
  for (int q = 0; q < this->nQ(); ++q) {
    for (int i = 0; i < this->size(); ++i) {
      this->SetFaceValue(q, i, this->shape(this->QLocal(q), i));
      this->gradient[q][i] =
          this->GetTransformation(q) * this->shape.LocalGradient(this->QLocal(q), i);
    }
  }
}

template class VectorFieldFaceElementT<>;

#ifdef BUILD_IA

template class VectorFieldFaceElementT<IAInterval, SpaceDimension, TimeDimension>;

#endif