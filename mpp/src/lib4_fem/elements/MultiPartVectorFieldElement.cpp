#include "MultiPartVectorFieldElement.hpp"

template<typename TT, int sDim, int tDim>
MultiPartVectorFieldElementT<TT, sDim, tDim>::MultiPartVectorFieldElementT(
    const VectorMatrixBase &base, const Cell &c, int index) :
    MultiPartScalarElementT<TT, sDim, tDim>(base, c, index), dim(c.dim()) {
  this->discSize = dim;
}

template<typename TT, int sDim, int tDim>
MultiPartVectorFieldElementT<TT, sDim, tDim>::MultiPartVectorFieldElementT(
    const VectorMatrixBase &base, const BasicElementT<TT, sDim, tDim> &baseElement, int index) :
    MultiPartScalarElementT<TT, sDim, tDim>(base, baseElement, index), dim(this->c.dim()) {
  this->discSize = dim;
}

template<typename TT, int sDim, int tDim>
VectorFieldComponentT<TT>
MultiPartVectorFieldElementT<TT, sDim, tDim>::VectorComponentValue(int q, int i, int k) const {
  return VectorFieldComponentT<TT>(k, this->Value(q, i));
}

template<typename TT, int sDim, int tDim>
VectorFieldT<TT, sDim, tDim>
MultiPartVectorFieldElementT<TT, sDim, tDim>::VectorValue(int q, const Vector &u) const {
  VectorFieldT<TT, sDim> V;
  for (int i = 0; i < this->size(); ++i) {
    TT s = this->Value(q, i);
    for (int k = 0; k < dim; ++k)
      V[k] += s * u(this->r(i), this->rowIndex(k, this->n(i)));
  }
  return V;
}

template<typename TT, int sDim, int tDim>
VectorFieldT<TT, sDim, tDim>
MultiPartVectorFieldElementT<TT, sDim, tDim>::VectorValue(const PointT<TT, sDim, tDim> &z,
                                                          const Vector &u) const {
  VectorFieldT<TT, sDim> V;
  for (int k = 0; k < dim; ++k)
    V[k] = this->Value(z, u, k);
  return V;
}

template<typename TT, int sDim, int tDim>
VectorFieldT<TT, sDim, tDim>
MultiPartVectorFieldElementT<TT, sDim, tDim>::VectorValue(const PointT<TT, sDim, tDim> &z,
                                                          int i) const {
  VectorFieldT<TT, sDim> V;
  return VectorFieldT<TT, sDim, tDim>(this->Value(z, i));
}

template<typename TT, int sDim, int tDim>
TensorRowT<TT, sDim> MultiPartVectorFieldElementT<TT, sDim, tDim>::VectorRowGradient(int q, int i,
                                                                                     int k) const {
  return TensorRowT<TT, sDim>(k, this->Derivative(q, i));
}

template<typename TT, int sDim, int tDim>
TensorT<TT, sDim>
MultiPartVectorFieldElementT<TT, sDim, tDim>::VectorGradient(int q, const Vector &u) const {
  TensorT<TT, sDim> T;
  for (int i = 0; i < this->size(); ++i)
    for (int k = 0; k < dim; ++k) {
      VectorFieldT<TT, sDim> G =
          u(this->r(i), this->rowIndex(k, this->n(i))) * this->Derivative(q, i);
      T[k] += G;
    }
  return T;
}

template<typename TT, int sDim, int tDim>
TensorT<TT, sDim>
MultiPartVectorFieldElementT<TT, sDim, tDim>::VectorGradient(const PointT<TT, sDim, tDim> &z,
                                                             const Vector &u) const {
  TensorT<TT, sDim> T{};
  for (int k = 0; k < dim; ++k) {
    T[k] += this->Derivative(z, u, k);
  }
  return T;
}


template class MultiPartVectorFieldElementT<>;

#ifdef BUILD_IA

template class MultiPartVectorFieldElementT<IAInterval, SpaceDimension, TimeDimension>;

#endif