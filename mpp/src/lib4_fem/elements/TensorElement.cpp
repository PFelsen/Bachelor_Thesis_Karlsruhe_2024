#include "TensorElement.hpp"

template<typename TT, int sDim, int tDim>
TensorElementT<TT, sDim, tDim>::TensorElementT(const VectorMatrixBase &base, const Cell &c) :
    VectorFieldElementT<TT, sDim, tDim>(base, c) {}

template<typename TT, int sDim, int tDim>
TensorElementT<TT, sDim, tDim>::TensorElementT(const VectorMatrixBase &base,
                                               const BasicElementT<TT, sDim, tDim> &baseElement) :
    VectorFieldElementT<TT, sDim, tDim>(base, baseElement) {}

template<typename TT, int sDim, int tDim>
TensorComponentT<TT> TensorElementT<TT, sDim, tDim>::TensorComponentValue(int q, int i,
                                                                          int k) const {
  if (k < this->Dim()) return TensorComponentT<TT>(0, k, this->Value(q, i));
  if (k >= 2 * this->Dim()) return TensorComponentT<TT>(2, k - 2 * this->Dim(), this->Value(q, i));
  return TensorComponentT<TT>(1, k - this->Dim(), this->Value(q, i));
}

template<typename TT, int sDim, int tDim>
TensorT<TT, sDim> TensorElementT<TT, sDim, tDim>::TensorValue(int q, const Vector &u) const {
  TensorT<TT, sDim> T;
  if (this->Dim() == 3) {
    for (int i = 0; i < this->size(); ++i) {
      TT s = this->Value(q, i);
      for (int k = 0; k < this->Dim(); ++k) {
        T[0][k] += s * u(this->r(i), k);
        T[1][k] += s * u(this->r(i), k + this->Dim());
        T[2][k] += s * u(this->r(i), k + 2 * this->Dim());
      }
    }
  } else {
    for (int i = 0; i < this->size(); ++i) {
      TT s = this->Value(q, i);
      for (int k = 0; k < this->Dim(); ++k) {
        T[0][k] += s * u(this->r(i), k);
        T[1][k] += s * u(this->r(i), k + this->Dim());
      }
    }
  }
  return T;
}

template<typename TT, int sDim, int tDim>
VectorFieldT<TT, sDim> TensorElementT<TT, sDim, tDim>::DivergenceVector(int q, int i, int k) const {
  if (k < this->Dim()) return VectorFieldT<TT, sDim>(0, this->Divergence(q, i, k));
  if (k >= 2 * this->Dim())
    return VectorFieldT<TT, sDim>(2, this->Divergence(q, i, k - 2 * this->Dim()));
  return VectorFieldT<TT, sDim>(1, this->Divergence(q, i, k - this->Dim()));
}

template<typename TT, int sDim, int tDim>
VectorFieldT<TT, sDim> TensorElementT<TT, sDim, tDim>::DivergenceVector(int q,
                                                                        const Vector &u) const {
  VectorFieldT<TT, sDim> V;
  if (this->Dim() == 3) {
    for (int i = 0; i < this->size(); ++i) {
      VectorFieldT<TT, sDim> G = this->Derivative(q, i);
      for (int k = 0; k < this->Dim(); ++k) {
        V[0] += u(this->r(i), k) * G[k];
        V[1] += u(this->r(i), k + this->Dim()) * G[k];
        V[2] += u(this->r(i), k + 2 * this->Dim()) * G[k];
      }
    }
  } else {
    for (int i = 0; i < this->size(); ++i) {
      VectorFieldT<TT, sDim> G = this->Derivative(q, i);
      for (int k = 0; k < this->Dim(); ++k) {
        V[0] += u(this->r(i), k) * G[k];
        V[1] += u(this->r(i), k + this->Dim()) * G[k];
      }
    }
  }
  return V;
}


template class TensorElementT<>;

#ifdef BUILD_IA

template class TensorElementT<IAInterval, SpaceDimension, TimeDimension>;

#endif