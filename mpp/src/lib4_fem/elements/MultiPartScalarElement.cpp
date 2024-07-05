#include "MultiPartScalarElement.hpp"

template<typename TT, int sDim, int tDim>
MultiPartScalarElementT<TT, sDim, tDim>::MultiPartScalarElementT(const VectorMatrixBase &base,
                                                                 const Cell &c, int index) :
    ElementT<TT, sDim, tDim>(base, c), gradient(this->nQ(), this->shape.size()), mpIndex(index) {
  for (int q = 0; q < this->nQ(); ++q)
    for (int i = 0; i < this->size(); ++i)
      gradient[q][i] = this->GetTransformation(q) * this->shape.LocalGradient(q, i);
}

template<typename TT, int sDim, int tDim>
MultiPartScalarElementT<TT, sDim, tDim>::MultiPartScalarElementT(
    const VectorMatrixBase &base, const BasicElementT<TT, sDim, tDim> &baseElement, int index) :
    ElementT<TT, sDim, tDim>(base, baseElement), gradient(this->nQ(), this->shape.size()),
    mpIndex(index) {
  for (int q = 0; q < this->nQ(); ++q)
    for (int i = 0; i < this->size(); ++i)
      gradient[q][i] = this->GetTransformation(q) * this->shape.LocalGradient(q, i);
}

template<typename TT, int sDim, int tDim>
int MultiPartScalarElementT<TT, sDim, tDim>::rowIndex(int k, int n) const {
  return (2 - (n <= discSize)) * k + indexPart(n);
}

template<typename TT, int sDim, int tDim>
int MultiPartScalarElementT<TT, sDim, tDim>::indexPart(int n) const {
  return n > discSize ? mpIndex : 0;
}

template<typename TT, int sDim, int tDim>
int MultiPartScalarElementT<TT, sDim, tDim>::IndexPart(int i) const {
  return indexPart(this->n(i));
}

template<typename TT, int sDim, int tDim>
TT MultiPartScalarElementT<TT, sDim, tDim>::Value(int q, int i) const {
  return this->shape(q, i);
}

template<typename TT, int sDim, int tDim>
TT MultiPartScalarElementT<TT, sDim, tDim>::Value(const PointT<TT, sDim, tDim> &z, int i) const {
  return this->shape(z, i);
}

template<typename TT, int sDim, int tDim>
TT MultiPartScalarElementT<TT, sDim, tDim>::Value(int q, const Vector &u, int k) const {
  TT U = TT(0.0);
  for (int i = 0; i < this->size(); ++i) {
    U += u(this->r(i), this->rowIndex(k, this->n(i))) * this->shape(q, i);
  }
  return U;
}

template<typename TT, int sDim, int tDim>
TT MultiPartScalarElementT<TT, sDim, tDim>::Value(const PointT<TT, sDim, tDim> &z, const Vector &u,
                                                  int k) const {
  TT U = TT(0.0);
  for (int i = 0; i < this->size(); ++i) {
    U += u(this->r(i), this->rowIndex(k, this->n(i))) * this->shape(z, i);
  }

  return U;
}

template<typename TT, int sDim, int tDim>
void MultiPartScalarElementT<TT, sDim, tDim>::Values(int q, const Vectors &u, std::vector<TT> &U,
                                                     int k) const {
  for (int j = 0; j < U.size(); ++j)
    U[j] = 0;
  for (int i = 0; i < this->size(); ++i) {
    int m = this->rowIndex(k, this->n(i));
    TT S = this->shape(q, i);
    for (int j = 0; j < U.size(); ++j)
      U[j] += u[j](this->r(i), m) * S;
  }
}

template<typename TT, int sDim, int tDim>
VectorFieldT<TT, sDim> MultiPartScalarElementT<TT, sDim, tDim>::Derivative(int q, int i) const {
  return gradient[q][i];
}

template<typename TT, int sDim, int tDim>
VectorFieldT<TT, sDim>
MultiPartScalarElementT<TT, sDim, tDim>::Derivative(const PointT<TT, sDim, tDim> &z, int i) const {
  return this->GetTransformation(z) * this->shape.LocalGradient(z, i);
}

template<typename TT, int sDim, int tDim>
VectorFieldT<TT, sDim> MultiPartScalarElementT<TT, sDim, tDim>::Derivative(int q, const Vector &u,
                                                                           int k) const {
  VectorFieldT<TT, sDim> Du;
  for (int i = 0; i < this->size(); ++i) {
    Du += u(this->r(i), this->rowIndex(k, this->n(i))) * gradient[q][i];
  }
  return Du;
}

template<typename TT, int sDim, int tDim>
VectorFieldT<TT, sDim>
MultiPartScalarElementT<TT, sDim, tDim>::Derivative(const PointT<TT, sDim, tDim> &z,
                                                    const Vector &u, int k) const {
  VectorFieldT<TT, sDim> Du;
  const TransformationT<TT, sDim> &T = this->GetTransformation(z);
  for (int i = 0; i < this->size(); ++i) {
    Du += u(this->r(i), this->rowIndex(k, this->n(i))) * (T * this->shape.LocalGradient(z, i));
  }
  return Du;
}


template class MultiPartScalarElementT<>;

#ifdef BUILD_IA

template class MultiPartScalarElementT<IAInterval, SpaceDimension, TimeDimension>;

#endif
