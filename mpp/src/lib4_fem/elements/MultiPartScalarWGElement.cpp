#include "MultiPartScalarWGElement.hpp"

template<typename TT, int sDim, int tDim>
MultiPartScalarElementWGT<TT, sDim, tDim>::MultiPartScalarElementWGT(const VectorMatrixBase &base,
                                                                     const Cell &c, int index) :
    ElementT<TT, sDim, tDim>(base, c), mpIndex(index) {}

template<typename TT, int sDim, int tDim>
int MultiPartScalarElementWGT<TT, sDim, tDim>::indexPart(int n) const {
  return n > 1 ? mpIndex : 0;
}

template<typename TT, int sDim, int tDim>
int MultiPartScalarElementWGT<TT, sDim, tDim>::IndexPart(int i) const {
  return indexPart(this->n(i));
}

template<typename TT, int sDim, int tDim>
TT MultiPartScalarElementWGT<TT, sDim, tDim>::Value(int q, int i) const {
  return this->shape(q, i);
}

template<typename TT, int sDim, int tDim>
TT MultiPartScalarElementWGT<TT, sDim, tDim>::Value(const PointT<TT, sDim, tDim> &z, int i) const {
  return this->shape(z, i);
}

template<typename TT, int sDim, int tDim>
TT MultiPartScalarElementWGT<TT, sDim, tDim>::Value(int q, const Vector &u, int k) const {
  TT U = TT(0.0);
  for (int i = 0; i < this->size(); ++i) {
    int n = this->n(i);
    U += u(this->r(i), k * n + indexPart(n)) * this->shape(q, i);
  }
  return U;
}

template<typename TT, int sDim, int tDim>
void MultiPartScalarElementWGT<TT, sDim, tDim>::Values(int q, const Vectors &u, std::vector<TT> &U,
                                                       int k) const {
  for (int j = 0; j < U.size(); ++j)
    U[j] = 0;
  for (int i = 0; i < this->size(); ++i) {
    int n = this->n(i);
    int m = k * n + indexPart(n);
    TT S = this->shape(q, i);
    for (int j = 0; j < U.size(); ++j)
      U[j] += u[j](this->r(i), m) * S;
  }
}

template<typename TT, int sDim, int tDim>
TT MultiPartScalarElementWGT<TT, sDim, tDim>::Value(const PointT<TT, sDim, tDim> &z,
                                                    const Vector &u, int k) const {
  TT U = TT(0.0);
  for (int i = 0; i < this->size(); ++i) {
    int n = this->n(i);
    U += u(this->r(i), k * n + indexPart(n)) * this->shape(z, i);
  }

  return U;
}
template class MultiPartScalarElementWGT<>;
#ifdef BUILD_IA

template class MultiPartScalarElementWGT<IAInterval, SpaceDimension, TimeDimension>;

#endif
