#include "ScalarElement.hpp"

template<typename TT, int sDim, int tDim>
void ScalarElementT<TT, sDim, tDim>::H10BC(Vector &u) {
  u.ClearDirichletFlags();
  for (cell c = u.cells(); c != u.cells_end(); ++c) {
    H10BC(u, *c);
  }
  u.DirichletConsistent();
}

template<typename TT, int sDim, int tDim>
void ScalarElementT<TT, sDim, tDim>::H10BC(Vector &u, const Cell &C) {
  if (!u.OnBoundary(C)) return;
  for (int face = 0; face < C.Faces(); ++face) {
    if (!u.OnBoundary(C, face)) continue;
    rows R(u.GetMatrixGraph(), C, face);
    for (const row &r : R) {
      for (int k = 0; k < u.NumberOfDoFs(); ++k) {
        u(r, k) = 0.0;
        u.D(r, k) = true;
      }
    }
  }
}

template<typename TT, int sDim, int tDim>
void ScalarElementT<TT, sDim, tDim>::H10BC(int k, Vector &u) {
  u.ClearDirichletFlags();
  for (cell c = u.cells(); c != u.cells_end(); ++c) {
    H10BC(k, u, *c);
  }
  u.DirichletConsistent();
}

template<typename TT, int sDim, int tDim>
void ScalarElementT<TT, sDim, tDim>::H10BC(int k, Vector &u, const Cell &C) {
  if (!u.OnBoundary(C)) return;
  for (int face = 0; face < C.Faces(); ++face) {
    if (!u.OnBoundary(C, face)) continue;
    rows R(u.GetMatrixGraph(), C, face);
    for (const row &r : R) {
      u(r, k) = 0.0;
      u.D(r, k) = true;
    }
  }
}

template<typename TT, int sDim, int tDim>
void ScalarElementT<TT, sDim, tDim>::init() {
  for (int q = 0; q < this->nQ(); ++q)
    for (int i = 0; i < this->size(); ++i)
      gradient[q][i] = this->GetTransformation(q) * this->shape.LocalGradient(q, i);
}

template<typename TT, int sDim, int tDim>
ScalarElementT<TT, sDim, tDim>::ScalarElementT(const VectorMatrixBase &base, const Cell &c) :
    ElementT<TT, sDim, tDim>(base, c), value(this->shape.values()),
    gradient(this->nQ(), this->shape.size()) {
  init();
}

template<typename TT, int sDim, int tDim>
ScalarElementT<TT, sDim, tDim>::ScalarElementT(const VectorMatrixBase &base,
                                               const BasicElementT<TT, sDim, tDim> &baseElement) :
    ElementT<TT, sDim, tDim>(base, baseElement), value(this->shape.values()),
    gradient(this->nQ(), this->shape.size()) {
  init();
}

template<typename TT, int sDim, int tDim>
ScalarElementT<TT, sDim, tDim>::ScalarElementT(const VectorMatrixBase &base, const Cell &c,
                                               int face, const ShapeValues<TT> &faceValues) :
    ElementT<TT, sDim, tDim>(base, c, face), value(faceValues),
    gradient(this->nQ(), this->shape.size()) {}

template<typename TT, int sDim, int tDim>
TT ScalarElementT<TT, sDim, tDim>::Value(int q, const Vector &u, int k) const {
  TT U = TT(0.0);
  for (int i = 0; i < this->size(); ++i)
    U += u(this->r(i), k) * this->Value(q, i);
  return U;
}

template<typename TT, int sDim, int tDim>
TT ScalarElementT<TT, sDim, tDim>::Value(const PointT<TT, sDim, tDim> &z, const Vector &u,
                                         int k) const {
  TT U = TT(0.0);
  for (int i = 0; i < this->size(); ++i)
    U += u(this->r(i), k) * this->shape(z, i);
  return U;
}

template<typename TT, int sDim, int tDim>
VectorFieldT<TT, sDim> ScalarElementT<TT, sDim, tDim>::Derivative(const PointT<TT, sDim, tDim> &z,
                                                                  int i) const {
  return this->GetTransformation(z) * this->shape.LocalGradient(z, i);
}

template<typename TT, int sDim, int tDim>
VectorFieldT<TT, sDim> ScalarElementT<TT, sDim, tDim>::Derivative(int q, const Vector &u,
                                                                  int k) const {
  VectorFieldT<TT, sDim> Du;
  for (int i = 0; i < this->size(); ++i)
    Du += u(this->r(i), k) * gradient[q][i];
  return Du;
}

template<typename TT, int sDim, int tDim>
VectorFieldT<TT, sDim> ScalarElementT<TT, sDim, tDim>::Derivative(const PointT<TT, sDim, tDim> &z,
                                                                  const Vector &u, int k) const {
  VectorFieldT<TT, sDim> Du;
  const TransformationT<TT, sDim> &T = this->GetTransformation(z);
  for (int i = 0; i < this->size(); ++i)
    Du += u(this->r(i), k) * (T * this->shape.LocalGradient(z, i));
  return Du;
}


template class ScalarElementT<>;

#ifdef BUILD_IA

template class ScalarElementT<IAInterval, SpaceDimension, TimeDimension>;

#endif

template<typename TT, int sDim, int tDim>
ScalarFaceElementT<TT, sDim, tDim>::ScalarFaceElementT(const VectorMatrixBase &base, const Cell &c,
                                                       int face) :
    IFaceElementT<TT, sDim, tDim>(), ScalarElementT<TT, sDim, tDim>(base, c, face, this->values()) {
  this->ResizeFaceValues(this->nQ(), this->size());
  for (int q = 0; q < this->nQ(); ++q) {
    for (int i = 0; i < this->size(); ++i) {
      this->SetFaceValue(q, i, this->shape(this->QLocal(q), i));
      this->gradient[q][i] =
          this->GetTransformation(q) * this->shape.LocalGradient(this->QLocal(q), i);
    }
  }
}

template class ScalarFaceElementT<>;

#ifdef BUILD_IA

template class ScalarFaceElementT<IAInterval, SpaceDimension, TimeDimension>;

#endif