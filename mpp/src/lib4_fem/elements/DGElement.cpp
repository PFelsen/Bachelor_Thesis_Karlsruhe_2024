#include "DGElement.hpp"

template<typename TT, int sDim, int tDim>
void DGElementT<TT, sDim, tDim>::init() {
  this->shape.NodalPoints(this->c, nodalPoints);

  for (int q = 0; q < this->nQ(); ++q) {
    for (int i = 0; i < this->shape.size(); ++i) {
      gradient[q][i] = this->GetTransformation(q) * this->shape.LocalGradient(q, i);
    }
  }
}

template<typename TT, int sDim, int tDim>
DGElementT<TT, sDim, tDim>::DGElementT(const VectorMatrixBase &base, const Cell &c) :
    ElementT<TT, sDim, tDim>(base, c), value(this->shape.values()),
    gradient(this->nQ(), this->shape.size()), nodalPoints(0) {
  init();
}

template<typename TT, int sDim, int tDim>
DGElementT<TT, sDim, tDim>::DGElementT(const VectorMatrixBase &base,
                                       const BasicElementT<TT, sDim, tDim> &baseElement) :
    ElementT<TT, sDim, tDim>(base, baseElement), value(this->shape.values()),
    gradient(this->nQ(), this->shape.size()), nodalPoints(0) {
  init();
}

template<typename TT, int sDim, int tDim>
DGElementT<TT, sDim, tDim>::DGElementT(const VectorMatrixBase &base, const Cell &c, int face,
                                       const ShapeValues<TT> &faceValues) :
    ElementT<TT, sDim, tDim>(base, c, face), value(faceValues),
    gradient(this->nQ(), this->shape.size()) {}

template<typename TT, int sDim, int tDim>
int DGElementT<TT, sDim, tDim>::shape_size() const {
  return this->shape.size();
}

template<typename TT, int sDim, int tDim>
int DGElementT<TT, sDim, tDim>::NodalPoints() const {
  return nodalPoints.size();
}

template<typename TT, int sDim, int tDim>
const PointT<TT, sDim, tDim> &DGElementT<TT, sDim, tDim>::NodalPoint(int i) const {
  return nodalPoints[i];
}

template<typename TT, int sDim, int tDim>
TT DGElementT<TT, sDim, tDim>::Value(int q, int i) const {
  return value[q][i];
}

template<typename TT, int sDim, int tDim>
TT DGElementT<TT, sDim, tDim>::Value(int q, const Vector &u) const {
  TT p = 0;
  for (int i = 0; i < this->shape.size(); ++i) {
    p += u(this->r(0), i) * value[q][i];
  }
  return p;
}

template<typename TT, int sDim, int tDim>
TT DGElementT<TT, sDim, tDim>::Value(int q, const Vector &u, int k) const {
  TT p = 0.0;
  for (int i = 0; i < this->size(); ++i) {
    p += u(this->r(0), i + k * this->size()) * value[q][i];
  }
  return p;
}

template<typename TT, int sDim, int tDim>
VectorFieldT<TT, sDim> DGElementT<TT, sDim, tDim>::Derivative(int q, int i) const {
  return gradient[q][i];
}

template<typename TT, int sDim, int tDim>
VectorFieldT<TT, sDim> DGElementT<TT, sDim, tDim>::Derivative(int q, const Vector &u) const {
  VectorFieldT<TT, sDim> Du{};
  for (int i = 0; i < this->shape.size(); ++i) {
    Du += u(this->r(0), i) * gradient[q][i];
  }
  return Du;
}

template class DGElementT<>;


#ifdef BUILD_IA

template class DGElementT<IAInterval, SpaceDimension, TimeDimension>;

#endif


template<typename TT, int sDim, int tDim>
DGFaceElementT<TT, sDim, tDim>::DGFaceElementT(const VectorMatrixBase &base, const Cell &c,
                                               int face) :
    IFaceElementT<TT, sDim, tDim>(), DGElementT<TT, sDim, tDim>(base, c, face, this->values()) {
  this->ResizeFaceValues(this->nQ(), this->shape.size());

  this->shape.NodalPoints(this->c, this->nodalPoints);

  for (int q = 0; q < this->nQ(); ++q) {
    for (int i = 0; i < this->shape_size(); ++i) {
      this->SetFaceValue(q, i, this->shape(this->QLocal(q), i));
      this->gradient[q][i] =
          this->GetTransformation(q) * this->shape.LocalGradient(this->QLocal(q), i);
    }
  }
}

template<typename TT, int sDim, int tDim>
int DGFaceElementT<TT, sDim, tDim>::findQPointID(
    const DGFaceElementT<TT, sDim, tDim> &otherFaceElem, const PointT<TT, sDim, tDim> &Qf_c) const {
  for (int q = 0; q < this->nQ(); ++q)
    if (Qf_c == otherFaceElem.QPoint(q)) return q;
  THROW("Error: no qPoint found")
}


template class DGFaceElementT<>;

#ifdef BUILD_IA

template class DGFaceElementT<IAInterval, SpaceDimension, TimeDimension>;

#endif