#include "MixedEGElement.hpp"

template<typename TT, int sDim, int tDim>
MixedEGElementT<TT, sDim, tDim>::MixedEGElementT(const VectorMatrixBase &base, const Cell &c) :
    MixedScalarElementT<TT, sDim, tDim>(base, c) {}

template<typename TT, int sDim, int tDim>
MixedEGElementT<TT, sDim, tDim>::MixedEGElementT(const VectorMatrixBase &base, const Cell &c,
                                                 int face,
                                                 const vector<ShapeValues<TT>> &faceValues) :
    MixedScalarElementT<TT, sDim, tDim>(base, c, face, faceValues) {}

template<typename TT, int sDim, int tDim>
TT MixedEGElementT<TT, sDim, tDim>::PenaltyValue(int q) const {
  return TT(1.0);
}

template<typename TT, int sDim, int tDim>
TT MixedEGElementT<TT, sDim, tDim>::PenaltyValue(int q, const Vector &u) const {
  return MixedScalarElementT<TT, sDim, tDim>::Value(q, u, 1);
}

template<typename TT, int sDim, int tDim>
VectorFieldT<TT, sDim> MixedEGElementT<TT, sDim, tDim>::PenaltyGradient(int q) const {
  return zero;
}

template<typename TT, int sDim, int tDim>
VectorFieldT<TT, sDim> MixedEGElementT<TT, sDim, tDim>::PenaltyGradient(int q,
                                                                        const Vector &u) const {
  return MixedScalarElementT<TT, sDim, tDim>::Derivative(q, u, 1);
}

template<typename TT, int sDim, int tDim>
TT MixedEGElementT<TT, sDim, tDim>::Value(int q, int i) const {
  return MixedScalarElementT<TT, sDim, tDim>::Value(q, 0, i);
}

template<typename TT, int sDim, int tDim>
TT MixedEGElementT<TT, sDim, tDim>::Value(int q, const Vector &u) const {
  return MixedScalarElementT<TT, sDim, tDim>::Value(q, u, 0)
         + MixedScalarElementT<TT, sDim, tDim>::Value(q, u, 1);
}

template<typename TT, int sDim, int tDim>
TT MixedEGElementT<TT, sDim, tDim>::Value(const PointT<TT, sDim, tDim> &z, const Vector &u) const {
  return MixedScalarElementT<TT, sDim, tDim>::Value(z, u, 0)
         + MixedScalarElementT<TT, sDim, tDim>::Value(z, u, 1);
}

template<typename TT, int sDim, int tDim>
VectorFieldT<TT, sDim> MixedEGElementT<TT, sDim, tDim>::Derivative(int q, int i) const {
  return MixedScalarElementT<TT, sDim, tDim>::Derivative(q, 0, i);
}

template<typename TT, int sDim, int tDim>
VectorFieldT<TT, sDim> MixedEGElementT<TT, sDim, tDim>::Derivative(int q, const Vector &u) const {
  return MixedScalarElementT<TT, sDim, tDim>::Derivative(q, u, 0)
         + MixedScalarElementT<TT, sDim, tDim>::Derivative(q, u, 1);
}

template<typename TT, int sDim, int tDim>
VectorFieldT<TT, sDim> MixedEGElementT<TT, sDim, tDim>::Derivative(const PointT<TT, sDim, tDim> &z,
                                                                   const Vector &u) const {
  return MixedScalarElementT<TT, sDim, tDim>::Derivative(z, u, 0);
}


template class MixedEGElementT<>;

#ifdef BUILD_IA

template class MixedEGElementT<IAInterval, SpaceDimension, TimeDimension>;

#endif


template<typename TT, int sDim, int tDim>
MixedEGFaceElementT<TT, sDim, tDim>::MixedEGFaceElementT(const VectorMatrixBase &base,
                                                         const Cell &c, int face) :
    IMixedFaceElementT<TT, sDim, tDim>(base.GetShapeDoF().NumberOfDoFs()),
    MixedEGElementT<TT, sDim, tDim>(base, c, face, this->values()) {
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

template<typename TT, int sDim, int tDim>
int MixedEGFaceElementT<TT, sDim, tDim>::FindQPointID(
    const MixedEGFaceElementT<TT, sDim, tDim> &otherFaceElem,
    const PointT<TT, sDim, tDim> &Qf_c) const {
  for (int q = 0; q < this->nQ(); ++q)
    if (Qf_c == otherFaceElem.QPoint(q)) return q;
  THROW("Error: no qPoint found")
}

template class MixedEGFaceElementT<>;

#ifdef BUILD_IA

template class MixedEGFaceElementT<IAInterval, SpaceDimension, TimeDimension>;

#endif
