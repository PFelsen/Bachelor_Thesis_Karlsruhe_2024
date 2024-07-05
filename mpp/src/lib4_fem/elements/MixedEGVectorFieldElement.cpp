#include "MixedEGVectorFieldElement.hpp"

template<typename TT, int sDim, int tDim>
MixedEGVectorFieldElementT<TT, sDim, tDim>::MixedEGVectorFieldElementT(
    const VectorMatrixBase &base, const Cell &c, int face,
    const vector<ShapeValues<TT>> &faceValues) :
    MixedVectorFieldElementT<TT, sDim, tDim>(base, c, face, faceValues) {}

template<typename TT, int sDim, int tDim>
VectorFieldT<TT, sDim> MixedEGVectorFieldElementT<TT, sDim, tDim>::PenaltyVectorValue(int q) const {
  return VectorFieldT<TT, sDim>(this->QPoint(q) - this->GetCell().Center());
}

template<typename TT, int sDim, int tDim>

VectorFieldT<TT, sDim>
MixedEGVectorFieldElementT<TT, sDim, tDim>::PenaltyVectorValue(int q, const Vector &u) const {
  return MixedVectorFieldElementT<TT, sDim, tDim>::VectorValue(q, u, 1);
}

template<typename TT, int sDim, int tDim>
TensorT<TT, sDim>
MixedEGVectorFieldElementT<TT, sDim, tDim>::PenaltyVectorGradient(int q, const Vector &u) const {
  return MixedVectorFieldElementT<TT, sDim, tDim>::VectorGradient(q, u, 1);
}

template<typename TT, int sDim, int tDim>
VectorFieldComponentT<TT>
MixedEGVectorFieldElementT<TT, sDim, tDim>::VectorComponentValue(int q, int i, int k) const {
  return MixedVectorFieldElementT<TT, sDim, tDim>::VectorComponentValue(q, 0, i, k);
}

template<typename TT, int sDim, int tDim>

VectorFieldT<TT, sDim>
MixedEGVectorFieldElementT<TT, sDim, tDim>::VectorValue(int q, const Vector &u) const {
  return MixedVectorFieldElementT<TT, sDim, tDim>::VectorValue(q, u, 0)
         + u(this->r(this->nodalPointData[1][0].index), this->nodalPointData[1][0].shift)
               * (this->QPoint(q) - this->GetCell().Center());
}

template<typename TT, int sDim, int tDim>
VectorFieldT<TT, sDim>
MixedEGVectorFieldElementT<TT, sDim, tDim>::VectorValue(const PointT<TT, sDim, tDim> &z,
                                                        const Vector &u) const {
  return MixedVectorFieldElementT<TT, sDim, tDim>::VectorValue(z, u, 0);
}

template<typename TT, int sDim, int tDim>
TensorRowT<TT, sDim> MixedEGVectorFieldElementT<TT, sDim, tDim>::VectorRowGradient(int q, int i,
                                                                                   int k) const {
  return MixedVectorFieldElementT<TT, sDim, tDim>::VectorRowGradient(q, 0, i, k);
}

template<typename TT, int sDim, int tDim>
TensorT<TT, sDim>
MixedEGVectorFieldElementT<TT, sDim, tDim>::VectorGradient(int q, const Vector &u) const {
  return MixedVectorFieldElementT<TT, sDim, tDim>::VectorGradient(q, u, 0)
         + u(this->r(this->nodalPointData[1][0].index), this->nodalPointData[1][0].shift) * One;
}

template<typename TT, int sDim, int tDim>
TensorT<TT, sDim>
MixedEGVectorFieldElementT<TT, sDim, tDim>::VectorGradient(const PointT<TT, sDim, tDim> &z,
                                                           const Vector &u) const {
  return MixedVectorFieldElementT<TT, sDim, tDim>::VectorGradient(z, u, 0);
}


template class MixedEGVectorFieldElementT<>;

#ifdef BUILD_IA

template class MixedEGVectorFieldElementT<IAInterval, SpaceDimension, TimeDimension>;

#endif


template<typename TT, int sDim, int tDim>
MixedEGVectorFieldFaceElementT<TT, sDim, tDim>::MixedEGVectorFieldFaceElementT(
    const VectorMatrixBase &base, const Cell &c, int face) :
    IMixedFaceElementT<TT, sDim, tDim>(base.GetShapeDoF().NumberOfDoFs()),
    MixedEGVectorFieldElementT<TT, sDim, tDim>(base, c, face, this->values()) {
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

template<typename TT, int sDim, int tDim>
int MixedEGVectorFieldFaceElementT<TT, sDim, tDim>::FindQPointID(
    const MixedEGVectorFieldFaceElementT<TT, sDim, tDim> &otherFaceElem,
    const PointT<TT, sDim, tDim> &Qf_c) const {
  for (int q = 0; q < this->nQ(); ++q)
    if (Qf_c == otherFaceElem.QPoint(q)) return q;
  THROW("Error: no qPoint found")
}

template class MixedEGVectorFieldFaceElementT<>;

#ifdef BUILD_IA

template class MixedEGVectorFieldFaceElementT<IAInterval, SpaceDimension, TimeDimension>;

#endif
