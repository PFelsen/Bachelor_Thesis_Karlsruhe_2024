#include "DGVectorFieldElement.hpp"

template<typename TT, int sDim, int tDim>
DGVectorFieldElementT<TT, sDim, tDim>::DGVectorFieldElementT(const VectorMatrixBase &base,
                                                             const Cell &c, int face,
                                                             const ShapeValues<TT> &faceValues) :
    DGElementT<TT, sDim, tDim>(base, c, face, faceValues), dim(c.dim()) {}

template<typename TT, int sDim, int tDim>
DGVectorFieldElementT<TT, sDim, tDim>::DGVectorFieldElementT(const VectorMatrixBase &base,
                                                             const Cell &c) :
    DGElementT<TT, sDim, tDim>(base, c), dim(c.dim()) {}

template<typename TT, int sDim, int tDim>
DGVectorFieldElementT<TT, sDim, tDim>::DGVectorFieldElementT(
    const VectorMatrixBase &base, const BasicElementT<TT, sDim, tDim> &baseElement) :
    DGElementT<TT, sDim, tDim>(base, baseElement), dim(this->c.dim()) {}

template<typename TT, int sDim, int tDim>
VectorFieldComponentT<TT> DGVectorFieldElementT<TT, sDim, tDim>::VectorComponentValue(int q, int i,
                                                                                      int k) const {
  return VectorFieldComponentT<TT>(k, this->Value(q, i));
}

template<typename TT, int sDim, int tDim>
VectorFieldT<TT, sDim> DGVectorFieldElementT<TT, sDim, tDim>::VectorValue(int q,
                                                                          const Vector &u) const {
  VectorFieldT<TT, sDim> V;
  for (int i = 0; i < this->shape_size(); ++i) {
    TT s = this->Value(q, i);
    for (int k = 0; k < dim; ++k) {
      V[k] += s * u(this->r(0), k * this->shape_size() + i);
    }
  }
  return V;
}

template<typename TT, int sDim, int tDim>
VectorFieldT<TT, sDim>
DGVectorFieldElementT<TT, sDim, tDim>::VectorValue(const PointT<TT, sDim, tDim> &z,
                                                   const Vector &u) const {
  VectorFieldT<TT, sDim> V;
  for (int k = 0; k < dim; ++k) {
    for (int i = 0; i < this->shape_size(); ++i) {
      V[k] += u(this->r(0), k * this->shape_size() + i) * this->shape(z, i);
    }
  }
  return V;
}

template<typename TT, int sDim, int tDim>
TensorRowT<TT, sDim> DGVectorFieldElementT<TT, sDim, tDim>::VectorRowGradient(int q, int i,
                                                                              int k) const {
  return TensorRowT<TT, sDim>(k, this->Derivative(q, i));
}

template<typename TT, int sDim, int tDim>
TensorT<TT, sDim> DGVectorFieldElementT<TT, sDim, tDim>::VectorGradient(int q,
                                                                        const Vector &u) const {
  TensorT<TT, sDim> T;
  for (int i = 0; i < this->shape_size(); ++i)
    for (int k = 0; k < dim; ++k) {
      VectorFieldT<TT, sDim> G = u(this->r(0), k * this->shape_size() + i) * this->Derivative(q, i);
      for (int l = 0; l < dim; ++l)
        T[k][l] += G[l];
    }
  return T;
}

template<typename TT, int sDim, int tDim>
TensorT<TT, sDim>
DGVectorFieldElementT<TT, sDim, tDim>::VectorGradient(const PointT<TT, sDim, tDim> &z,
                                                      const Vector &u) const {
  TensorT<TT, sDim> T;
  const auto &trans = this->GetCell().GetTransformation(z);
  for (int k = 0; k < dim; ++k) {
    VectorFieldT<TT, sDim> G{};
    for (int i = 0; i < (this->shape).size(); ++i) {
      G +=
          u(this->r(0), k * (this->shape).size() + i) * (trans * (this->shape).LocalGradient(z, i));
      // u(this->r(i), k * this->shape.size() + i) * (this->GetTransformation(z)
      //                                             * (this->shape).LocalGradient(z, i));
    }
    for (int l = 0; l < dim; ++l) {
      T[k][l] += G[l];
    }
  }
  return T;
}

template<typename TT, int sDim, int tDim>
TT DGVectorFieldElementT<TT, sDim, tDim>::Divergence(int q, int i, int k) const {
  VectorFieldT<TT, sDim> G = this->Derivative(q, i);
  return G[k];
}

template<typename TT, int sDim, int tDim>
TT DGVectorFieldElementT<TT, sDim, tDim>::Divergence(int q, const Vector &u) const {
  TT s = TT(0.0);
  for (int i = 0; i < this->size(); ++i) {
    VectorFieldT<TT, sDim> G = this->Derivative(q, i);
    for (int k = 0; k < dim; ++k)
      s += u(this->r(0), k * (this->shape).size() + i) * G[k];
  }
  return s;
}


template class DGVectorFieldElementT<>;

#ifdef BUILD_IA

template class DGVectorFieldElementT<IAInterval, SpaceDimension, TimeDimension>;

#endif

template<typename TT, int sDim, int tDim>
DGVectorFieldFaceElementT<TT, sDim, tDim>::DGVectorFieldFaceElementT(const VectorMatrixBase &base,
                                                                     const Cell &c, int face) :
    IFaceElementT<TT, sDim, tDim>(),
    DGVectorFieldElementT<TT, sDim, tDim>(base, c, face, this->values()) {
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
int DGVectorFieldFaceElementT<TT, sDim, tDim>::findQPointID(
    const DGVectorFieldFaceElementT<TT, sDim, tDim> &otherFaceElem,
    const PointT<TT, sDim, tDim> &Qf_c) const {
  for (int q = 0; q < this->nQ(); ++q)
    if (Qf_c == otherFaceElem.QPoint(q)) return q;
  THROW("Error: no qPoint found")
}

template class DGVectorFieldFaceElementT<>;

#ifdef BUILD_IA

template class DGVectorFieldFaceElementT<IAInterval, SpaceDimension, TimeDimension>;

#endif