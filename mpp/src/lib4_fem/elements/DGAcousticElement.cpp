#include "DGAcousticElement.hpp"

template<typename TT, int sDim, int tDim>
DGAcousticElementT<TT, sDim, tDim>::DGAcousticElementT(const VectorMatrixBase &base, const Cell &c,
                                                       int _numL) :
    DGElementT<TT, sDim, tDim>(base, c), dim(c.dim()), degree(base.GetDoF().get_cell_deg(this->c)),
    numL(_numL), velocity(this->nQ(), dim * this->shape.size()),
    pressure(this->nQ(), this->shape.size()), divergence(this->nQ(), dim * this->shape.size()) {

  for (int q = 0; q < this->nQ(); ++q) {
    for (int k = 0; k < dim; ++k) {
      VectorField V = zero;
      V[k] = 1.0;
      for (int i = 0; i < this->shape.size(); ++i) {
        divergence[q][i + k * this->shape.size()] = this->gradient[q][i][k];
        velocity[q][i + k * this->shape.size()] = this->value[q][i] * V;
      }
    }
  }
}

template<typename TT, int sDim, int tDim>
DGAcousticElementT<TT, sDim, tDim>::DGAcousticElementT(const VectorMatrixBase &base, const Cell &c,
                                                       int face, int _numL,
                                                       const ShapeValues<TT> &faceValues) :
    DGElementT<TT, sDim, tDim>(base, c, face, faceValues), dim(c.dim()),
    degree(base.GetDoF().get_cell_deg(this->c)), numL(_numL),
    velocity(this->nQ(), dim * this->shape.size()), pressure(this->nQ(), this->shape.size()),
    divergence(this->nQ(), dim * this->shape.size()) {}

template<typename TT, int sDim, int tDim>
VectorFieldT<TT, sDim> DGAcousticElementT<TT, sDim, tDim>::Velocity(int q, const Vector &u) const {
  VectorFieldT<TT, sDim> V{};
  for (int i = 0; i < V_dimension(); ++i) {
    V += u(this->r(0), i) * velocity[q][i];
  }
  return V;
}

template<typename TT, int sDim, int tDim>
TT DGAcousticElementT<TT, sDim, tDim>::Pressure(int q, const Vector &u) const {
  double p = 0;
  for (int n = 0; n < p_dimension(); ++n)
    p += Pressure_i(q, u, n);
  return p;
}

template<typename TT, int sDim, int tDim>
TT DGAcousticElementT<TT, sDim, tDim>::Pressure_i(int q, const Vector &u, int j) const {
  double p = 0;
  for (int i = 0; i < this->shape.size(); ++i)
    p += u(this->r(0), V_dimension() + j * this->shape.size() + i) * this->value[q][i];
  return p;
}

template<typename TT, int sDim, int tDim>
TT DGAcousticElementT<TT, sDim, tDim>::SVValue(const PointT<TT, sDim, tDim> &lz, const Vector &u,
                                               int k) const {
  TT v{0.0};
  for (int i = 0; i < this->shape.size(); ++i)
    v += u(this->r(0), k * this->shape.size() + i) * this->shape(lz, i);
  return v;
}

template<typename TT, int sDim, int tDim>
VectorFieldT<TT, sDim>
DGAcousticElementT<TT, sDim, tDim>::evaluateVelocityShape(const PointT<TT, sDim, tDim> &lz,
                                                          int i) const {
  int j = i % this->shape.size(), k = i / this->shape.size();
  VectorFieldT<TT, sDim> V(zero);
  V[k] = 1.0;
  return this->shape(lz, j) * V;
}


template class DGAcousticElementT<>;

template<typename TT, int sDim, int tDim>
DGAcousticFaceElementT<TT, sDim, tDim>::DGAcousticFaceElementT(const VectorMatrixBase &base,
                                                               const Cell &c, int face, int _numL) :
    IFaceElementT<TT, sDim, tDim>(),
    DGAcousticElementT<TT, sDim, tDim>(base, c, face, _numL, this->values()) {
  this->ResizeFaceValues(this->nQ(), this->shape.size());
  this->shape.NodalPoints(this->c, this->nodalPoints);

  for (int q = 0; q < this->nQ(); ++q) {
    for (int i = 0; i < this->shape_size(); ++i) {
      this->SetFaceValue(q, i, this->shape(this->QLocal(q), i));
      this->gradient[q][i] =
          this->GetTransformation(q) * this->shape.LocalGradient(this->QLocal(q), i);
      this->pressure[q][i] = this->shape(this->QLocal(q), i);
      for (int k = 0; k < this->dim; ++k) {
        VectorFieldT<TT, sDim> V{};
        V[k] = 1.0;
        for (int i = 0; i < this->shape.size(); ++i) {
          this->velocity[q][i + k * this->shape.size()] = this->pressure[q][i] * V;
        }
      }
    }
  }
}

template class DGAcousticFaceElementT<>;