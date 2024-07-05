#ifndef DGACOUSTICELEMENT_HPP
#define DGACOUSTICELEMENT_HPP

#include "DGElement.hpp"

template<typename TT = double, int sDim = SpaceDimension, int tDim = TimeDimension>
class DGAcousticElementT : public DGElementT<TT, sDim, tDim> {
protected:
  int degree;
  int dim;

  mutable Tensor *IT = nullptr;
  mutable bool inverseTensor = false;
  ShapeValues<TT> divergence{};
  ShapeValues<TT> pressure{};

  ShapeValues<VectorFieldT<TT, sDim>> velocity{};
  int numL;

  DGAcousticElementT(const VectorMatrixBase &base, const Cell &c, int f_id, int _numL,
                     const ShapeValues<TT> &faceValues);
public:
  DGAcousticElementT(const VectorMatrixBase &base, const Cell &c, int _numL = 0);

  int getDim() const { return dim; }

  int V_dimension() const { return dim * this->shape.size(); }

  int dimension() const { return (dim + 1 + numL) * this->shape.size(); }

  int u_dimension() const { return dim + p_dimension(); }

  int p_dimension() const { return 1 + numL; }

  const VectorField &Gradpressure(int q, int i) const { return this->gradient[q][i]; }

  TT Divergence(int q, int i) const { return divergence[q][i]; }

  TT Pressure(int q, int i) const { return this->value[q][i]; }

  const VectorFieldT<TT, sDim> &Velocity(int q, int i) const { return velocity[q][i]; }

  VectorFieldT<TT, sDim> Velocity(int q, const Vector &u) const;

  TT Pressure(int q, const Vector &u) const;

  TT Pressure_i(int q, const Vector &u, int j) const;

  TT SVValue(const PointT<TT, sDim, tDim> &lz, const Vector &u, int k) const;

  TT evaluateShape(const PointT<TT, sDim, tDim> &lz, int i) const { return this->shape(lz, i); }

  VectorFieldT<TT, sDim> evaluateVelocityShape(const PointT<TT, sDim, tDim> &lz, int i) const;

  const row &Row() const { return this->r(0); }

  ~DGAcousticElementT() {}
};

using DGAcousticElement = DGAcousticElementT<>;

// Todo derive from DGFaceElement and cleanup DGFaceElement
template<typename TT = double, int sDim = SpaceDimension, int tDim = TimeDimension>
class DGAcousticFaceElementT :
    public IFaceElementT<TT, sDim, tDim>,
    public DGAcousticElementT<TT, sDim, tDim> {
public:
  DGAcousticFaceElementT(const VectorMatrixBase &base, const Cell &c, int face, int _numL);
};

using DGAcousticFaceElement = DGAcousticFaceElementT<>;

#endif // DGACOUSTICELEMENT_HPP
