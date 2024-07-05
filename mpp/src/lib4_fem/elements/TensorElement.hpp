#ifndef TENSORELEMENT_H
#define TENSORELEMENT_H

#include "TensorComponent.hpp"
#include "VectorFieldElement.hpp"

template<typename TT = double, int sDim = SpaceDimension, int tDim = TimeDimension>
class TensorElementT : public VectorFieldElementT<TT, sDim, tDim> {
public:
  TensorElementT(const VectorMatrixBase &base, const Cell &c);

  TensorElementT(const VectorMatrixBase &base, const BasicElementT<TT, sDim, tDim> &baseElement);

  /// Returns the i,k-th nodal function at q-th quadrature point.
  /**
   * If \f$\varphi\f$ is the i-th nodal function of the ScalarElement and n the dimension, then
   * TensorValue returns a matrix with phi as (k/n),(k%n)-th entry, e.g for \f$ n=2,\ k=1 \f$ one
   * gets \f[ V_{i,k}(x,y)=\begin{pmatrix}0&\varphi_i\\0&0\end{pmatrix} \f]
   */

  TensorComponentT<TT> TensorComponentValue(int q, int i, int k) const;

  TensorT<TT, sDim> TensorValue(int q, const Vector &u) const;

  VectorFieldT<TT, sDim> DivergenceVector(int q, int i, int k) const;

  VectorFieldT<TT, sDim> DivergenceVector(int q, const Vector &u) const;
};

using TensorElement = TensorElementT<>;

#ifdef BUILD_IA

using IATensorElement = TensorElementT<IAInterval, SpaceDimension, TimeDimension>;

#endif

#endif // TENSORELEMENT_H
