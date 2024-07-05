#ifndef SCALARELEMENT_H
#define SCALARELEMENT_H

#include "Element.hpp"

template<typename TT = double, int sDim = SpaceDimension, int tDim = TimeDimension>
class ScalarElementT : public ElementT<TT, sDim, tDim> {
public:
  // All components are set to zero
  static void H10BC(Vector &u);

  static void H10BC(Vector &u, const Cell &C);

  // Only k-th component is set to zero
  static void H10BC(int k, Vector &u);

  static void H10BC(int k, Vector &u, const Cell &C);
private:
  void init();
protected:
  const ShapeValues<TT> &value;
  ShapeValues<VectorFieldT<TT, sDim>> gradient{};

  ScalarElementT(const VectorMatrixBase &base, const Cell &c, int face,
                 const ShapeValues<TT> &faceValues);
public:
  ScalarElementT(const VectorMatrixBase &base, const Cell &c);

  ScalarElementT(const VectorMatrixBase &base, const BasicElementT<TT, sDim, tDim> &baseElement);

  TT Value(int q, int i) const { return this->value[q][i]; }

  TT Value(const PointT<TT, sDim, tDim> &z, int i) const { return this->shape(z, i); }

  TT Value(int q, const Vector &u, int k = 0) const;

  TT Value(const PointT<TT, sDim, tDim> &z, const Vector &u, int k = 0) const;

  VectorFieldT<TT, sDim> Derivative(int q, int i) const { return gradient[q][i]; }

  VectorFieldT<TT, sDim> Derivative(const PointT<TT, sDim, tDim> &z, int i) const;

  VectorFieldT<TT, sDim> Derivative(int q, const Vector &u, int k = 0) const;

  VectorFieldT<TT, sDim> Derivative(const PointT<TT, sDim, tDim> &z, const Vector &u,
                                    int k = 0) const;

  friend std::ostream &operator<<(std::ostream &s, const ScalarElementT<TT, sDim, tDim> &elem) {
    return s << elem.shape;
  }
};

using ScalarElement = ScalarElementT<>;

#ifdef BUILD_IA

using IAScalarElement = ScalarElementT<IAInterval, SpaceDimension, TimeDimension>;

#endif


template<typename TT = double, int sDim = SpaceDimension, int tDim = TimeDimension>
class ScalarFaceElementT :
    public IFaceElementT<TT, sDim, tDim>,
    public ScalarElementT<TT, sDim, tDim> {
public:
  ScalarFaceElementT(const VectorMatrixBase &base, const Cell &c, int face);
};

using ScalarFaceElement = ScalarFaceElementT<>;

#ifdef BUILD_IA

using IAScalarFaceElement = ScalarFaceElementT<IAInterval, SpaceDimension, TimeDimension>;

#endif

#endif // SCALARELEMENT_H
