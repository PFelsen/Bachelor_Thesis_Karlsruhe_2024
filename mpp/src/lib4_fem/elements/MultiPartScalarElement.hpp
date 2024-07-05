#ifndef MULTIPARTSCALARELEMENT_HPP
#define MULTIPARTSCALARELEMENT_HPP

#include "Element.hpp"

template<typename TT = double, int sDim = SpaceDimension, int tDim = TimeDimension>
class MultiPartScalarElementT : public ElementT<TT, sDim, tDim> {
  ShapeValues<VectorFieldT<TT, sDim>> gradient;
protected:
  int discSize{1};
  int mpIndex{0};

  int indexPart(int n) const;
  int rowIndex(int k, int n) const;
public:
  MultiPartScalarElementT(const VectorMatrixBase &base, const Cell &c, int index = 0);

  MultiPartScalarElementT(const VectorMatrixBase &base,
                          const BasicElementT<TT, sDim, tDim> &baseElement, int index = 0);

  int IndexPart(int i) const;

  TT Value(int q, int i) const;

  TT Value(const PointT<TT, sDim, tDim> &z, int i) const;

  TT Value(int q, const Vector &u, int k = 0) const;

  TT Value(const PointT<TT, sDim, tDim> &z, const Vector &u, int k = 0) const;

  void Values(int q, const Vectors &u, std::vector<TT> &U, int k = 0) const;


  VectorFieldT<TT, sDim> Derivative(int q, int i) const;

  VectorFieldT<TT, sDim> Derivative(const PointT<TT, sDim, tDim> &z, int i) const;

  VectorFieldT<TT, sDim> Derivative(int q, const Vector &u, int k = 0) const;

  VectorFieldT<TT, sDim> Derivative(const PointT<TT, sDim, tDim> &z, const Vector &u,
                                    int k = 0) const;

  friend std::ostream &operator<<(std::ostream &s,
                                  const MultiPartScalarElementT<TT, sDim, tDim> &elem) {
    return s << elem.shape;
  };
};

using MultiPartScalarElement = MultiPartScalarElementT<>;
#ifdef BUILD_IA

using IAMultiPartScalarElement = MultiPartScalarElementT<IAInterval, SpaceDimension, TimeDimension>;

#endif

#endif // MULTIPARTSCALARELEMENT_HPP
