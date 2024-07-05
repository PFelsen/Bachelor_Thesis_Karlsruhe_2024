#ifndef MULTIPARTSCALARWGELEMENT_HPP
#define MULTIPARTSCALARWGELEMENT_HPP


#include "Element.hpp"

template<typename TT = double, int sDim = SpaceDimension, int tDim = TimeDimension>
class MultiPartScalarElementWGT : public ElementT<TT, sDim, tDim> {
  int mpIndex{0};

  int indexPart(int n) const;
public:
  MultiPartScalarElementWGT(const VectorMatrixBase &base, const Cell &c, int index = 0);

  int IndexPart(int i) const;

  TT Value(int q, int i) const;

  TT Value(const PointT<TT, sDim, tDim> &z, int i) const;

  TT Value(int q, const Vector &u, int k = 0) const;

  TT Value(const PointT<TT, sDim, tDim> &z, const Vector &u, int k = 0) const;

  void Values(int q, const Vectors &u, std::vector<TT> &U, int k = 0) const;

  friend std::ostream &operator<<(std::ostream &s,
                                  const MultiPartScalarElementWGT<TT, sDim, tDim> &elem) {
    return s << elem.shape;
  };
};

using MultiPartScalarElementWG = MultiPartScalarElementWGT<>;
#ifdef BUILD_IA

using IAMultiPartScalarElementWG =
    MultiPartScalarElementWGT<IAInterval, SpaceDimension, TimeDimension>;

#endif


#endif // MULTIPARTSCALARELEMENT_HPP
