#ifndef THSERENDIPITYDISCRETIZATION_HPP
#define THSERENDIPITYDISCRETIZATION_HPP

#include "MixedSerendipityDiscretization.hpp"

template<typename T = double, int sDim = SpaceDimension, int tDim = TimeDimension>
class THSerendipityDiscretizationT : public MixedSerendipityDiscretizationT<T, sDim, tDim> {
public:
  explicit THSerendipityDiscretizationT(const Meshes &meshes, int degree = 2, int size = 1) :
      MixedSerendipityDiscretizationT<T, sDim, tDim>(meshes, degree, degree - 1, size) {
    IDiscretizationT<T, sDim, tDim>::discName = this->name(degree, degree - 1, "TaylorHood");
  }

  THSerendipityDiscretizationT(const Meshes &meshes, int degree, int quadExactUpTo, int size) :
      MixedSerendipityDiscretizationT<T, sDim, tDim>(meshes, degree, degree - 1, quadExactUpTo,
                                                     size) {
    IDiscretizationT<T, sDim, tDim>::discName = this->name(degree, degree - 1, "TaylorHood");
  }

  THSerendipityDiscretizationT(const IDiscretizationT<T, sDim, tDim> &disc, int dim, int degree = 2,
                               int size = 1) :
      MixedSerendipityDiscretizationT<T, sDim, tDim>(disc, dim, degree, degree - 1, size) {
    IDiscretizationT<T, sDim, tDim>::discName = this->name(degree, degree - 1, "TaylorHood");
  }
};

using THSerendipityDiscretization = THSerendipityDiscretizationT<>;

#ifdef BUILD_IA

using IATHSerendipityDiscretization =
    THSerendipityDiscretizationT<IAInterval, SpaceDimension, TimeDimension>;

#endif

#endif // SERENDIPITYTHDISCRETIZATION_HPP
