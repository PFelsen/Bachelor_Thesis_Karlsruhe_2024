#ifndef THDISCRETIZATION_HPP
#define THDISCRETIZATION_HPP

#include "MixedLagrangeDiscretization.hpp"

template<typename T = double, int sDim = SpaceDimension, int tDim = TimeDimension>
class THDiscretizationT : public MixedLagrangeDiscretizationT<T, sDim, tDim> {
public:
  explicit THDiscretizationT(const Meshes &meshes, int degree = 2, int size = 1) :
      MixedLagrangeDiscretizationT<T, sDim, tDim>(meshes, degree, degree - 1, size) {
    IDiscretizationT<T, sDim, tDim>::discName = this->name(degree, degree - 1, "TaylorHood");
  }

  THDiscretizationT(const Meshes &meshes, int degree, int quadExactUpTo, int size) :
      MixedLagrangeDiscretizationT<T, sDim, tDim>(meshes, degree, degree - 1, quadExactUpTo, size) {
    IDiscretizationT<T, sDim, tDim>::discName = this->name(degree, degree - 1, "TaylorHood");
  }

  THDiscretizationT(const IDiscretizationT<T, sDim, tDim> &disc, int dim, int degree = 2,
                    int size = 1) :
      MixedLagrangeDiscretizationT<T, sDim, tDim>(disc, dim, degree, degree - 1, size) {
    IDiscretizationT<T, sDim, tDim>::discName = this->name(degree, degree - 1, "TaylorHood");
  }
};

using THDiscretization = THDiscretizationT<>;

#ifdef BUILD_IA

using IATHDiscretization = THDiscretizationT<IAInterval, SpaceDimension, TimeDimension>;

#endif

#endif // THDISCRETIZATION_HPP
