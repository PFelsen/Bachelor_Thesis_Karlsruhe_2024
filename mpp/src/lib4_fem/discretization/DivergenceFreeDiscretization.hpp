#ifndef DIVERGENCEFREEDISCRETIZATION_HPP
#define DIVERGENCEFREEDISCRETIZATION_HPP

#include "ArgyrisDoF.hpp"
#include "ArgyrisShapes.hpp"
#include "IDiscretization.hpp"
#include "Meshes.hpp"

void transformDivergenceFreeData(DataSynchronization &synchronization,
                                 const std::vector<int> &indices, const Vector &data);

template<typename T = double, int sDim = SpaceDimension, int tDim = TimeDimension>
class DivergenceFreeDiscretizationT : public NonAdaptiveIDiscretizationT<T, sDim, tDim> {
  int size{};
protected:
  std::unique_ptr<IDoF> createDoF(LevelPair levels) const override {
    return std::make_unique<ArgyrisDoF>(size);
  }

  std::unique_ptr<ShapeT<T, sDim, tDim>> createShape(CELLTYPE type, int n, int) const override {
    if (type != TRIANGLE) THROW("DivergenceFreeDiscretization only implemented for triangles")
    if (n == 0) return std::unique_ptr<ShapeT<T, sDim, tDim>>(new ArgyrisTriT<T, sDim, tDim>());
    THROW("Shape not implemented in DivergenceFreeDiscretization")
  }
public:
  DivergenceFreeDiscretizationT(const Meshes &meshes, int size = 1) :
      DivergenceFreeDiscretizationT(meshes, 2 * 4, size) {}

  DivergenceFreeDiscretizationT(const Meshes &meshes, int quadExactUpTo, int size) :
      NonAdaptiveIDiscretizationT<T, sDim, tDim>(meshes, "DivergenceFree", quadExactUpTo),
      size(size) {}

  DivergenceFreeDiscretizationT(const NonAdaptiveIDiscretizationT<T, sDim, tDim> &disc,
                                int size = 1) :
      NonAdaptiveIDiscretizationT<T, sDim, tDim>(disc, "DivergenceFree"), size(size) {}

  void transformData(DataSynchronization &synchronization, const std::vector<int> &indices,
                     const Vector &data) const override {
    transformDivergenceFreeData(synchronization, indices, data);
  }
};

using DivergenceFreeDiscretization = DivergenceFreeDiscretizationT<>;

#ifdef BUILD_IA

using IADivergenceFreeDiscretization =
    DivergenceFreeDiscretizationT<IAInterval, SpaceDimension, TimeDimension>;

#endif

#endif // DIVERGENCEFREEDISCRETIZATION_HPP
