#ifndef ARGYRISDISCRETIZATION_HPP
#define ARGYRISDISCRETIZATION_HPP

#include "ArgyrisDoF.hpp"
#include "ArgyrisShapes.hpp"
#include "IDiscretization.hpp"
#include "Meshes.hpp"

void transformArgyrisData(DataSynchronization &synchronization, const std::vector<int> &indices,
                          const Vector &data);

template<typename T = double, int sDim = SpaceDimension, int tDim = TimeDimension>
class ArgyrisDiscretizationT : public NonAdaptiveIDiscretizationT<T, sDim, tDim> {
  int size{};
protected:
  std::unique_ptr<IDoF> createDoF(LevelPair levels) const override {
    return std::make_unique<ArgyrisDoF>(size);
  }

  std::unique_ptr<ShapeT<T, sDim, tDim>> createShape(CELLTYPE type, int n, int) const override {
    if (type != TRIANGLE) THROW("ArgyrisDiscretization only implemented for triangles")
    if (n == 0) return std::unique_ptr<ShapeT<T, sDim, tDim>>(new ArgyrisTriT<T, sDim, tDim>());
    THROW("Shape not implemented in ArgyrisDiscretization")
  }
public:
  ArgyrisDiscretizationT(const Meshes &meshes, int size = 1) :
      ArgyrisDiscretizationT(meshes, 2 * 5, size) {}

  ArgyrisDiscretizationT(const Meshes &meshes, int quadExactUpTo, int size) :
      NonAdaptiveIDiscretizationT<T, sDim, tDim>(meshes, "Argyris", quadExactUpTo), size(size) {}

  ArgyrisDiscretizationT(const NonAdaptiveIDiscretizationT<T, sDim, tDim> &disc, int size = 1) :
      NonAdaptiveIDiscretizationT<T, sDim, tDim>(disc, "Argyris"), size(size) {}

  void transformData(DataSynchronization &synchronization, const std::vector<int> &indices,
                     const Vector &data) const override {
    transformArgyrisData(synchronization, indices, data);
  }
};

using ArgyrisDiscretization = ArgyrisDiscretizationT<>;

#ifdef BUILD_IA

using IAArgyrisDiscretization = ArgyrisDiscretizationT<IAInterval, SpaceDimension, TimeDimension>;

#endif

#endif // ARGYRISDISCRETIZATION_HPP
