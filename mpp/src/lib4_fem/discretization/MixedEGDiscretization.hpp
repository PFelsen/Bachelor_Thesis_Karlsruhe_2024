#ifndef MIXEDEGDISCRETIZATION_HPP
#define MIXEDEGDISCRETIZATION_HPP

#include "EGMatrixGraph.hpp"
#include "IDiscretization.hpp"
#include "LagrangeDisplacement.hpp"
#include "LagrangeDoF.hpp"
#include "LagrangeShapes.hpp"
#include "MixedDoF.hpp"

template<typename T = double, int sDim = SpaceDimension, int tDim = TimeDimension>
class MixedEGDiscretizationT : public NonAdaptiveIDiscretizationT<T, sDim, tDim> {
  int degree{};
  int size{};
protected:
  std::unique_ptr<IMatrixGraph> createMatrixGraph(LevelPair levels) const override {
    return std::make_unique<EGMatrixGraph>((this->meshes)[levels], createDoF(levels), 1);
  }

  std::unique_ptr<IDoF> createDoF(LevelPair levels) const override {
    std::vector<std::unique_ptr<IDoF>> dofs(2);
    dofs[0] = std::make_unique<LagrangeDoF>(degree, size);
    dofs[1] = std::make_unique<LagrangeDoF>(0, 1);
    return std::make_unique<MixedDoF>(std::move(dofs));
  }

  std::unique_ptr<ShapeT<T, sDim, tDim>> createShape(CELLTYPE type, int n, int) const override {
    if (n == 0)
      return std::unique_ptr<ShapeT<T, sDim, tDim>>(
          createLagrangeShape<T, sDim, tDim>(type, degree));
    if (n == 1)
      return std::unique_ptr<ShapeT<T, sDim, tDim>>(createLagrangeShape<T, sDim, tDim>(type, 0));
    THROW("Shape not implemented in MixedEGDiscretization")
  }
public:
  MixedEGDiscretizationT(const Meshes &meshes, int degree, int size = 1) :
      MixedEGDiscretizationT(meshes, degree, 2 * degree, size) {}

  MixedEGDiscretizationT(const Meshes &meshes, int degree, int quadExactUpTo, int size) :
      NonAdaptiveIDiscretizationT<T, sDim, tDim>(meshes, "EG_P" + std::to_string(degree),
                                                 quadExactUpTo),
      degree(degree), size(size) {}

  MixedEGDiscretizationT(const NonAdaptiveIDiscretizationT<T, sDim, tDim> &disc, int degree,
                         int size = 1) :
      NonAdaptiveIDiscretizationT<T, sDim, tDim>(disc, "EG_P" + std::to_string(degree), degree),
      degree(degree), size(size) {}

  int Degree() const { return degree; }

  int Size() const { return size; }

  void transformData(DataSynchronization &synchronization, const std::vector<int> &indices,
                     const Vector &data) const override {
    //    transformGenericCellData(synchronization, indices, data);
    transformGenericPointData(synchronization, indices, data);
  }

  void transformDisplacement(MeshSynchronization &synchronization,
                             const Vector &displacement) const override {
    transformLagrangeDisplacement(synchronization, displacement);
  }
};

using MixedEGDiscretization = MixedEGDiscretizationT<>;

#ifdef BUILD_IA

using IAMixedEGDiscretization = MixedEGDiscretizationT<IAInterval, SpaceDimension, TimeDimension>;

#endif

#endif // MIXEDEGDISCRETIZATION_HPP
