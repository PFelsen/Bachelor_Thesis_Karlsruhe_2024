#ifndef DGDISCRETIZATION_H
#define DGDISCRETIZATION_H

#include <map>
#include "DGDoF.hpp"
#include "DGMatrixGraph.hpp"
#include "IDiscretization.hpp"
#include "LagrangeShapes.hpp"
#include "Meshes.hpp"

void transformDGDisplacement(MeshSynchronization &synchronization, const Vector &displacement);

template<typename T = double, int sDim = SpaceDimension, int tDim = TimeDimension>
class DGDiscretizationT : public NonAdaptiveIDiscretizationT<T, sDim, tDim> {
  int degree{};
  int size{};
protected:
  std::unique_ptr<IMatrixGraph> createMatrixGraph(LevelPair levels) const override {
    return std::make_unique<DGMatrixGraph>((this->meshes)[levels], createDoF(levels));
  }

  std::unique_ptr<IDoF> createDoF(LevelPair levels) const override {
    return std::make_unique<DGDoF>(degree, size);
  }

  std::unique_ptr<ShapeT<T, sDim, tDim>> createShape(CELLTYPE type, int n, int) const override {
    if (n == 0)
      return std::unique_ptr<ShapeT<T, sDim, tDim>>(
          createLagrangeShape<T, sDim, tDim>(type, degree));
    THROW("Shape not implemented in DGDiscretization")
  }
public:
  DGDiscretizationT(const Meshes &meshes, int degree, int size = 1) :
      DGDiscretizationT(meshes, degree, 2 * degree, size) {}

  DGDiscretizationT(const Meshes &meshes, int degree, int quadExactUpTo, int size) :
      NonAdaptiveIDiscretizationT<T, sDim, tDim>(meshes, "DG_P" + std::to_string(degree),
                                                 quadExactUpTo),
      degree(degree), size(size) {}

  DGDiscretizationT(const NonAdaptiveIDiscretizationT<T, sDim, tDim> &disc, int degree,
                    int size = 1) :
      NonAdaptiveIDiscretizationT<T, sDim, tDim>(disc.GetMeshes(), "DG_P" + std::to_string(degree),
                                                 disc.quadExactUpTo(0, degree)) {}

  int Degree() const { return degree; }

  void transformData(DataSynchronization &synchronization, const std::vector<int> &indices,
                     const Vector &data) const override {
    transformGenericCellData(synchronization, indices, data);
  }

  void transformDisplacement(MeshSynchronization &synchronization,
                             const Vector &displacement) const override {
    transformDGDisplacement(synchronization, displacement);
  }
};

typedef DGDiscretizationT<> DGDiscretization;

#endif // TUTORIAL_DGDISCRETIZATION_H
