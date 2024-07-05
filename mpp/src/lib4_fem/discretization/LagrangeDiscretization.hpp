#ifndef LAGRANGEDISCRETIZATIONT_HPP
#define LAGRANGEDISCRETIZATIONT_HPP

#include "IDiscretization.hpp"
#include "LagrangeDisplacement.hpp"
#include "LagrangeDoF.hpp"
#include "LagrangeShapes.hpp"
#include "Meshes.hpp"

template<typename T = double, int sDim = SpaceDimension, int tDim = TimeDimension>
class LagrangeDiscretizationT : public NonAdaptiveIDiscretizationT<T, sDim, tDim> {
  int degree{0};
  int size{};
protected:
  std::unique_ptr<IDoF> createDoF(LevelPair levels) const override {
    return std::make_unique<LagrangeDoF>(degree, size);
  }

  std::unique_ptr<ShapeT<T, sDim, tDim>> createShape(CELLTYPE type, int n, int) const override {
    if (n == 0)
      return std::unique_ptr<ShapeT<T, sDim, tDim>>(
          createLagrangeShape<T, sDim, tDim>(type, degree));
    THROW("Shape not implemented in LagrangeDiscretization")
  }

  std::string name(int degree, int size) const {
    if (size == 1) return "P" + std::to_string(degree);
    else return "(P" + std::to_string(degree) + ")^" + std::to_string(size);
  }
public:
  LagrangeDiscretizationT(const Meshes &meshes, int degree, int size = 1) :
      LagrangeDiscretizationT(meshes, degree, 2 * degree, size) {}

  LagrangeDiscretizationT(const Meshes &meshes, int degree, int quadExactUpTo, int size) :
      NonAdaptiveIDiscretizationT<T, sDim, tDim>(meshes, name(degree, size), quadExactUpTo),
      degree(degree), size(size) {}

  LagrangeDiscretizationT(const NonAdaptiveIDiscretizationT<T, sDim, tDim> &disc, int degree,
                          int size = 1) :
      NonAdaptiveIDiscretizationT<T, sDim, tDim>(disc, name(degree, size)), degree(degree),
      size(size) {}

  int Degree() const { return degree; }

  int Size() const { return size; }

  void transformData(DataSynchronization &synchronization, const std::vector<int> &indices,
                     const Vector &data) const override {
    if (degree > 0) {
      transformGenericPointData(synchronization, indices, data);
    } else {
      transformGenericCellData(synchronization, indices, data);
    }
  }

  void transformDisplacement(MeshSynchronization &synchronization,
                             const Vector &displacement) const override {
    transformLagrangeDisplacement(synchronization, displacement);
  }
};

using LagrangeDiscretization = LagrangeDiscretizationT<>;

#ifdef BUILD_IA

using IALagrangeDiscretization = LagrangeDiscretizationT<IAInterval, SpaceDimension, TimeDimension>;

#endif

#endif // LAGRANGEDISCRETIZATIONT_HPP
