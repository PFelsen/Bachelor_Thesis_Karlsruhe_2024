#ifndef CRDISCRETIZATION_HPP
#define CRDISCRETIZATION_HPP

#include "BasicDoFs.hpp"
#include "CRShapes.hpp"
#include "IDiscretization.hpp"
#include "Meshes.hpp"

template<typename T = double, int sDim = SpaceDimension, int tDim = TimeDimension>
class CRDiscretizationT : public NonAdaptiveIDiscretizationT<T, sDim, tDim> {
  int order{};
  int size{};
protected:
  std::unique_ptr<IDoF> createDoF(LevelPair levels) const override {
    return std::make_unique<FaceDoF>(order, size); // TODO: Check
  }

  std::unique_ptr<ShapeT<T, sDim, tDim>> createShape(CELLTYPE type, int n, int) const override {
    if (n == 0)
      return std::unique_ptr<ShapeT<T, sDim, tDim>>(createCRShape<T, sDim, tDim>(type, order));
    THROW("Shape not implemented in CRDiscretization")
  }
public:
  CRDiscretizationT(const Meshes &meshes, int order, int size = 1) :
      CRDiscretizationT(meshes, order, 2 * order, size) {}

  CRDiscretizationT(const Meshes &meshes, int order, int quadExactUpTo, int size) :
      NonAdaptiveIDiscretizationT<T, sDim, tDim>(meshes, "CR" + std::to_string(order),
                                                 quadExactUpTo),
      order(order), size(size) {}

  CRDiscretizationT(const NonAdaptiveIDiscretizationT<T, sDim, tDim> &disc, int order,
                    int size = 1) :
      NonAdaptiveIDiscretizationT<T, sDim, tDim>(disc, "CR" + std::to_string(order)), order(order),
      size(size) {}
};

using CRDiscretization = CRDiscretizationT<>;

#endif // CRDISCRETIZATION_HPP
