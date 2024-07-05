#ifndef RTDISCRETIZATION_HPP_
#define RTDISCRETIZATION_HPP_

#include "IDiscretization.hpp"
#include "RTDoF.hpp"
#include "RTShapes.hpp"

void transformRTData(DataSynchronization &synchronization, const std::vector<int> &indices,
                     const Vector &data);

template<typename T = double, int sDim = SpaceDimension, int tDim = TimeDimension>
class RTDiscretizationT : public NonAdaptiveIDiscretizationT<T, sDim, tDim> {
  int order{};
  int size{};
protected:
  std::unique_ptr<IDoF> createDoF(LevelPair levels) const override {
    return std::make_unique<RTDoF>(order, size);
  }

  std::unique_ptr<ShapeT<T, sDim, tDim>> createShape(CELLTYPE type, int n, int) const override {
    if (n == 0)
      return std::unique_ptr<ShapeT<T, sDim, tDim>>(createRTShape<T, sDim, tDim>(type, order));
    THROW("Shape not implemented in RTDiscretization")
  }

  std::string name(int order, int size) const {
    if (size == 1) return "RT" + std::to_string(order);
    else return "(RT" + std::to_string(order) + ")^" + std::to_string(size);
  }
public:
  RTDiscretizationT(const Meshes &meshes, int order, int size = 1) :
      RTDiscretizationT(meshes, order, 2 * (order + 1), size) {}

  RTDiscretizationT(const Meshes &meshes, int order, int quadExactUpTo, int size) :
      NonAdaptiveIDiscretizationT<T, sDim, tDim>(meshes, name(order, size), quadExactUpTo),
      order(order), size(size) {}

  RTDiscretizationT(const NonAdaptiveIDiscretizationT<T, sDim, tDim> &disc, int order,
                    int size = 1) :
      NonAdaptiveIDiscretizationT<T, sDim, tDim>(disc, name(order, size)), order(order),
      size(size) {}

  void transformData(DataSynchronization &synchronization, const std::vector<int> &indices,
                     const Vector &data) const override {
    transformRTData(synchronization, indices, data);
  }
};

using RTDiscretization = RTDiscretizationT<>;

#ifdef BUILD_IA

using IARTDiscretization = RTDiscretizationT<IAInterval, SpaceDimension, TimeDimension>;

#endif

#endif // RTDISCRETIZATION_HPP
