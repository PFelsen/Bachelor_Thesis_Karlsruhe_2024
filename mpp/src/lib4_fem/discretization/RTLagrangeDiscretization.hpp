#ifndef RTLAGRANGEDISCRETIZATION_HPP
#define RTLAGRANGEDISCRETIZATION_HPP

#include "BasicDoFs.hpp"
#include "HybridMatrixGraph.hpp"
#include "IDiscretization.hpp"
#include "LagrangeDoF.hpp"
#include "LagrangeShapes.hpp"
#include "Meshes.hpp"
#include "MixedDoF.hpp"
#include "RTDoF.hpp"
#include "RTShapes.hpp"

void transformRTLagrangeData(DataSynchronization &synchronization, const std::vector<int> &indices,
                             const Vector &data, int degree);

template<typename T = double, int sDim = SpaceDimension, int tDim = TimeDimension>
class RTLagrangeDiscretizationT : public NonAdaptiveIDiscretizationT<T, sDim, tDim> {
protected:
  int order{};
  int degree{};
  int size_rt{};
  int size_lagrange{};

  virtual std::unique_ptr<IDoF> createDoF(LevelPair levels) const override {
    std::vector<std::unique_ptr<IDoF>> dofs(2);
    dofs[0] = std::make_unique<RTDoF>(order, size_rt);
    dofs[1] = std::make_unique<LagrangeDoF>(degree, size_lagrange);
    return std::make_unique<MixedDoF>(std::move(dofs));
  }

  std::unique_ptr<ShapeT<T, sDim, tDim>> createShape(CELLTYPE type, int n, int) const override {
    switch (n) {
    case 0:
      return std::unique_ptr<ShapeT<T, sDim, tDim>>(createRTShape<T, sDim, tDim>(type, order));
    case 1:
      return std::unique_ptr<ShapeT<T, sDim, tDim>>(
          createLagrangeShape<T, sDim, tDim>(type, degree));
    default:
      THROW("Shape not implemented in RTLagrangeDiscretization")
    }
  }

  std::string name(int order, int degree, int size_rt, int size_lagrange) const {
    std::string name = "";
    if (size_rt == 1) name.append("RT" + std::to_string(order));
    else name.append("(RT" + std::to_string(order) + ")^" + std::to_string(size_rt));
    if (size_lagrange == 1) name.append("_P" + std::to_string(degree));
    else name.append("_(P" + std::to_string(degree) + ")^" + std::to_string(size_lagrange));
    return name;
  }
public:
  RTLagrangeDiscretizationT(const Meshes &meshes, int order, int degree) :
      RTLagrangeDiscretizationT(meshes, order, degree, std::max(2 * (order + 1), 2 * degree)) {}

  RTLagrangeDiscretizationT(const Meshes &meshes, int order, int degree, int size_rt,
                            int size_lagrange) :
      RTLagrangeDiscretizationT(meshes, order, degree, size_rt, size_lagrange,
                                std::max(2 * (order + 1), 2 * degree)) {}

  RTLagrangeDiscretizationT(const Meshes &meshes, int order, int degree, int quadExactUpTo) :
      RTLagrangeDiscretizationT(meshes, order, degree, 1, 1, quadExactUpTo) {}

  RTLagrangeDiscretizationT(const Meshes &meshes, int order, int degree, int size_rt,
                            int size_lagrange, int quadExactUpTo) :
      NonAdaptiveIDiscretizationT<T, sDim, tDim>(meshes,
                                                 name(order, degree, size_rt, size_lagrange),
                                                 quadExactUpTo),
      order(order), degree(degree), size_rt(size_rt), size_lagrange(size_lagrange) {}

  RTLagrangeDiscretizationT(const NonAdaptiveIDiscretizationT<T, sDim, tDim> &disc, int order,
                            int degree) : RTLagrangeDiscretizationT(disc, order, degree, 1, 1) {}

  RTLagrangeDiscretizationT(const NonAdaptiveIDiscretizationT<T, sDim, tDim> &disc, int order,
                            int degree, int size_rt, int size_lagrange) :
      NonAdaptiveIDiscretizationT<T, sDim, tDim>(disc, name(order, degree, size_rt, size_lagrange)),
      order(order), degree(degree), size_rt(size_rt), size_lagrange(size_lagrange) {}

  int RTOrder() const { return order; }

  int RTSize() const { return size_rt; }

  int LagrangeDegree() const { return degree; }

  int LagrangeSize() const { return size_lagrange; }

  void transformData(DataSynchronization &synchronization, const std::vector<int> &indices,
                     const Vector &data) const override {
    transformRTLagrangeData(synchronization, indices, data, degree);
  }
};

using RTLagrangeDiscretization = RTLagrangeDiscretizationT<>;

#ifdef BUILD_IA

using IARTLagrangeDiscretization =
    RTLagrangeDiscretizationT<IAInterval, SpaceDimension, TimeDimension>;

#endif

template<typename T = double, int sDim = SpaceDimension, int tDim = TimeDimension>
class HybridRTLagrangeDiscretizationT : public RTLagrangeDiscretizationT<T, sDim, tDim> {
protected:
  std::unique_ptr<IMatrixGraph> createMatrixGraph(LevelPair levels) const override {
    return std::make_unique<HybridMatrixGraph>(this->meshes[levels], std::make_unique<FaceDoF>(1),
                                               std::move(this->createDoF(levels)));
  }
public:
  HybridRTLagrangeDiscretizationT(const Meshes &meshes, int order, int degree) :
      RTLagrangeDiscretizationT<T, sDim, tDim>(meshes, order, degree) {
    this->discName = "HybridRTLagrangeDiscretization";
  }

  HybridRTLagrangeDiscretizationT(const NonAdaptiveIDiscretizationT<T, sDim, tDim> &disc, int order,
                                  int degree) :
      RTLagrangeDiscretizationT<T, sDim, tDim>(disc, order, degree) {
    this->discName = "HybridRTLagrangeDiscretization";
  }
};

using HybridRTLagrangeDiscretization = HybridRTLagrangeDiscretizationT<>;

#endif // RTLAGRANGEDISCRETIZATION_HPP
