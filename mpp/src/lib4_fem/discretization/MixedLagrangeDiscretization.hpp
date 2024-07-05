#ifndef MIXEDLAGRANGEDISCRETIZATION_HPP
#define MIXEDLAGRANGEDISCRETIZATION_HPP

#include "IDiscretization.hpp"
#include "LagrangeDoF.hpp"
#include "LagrangeShapes.hpp"
#include "Meshes.hpp"
#include "MixedDoF.hpp"

template<typename T = double, int sDim = SpaceDimension, int tDim = TimeDimension>
class MixedLagrangeDiscretizationT : public NonAdaptiveIDiscretizationT<T, sDim, tDim> {
  int degree1{};
  int degree2{};
  int dim{};
  int size{};
protected:
  std::string name(int degree1, int degree2, std::string prefix = "Mixed") {
    return prefix + "P" + std::to_string(degree1) + "P" + std::to_string(degree2);
  }

  std::unique_ptr<IDoF> createDoF(LevelPair levels) const override {
    std::vector<std::unique_ptr<IDoF>> dofs(2);
    dofs[0] = std::make_unique<LagrangeDoF>(degree1, dim * size);
    dofs[1] = std::make_unique<LagrangeDoF>(degree2, size);
    return std::make_unique<MixedDoF>(std::move(dofs));
  }

  std::unique_ptr<ShapeT<T, sDim, tDim>> createShape(CELLTYPE type, int n, int) const override {
    if (n == 0)
      return std::unique_ptr<ShapeT<T, sDim, tDim>>(
          createLagrangeShape<T, sDim, tDim>(type, degree1));
    if (n == 1)
      return std::unique_ptr<ShapeT<T, sDim, tDim>>(
          createLagrangeShape<T, sDim, tDim>(type, degree2));
    THROW("Shape not implemented in MixedLagrangeDiscretization")
  }
public:
  MixedLagrangeDiscretizationT(const Meshes &meshes, int degree1, int degree2, int size = 1) :
      MixedLagrangeDiscretizationT(meshes, degree1, degree2, std::max(2 * degree1, 2 * degree2),
                                   size) {}

  MixedLagrangeDiscretizationT(const Meshes &meshes, int degree1, int degree2, int quadExactUpTo,
                               int size) :
      NonAdaptiveIDiscretizationT<T, sDim, tDim>(meshes, name(degree1, degree2), quadExactUpTo),
      degree1(degree1), degree2(degree2), dim(meshes.dim()), size(size) {}

  MixedLagrangeDiscretizationT(const NonAdaptiveIDiscretizationT<T, sDim, tDim> &disc, int dim,
                               int degree1, int degree2, int size = 1) :
      NonAdaptiveIDiscretizationT<T, sDim, tDim>(disc, name(degree1, degree2)), degree1(degree1),
      degree2(degree2), dim(dim), size(size) {}

  int Degree1() const { return degree1; }

  int Degree2() const { return degree2; }

  int Dim() const { return dim; }

  int Size() const { return size; }
};

using MixedLagrangeDiscretization = MixedLagrangeDiscretizationT<>;

#ifdef BUILD_IA

using IAMixedLagrangeDiscretization =
    MixedLagrangeDiscretizationT<IAInterval, SpaceDimension, TimeDimension>;

#endif

#endif // MIXEDLAGRANGEDISCRETIZATION_HPP
