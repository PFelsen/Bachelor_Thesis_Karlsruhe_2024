#ifndef MIXEDSERENDIPITYDISCRETIZATION_HPP
#define MIXEDSERENDIPITYDISCRETIZATION_HPP

#include "IDiscretization.hpp"
#include "LagrangeDoF.hpp"
#include "Meshes.hpp"
#include "MixedDoF.hpp"
#include "SerendipityDoF.hpp"
#include "SerendipityShapes.hpp"

template<typename T = double, int sDim = SpaceDimension, int tDim = TimeDimension>
class MixedSerendipityDiscretizationT : public NonAdaptiveIDiscretizationT<T, sDim, tDim> {
  int degree1{};
  int degree2{};
  int dim{};
  int size{};
protected:
  std::string name(int degree1, int degree2, std::string prefix = "Mixed") {
    if (degree1 == 2) {
      if (degree2 == 2)
        return prefix + "SerendipityP" + std::to_string(degree1) + "SerendipityP"
               + std::to_string(degree2);
      else return prefix + "SerendipityP" + std::to_string(degree1) + "P" + std::to_string(degree2);
    } else {
      if (degree2 == 2)
        return prefix + "P" + std::to_string(degree1) + "SerendipityP" + std::to_string(degree2);
      else return prefix + "P" + std::to_string(degree1) + "P" + std::to_string(degree2);
    }
  }

  std::unique_ptr<IDoF> createDoF(LevelPair levels) const override {
    std::vector<std::unique_ptr<IDoF>> dofs(2);
    if (degree1 == 2) dofs[0] = std::make_unique<SerendipityDoF>(dim * size);
    else dofs[0] = std::make_unique<LagrangeDoF>(degree1, dim * size);
    if (degree2 == 2) dofs[1] = std::make_unique<SerendipityDoF>(size);
    else dofs[1] = std::make_unique<LagrangeDoF>(degree2, size);
    return std::make_unique<MixedDoF>(std::move(dofs));
  }

  std::unique_ptr<ShapeT<T, sDim, tDim>> createShape(CELLTYPE type, int n, int) const override {
    if (n == 0) {
      if (degree1 == 2)
        return std::unique_ptr<ShapeT<T, sDim, tDim>>(createSerendipityShape<T, sDim, tDim>(type));
      return std::unique_ptr<ShapeT<T, sDim, tDim>>(
          createLagrangeShape<T, sDim, tDim>(type, degree1));
    }
    if (n == 1) {
      if (degree2 == 2)
        return std::unique_ptr<ShapeT<T, sDim, tDim>>(createSerendipityShape<T, sDim, tDim>(type));
      return std::unique_ptr<ShapeT<T, sDim, tDim>>(
          createLagrangeShape<T, sDim, tDim>(type, degree2));
    }
    THROW("Shape not implemented in MixedSerendipityDiscretization")
  }
public:
  MixedSerendipityDiscretizationT(const Meshes &meshes, int degree1, int degree2, int size = 1) :
      MixedSerendipityDiscretizationT(meshes, degree1, degree2, std::max(2 * degree1, 2 * degree2),
                                      size) {}

  MixedSerendipityDiscretizationT(const Meshes &meshes, int degree1, int degree2, int quadExactUpTo,
                                  int size) :
      NonAdaptiveIDiscretizationT<T, sDim, tDim>(meshes, name(degree1, degree2), quadExactUpTo),
      degree1(degree1), degree2(degree2), dim(meshes.dim()), size(size) {}

  MixedSerendipityDiscretizationT(const NonAdaptiveIDiscretizationT<T, sDim, tDim> &disc, int dim,
                                  int degree1, int degree2, int size = 1) :
      NonAdaptiveIDiscretizationT<T, sDim, tDim>(disc, name(degree1, degree2)), degree1(degree1),
      degree2(degree2), dim(dim), size(size) {}
};

using MixedSerendipityDiscretization = MixedSerendipityDiscretizationT<>;

#ifdef BUILD_IA

using IAMixedSerendipityDiscretization =
    MixedSerendipityDiscretizationT<IAInterval, SpaceDimension, TimeDimension>;

#endif

#endif // MIXEDSERENDIPITYDISCRETIZATION_HPP
