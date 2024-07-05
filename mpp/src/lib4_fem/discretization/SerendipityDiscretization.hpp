#ifndef SERENDIPITYDISCRETIZATIONT_HPP
#define SERENDIPITYDISCRETIZATIONT_HPP

#include "IDiscretization.hpp"
#include "LagrangeDoF.hpp"
#include "Meshes.hpp"
#include "SerendipityDoF.hpp"
#include "SerendipityShapes.hpp"

template<typename T = double, int sDim = SpaceDimension, int tDim = TimeDimension>
class SerendipityDiscretizationT : public NonAdaptiveIDiscretizationT<T, sDim, tDim> {
  int degree{0};
  int size{};
protected:
  std::unique_ptr<IDoF> createDoF(LevelPair levels) const override {
    if (degree == 2) return std::make_unique<SerendipityDoF>(size);
    return std::make_unique<LagrangeDoF>(degree, size);
  }

  std::unique_ptr<ShapeT<T, sDim, tDim>> createShape(CELLTYPE type, int n, int) const override {
    if (n == 0) {
      if (degree == 2)
        return std::unique_ptr<ShapeT<T, sDim, tDim>>(createSerendipityShape<T, sDim, tDim>(type));
      return std::unique_ptr<ShapeT<T, sDim, tDim>>(
          createLagrangeShape<T, sDim, tDim>(type, degree));
    }

    THROW("Shape not implemented in LagrangeDiscretization")
  }
public:
  explicit SerendipityDiscretizationT(const Meshes &meshes, int degree = 2, int size = 1) :
      SerendipityDiscretizationT(meshes, degree, 2 * degree, size) {}

  SerendipityDiscretizationT(const Meshes &meshes, int degree, int quadExactUpTo, int size) :
      NonAdaptiveIDiscretizationT<T, sDim, tDim>(meshes,
                                                 (degree == 2) ? "SerendipityP2" :
                                                                 ("P" + std::to_string(degree)),
                                                 quadExactUpTo),
      degree(degree), size(size) {}

  SerendipityDiscretizationT(const NonAdaptiveIDiscretizationT<T, sDim, tDim> &disc, int degree = 2,
                             int size = 1) :
      NonAdaptiveIDiscretizationT<T, sDim, tDim>(disc, (degree == 2) ?
                                                           "SerendipityP2" :
                                                           ("P" + std::to_string(degree))),
      degree(degree), size(size) {}
};

using SerendipityDiscretization = SerendipityDiscretizationT<>;

#ifdef BUILD_IA

using IASerendipityDiscretization =
    SerendipityDiscretizationT<IAInterval, SpaceDimension, TimeDimension>;

#endif

#endif // SERENDIPITYDISCRETIZATIONT_HPP
