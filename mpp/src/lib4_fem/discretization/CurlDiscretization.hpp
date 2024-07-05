#ifndef CURLDISCRETIZATION_HPP
#define CURLDISCRETIZATION_HPP

#include "CurlDoF.hpp"
#include "CurlShapes.hpp"
#include "IDiscretization.hpp"
#include "Meshes.hpp"

template<typename T = double, int sDim = SpaceDimension, int tDim = TimeDimension>
class CurlDiscretizationT : public NonAdaptiveIDiscretizationT<T, sDim, tDim> {
  int degree{};
  int size{};
protected:
  std::unique_ptr<IDoF> createDoF(LevelPair levels) const override {
    return std::make_unique<CurlDoF>(degree, size);
  }

  std::unique_ptr<ShapeT<T, sDim, tDim>> createShape(CELLTYPE type, int n, int) const override {
    if (n == 0)
      return std::unique_ptr<ShapeT<T, sDim, tDim>>(createCurlShape<T, sDim, tDim>(type, degree));
    THROW("Shape not implemented in CurlDiscretization")
  }
public:
  CurlDiscretizationT(const Meshes &meshes, int degree, int size = 1) :
      CurlDiscretizationT(meshes, degree, 2 * (degree + 1), size) {}

  CurlDiscretizationT(const Meshes &meshes, int degree, int quadExactUpTo, int size) :
      NonAdaptiveIDiscretizationT<T, sDim, tDim>(meshes, "P" + std::to_string(degree) + "Curl",
                                                 quadExactUpTo),
      degree(degree), size(size) {}

  CurlDiscretizationT(const NonAdaptiveIDiscretizationT<T, sDim, tDim> &disc, int degree,
                      int size = 1) :
      NonAdaptiveIDiscretizationT<T, sDim, tDim>(disc, "P" + std::to_string(degree) + "Curl"),
      degree(degree), size(size) {}
};

using CurlDiscretization = CurlDiscretizationT<>;

#endif // CURLDISCRETIZATION_HPP
