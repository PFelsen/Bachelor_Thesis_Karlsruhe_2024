#ifndef MULTIPARTDISCRETIZATION_HPP
#define MULTIPARTDISCRETIZATION_HPP

#include <functional>
#include <utility>

#include <Vector.hpp>
#include "IDiscretization.hpp"
#include "LagrangeDiscretization.hpp"
#include "LagrangeShapes.hpp"
#include "Meshes.hpp"
#include "MultiPartDoF.hpp"

template<typename T = double, int sDim = SpaceDimension, int tDim = TimeDimension>
class MultiPartDiscretizationT : public NonAdaptiveIDiscretizationT<T, sDim, tDim> {
  int degree{};
  int size{};

  std::function<int(const Cell &)> cellIndex;
protected:
  std::unique_ptr<IDoF> createDoF(LevelPair levels) const override {
    auto multipartDoF = std::make_unique<MultiPartDoF>(degree, this->meshes.fine(), size);

    auto interfaceDisc =
        std::make_shared<const LagrangeDiscretization>((*this).GetMeshes(), degree, 2);

    Vector interface(0.0, interfaceDisc);
    interface.SetAccumulateFlag(false);
    const auto &M = (*this).GetMeshes().fine();

    vector<Point> z;
    for (cell c = M.cells(); c != M.cells_end(); ++c) {
      int sd = cellIndex(*c);
      z = multipartDoF->GetNodalPoints(*c);

      for (auto &nodalPoint : z)
        interface(nodalPoint, sd) = 1;
    }
    // TODO: Is accumulate needed here?
    interface.Accumulate();

    for (row r = interface.rows(); r != interface.rows_end(); ++r) {
      if (interface(r, 0) > 0 && interface(r, 1) > 0) multipartDoF->AddInterfacePoint(r());
    }
    multipartDoF->FinishInterfacePoints();

    return multipartDoF;
  }

  std::unique_ptr<ShapeT<T, sDim, tDim>> createShape(CELLTYPE type, int n, int) const override {
    if (n == 0)
      return std::unique_ptr<ShapeT<T, sDim, tDim>>(
          createLagrangeShape<T, sDim, tDim>(type, degree));
    THROW("Shape not implemented in MultiPartDiscretization")
  }
public:
  MultiPartDiscretizationT(const Meshes &meshes, int degree,
                           std::function<int(const Cell &)> mpIndex, int size = 1) :
      MultiPartDiscretizationT(meshes, degree, 2 * degree, std::move(mpIndex), size) {}

#pragma clang diagnostic push
#pragma ide diagnostic ignored "performance-unnecessary-value-param"

  MultiPartDiscretizationT(const Meshes &meshes, int degree, int quadExactUpTo,
                           std::function<int(const Cell &)> mpIndex, int size) :
      NonAdaptiveIDiscretizationT<T, sDim, tDim>(meshes, "MultiPartP" + std::to_string(degree),
                                                 quadExactUpTo),
      degree(degree), size(size), cellIndex(std::move(mpIndex)) {}

#pragma clang diagnostic pop

  MultiPartDiscretizationT(const MultiPartDiscretizationT<T, sDim, tDim> &disc, int size = 1) :
      NonAdaptiveIDiscretizationT<T, sDim, tDim>(disc,
                                                 "MultiPartP" + std::to_string(disc.Degree())),
      degree(disc.Degree()), size(size), cellIndex(disc.cellIndex) {}

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
};

using MultiPartDiscretization = MultiPartDiscretizationT<>;

#endif // MULTIPARTDISCRETIZATION_HPP
