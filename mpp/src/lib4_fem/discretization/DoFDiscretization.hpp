#ifndef DOFDISCRETIZATION_HPP
#define DOFDISCRETIZATION_HPP

#include <concepts>

#include "BasicDoFs.hpp"
#include "DGDoF.hpp"
#include "IDiscretization.hpp"
#include "Meshes.hpp"

/**
 * DoFDiscretization can be initialized with an arbitrary DoF specified via the additional
 * template parameter DOF
 */
template<typename DOF, typename T = double, int sDim = SpaceDimension, int tDim = TimeDimension>
class DoFDiscretizationT : public IDiscretizationT<T, sDim, tDim> {
protected:
  int size{};

  std::unique_ptr<IDoF> createDoF(LevelPair levels) const override {
    return std::make_unique<DOF>(size);
  }

  std::unique_ptr<ShapeT<T, sDim, tDim>> createShape(CELLTYPE type, int n, int) const override {
    THROW("Not implemented in DoFDiscretization")
  }

  int quadExactUpTo(int n, int degree) const override {
    THROW("Not implemented in DoFDiscretization")
  }

  int faceQuadExactUpTo(int n) const override { THROW("Not implemented in DoFDiscretization") }

  int maxDegree(int n) const override { THROW("Not implemented in DoFDiscretization") return 0; }
public:
  DoFDiscretizationT(const Meshes &meshes, int size = 1) :
      IDiscretizationT<T, sDim, tDim>(meshes, "DoFDiscretization"), size(size) {}

  void transformData(DataSynchronization &synchronization, const std::vector<int> &indices,
                     const Vector &data) const override {
    if constexpr (std::same_as<DOF, CellDoF> || std::same_as<DOF, DGDoF>) {
      transformGenericCellData(synchronization, indices, data);
    } else {
      static_assert(std::same_as<DOF, VertexDoF>);
      transformGenericPointData(synchronization, indices, data);
    }
  }
};

template<typename DOF>
using DoFDiscretization = DoFDiscretizationT<DOF>;

#ifdef BUILD_IA

template<typename DOF>
using IADoFDiscretization = DoFDiscretizationT<DOF, IAInterval, SpaceDimension, TimeDimension>;

#endif

#endif // DOFDISCRETIZATION_HPP
