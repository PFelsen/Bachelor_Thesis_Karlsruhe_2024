#ifndef SPACETIMEDISCRETIZATION_DGDG_HPP
#define SPACETIMEDISCRETIZATION_DGDG_HPP

#include "SpaceTimeDiscretization.hpp"
#include "Vector.hpp"
#include "VtuPlot.hpp"

template<typename T = double, int sDim = SpaceDimension, int tDim = TimeDimension>
class STDiscretizationT_DGDG : public STDiscretizationT<T, sDim, tDim> {
  void init();

  using STDiscretizationT<T, sDim, tDim>::GetCellShape;
protected:
  std::unique_ptr<IDoF> createDoF(LevelPair levels) const override;
public:
  STDiscretizationT_DGDG(const Meshes &meshes, DegreePair degree, int pDim,
                         bool singleEntries = false, bool blockDiagonal = false);

  STDiscretizationT_DGDG(const STDiscretizationT_DGDG &otherDisc, bool singleEntries = false,
                         bool blockDiagonal = false);

  std::shared_ptr<STDiscretizationT<>> create(const Meshes &meshes, DegreePair degree,
                                              int pDim) const override {
    return std::make_shared<STDiscretizationT_DGDG>(meshes, degree, pDim);
  };

  void adaptQuadrature(LevelPair levels, int quadDeg, bool adaptCellQuad = false) override;

  const ShapeT<T, sDim, tDim> &GetTimeShape(int deg, bool testSpace) const override {
    return this->GetShape(INTERVAL, this->timeIdx(), deg);
  }

  const QuadratureT<T, sDim, tDim> &GetTimeQuad(int deg, bool testSpace) const override {
    return this->GetQuad(INTERVAL, this->timeIdx(), deg);
  }

  const SpaceTimeQuadratureT<T, sDim, tDim> &GetSpaceTimeQuad(DegreePair degs) const override {
    auto quadIter = this->st_quad.find(degs);
    if (quadIter == this->st_quad.end()) {
      SpaceTimeQuadratureT<T, sDim, tDim> stquad(this->GetCellQuad(degs.space),
                                                 GetTimeQuad(degs.time, false));
      quadIter = this->st_quad.insert({degs, stquad}).first;
    }
    return quadIter->second;
  }

  const SpaceTimeShapeT<T, sDim, tDim> &GetSpaceTimeShape(DegreePair degs) const override {
    auto shapeIter = this->st_shapes.find(degs);
    if (shapeIter == this->st_shapes.end()) {
      SpaceTimeShapeT<T, sDim, tDim> stshape(this->GetCellShape(degs.space),
                                             GetTimeShape(degs.time, false));
      auto st_quad = GetSpaceTimeQuad(degs);
      stshape.fillValues(st_quad.QPoints());
      shapeIter = this->st_shapes.insert({degs, stshape}).first;
    }
    return shapeIter->second;
  }

  bool isDgInTime() const override { return true; }

  void transformData(DataSynchronization &synchronization, const std::vector<int> &indices,
                     const Vector &data) const override {
    const auto backend = &synchronization.GetBackend();
    if (dynamic_cast<const DefaultSpaceTimeBackend *>(backend) != nullptr) {
      transformGenericCellData(synchronization, indices, data);
    } else if (dynamic_cast<const SpaceTimeBackend *>(backend) != nullptr) {
      transformGenericCellData(synchronization, indices, data);
    } else {
      THROW("Wrong backend used with space time discretization!");
    }
  }
};

typedef STDiscretizationT_DGDG<> STDiscretization_DGDG;

template class STDiscretizationT_DGDG<>;

#endif // SPACETIMEDISCRETIZATION_HPP
