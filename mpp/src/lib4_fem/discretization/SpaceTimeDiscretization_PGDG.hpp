#ifndef SPACETIMEDISCRETIZATION_PGDG_HPP
#define SPACETIMEDISCRETIZATION_PGDG_HPP

#include "LagrangeShapes.hpp"
#include "SpaceTimeDiscretization.hpp"

template<typename T = double, int sDim = SpaceDimension, int tDim = TimeDimension>
class STDiscretizationT_PGDG : public STDiscretizationT<T, sDim, tDim> {
  void init();
protected:
  std::unique_ptr<IDoF> createDoF(LevelPair levels) const override;
public:
  STDiscretizationT_PGDG(const Meshes &meshes, DegreePair degree, int pDim,
                         bool singleEntries = false, bool blockDiagonal = false);

  STDiscretizationT_PGDG(const STDiscretizationT_PGDG &otherDisc, bool singleEntries = false,
                         bool blockDiagonal = false);

  std::shared_ptr<STDiscretizationT<>> create(const Meshes &meshes, DegreePair degree,
                                              int pDim) const override {
    return std::make_shared<STDiscretizationT_PGDG>(meshes, degree, pDim);
  };

  void adaptQuadrature(LevelPair levels, int quadDeg, bool adaptCellQuad = false) override;

  const ShapeT<T, sDim, tDim> &GetTimeShape(int deg, bool testSpace) const override {
    if (testSpace) { return this->GetShape(INTERVAL, this->testIdx(), deg); }
    return this->GetShape(INTERVAL, this->timeIdx(), deg);
  }

  const QuadratureT<T, sDim, tDim> &GetTimeQuad(int deg, bool testSpace) const override {
    if (testSpace) { return this->GetQuad(INTERVAL, this->testIdx(), deg); }
    return this->GetQuad(INTERVAL, this->timeIdx(), deg);
  }

  bool isDgInTime() const { return false; }

  const SpaceTimeQuadratureT<T, sDim, tDim> &GetSpaceTimeQuad(DegreePair degs) const override {
    THROW("Implement")
    auto quadIter = this->st_quad.find(degs);
    if (quadIter == this->st_quad.end()) {
      SpaceTimeQuadratureT<T, sDim, tDim> stquad(this->GetCellQuad(degs.space),
                                                 GetTimeQuad(degs.time, false));
      quadIter = this->st_quad.insert({degs, stquad}).first;
    }
    return quadIter->second;
  }

  const SpaceTimeShapeT<T, sDim, tDim> &GetSpaceTimeShape(DegreePair degs) const override {
    THROW("Implement")
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
};

typedef STDiscretizationT_PGDG<> STDiscretization_PGDG;

template class STDiscretizationT_PGDG<>;

#endif // SPACETIMEDISCRETIZATION_HPP
