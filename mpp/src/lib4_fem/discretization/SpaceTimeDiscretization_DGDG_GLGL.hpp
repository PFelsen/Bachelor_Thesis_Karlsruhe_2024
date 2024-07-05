#ifndef SPACETIMEDISCRETIZATION_DGDG_GLGL_HPP
#define SPACETIMEDISCRETIZATION_DGDG_GLGL_HPP

#include "SpaceTimeDiscretization.hpp"

template<typename T = double, int sDim = SpaceDimension, int tDim = TimeDimension>
class STDiscretizationT_DGDG_GLGL : public STDiscretizationT<T, sDim, tDim> {
  void init();
protected:
  std::unique_ptr<IDoF> createDoF(LevelPair levels) const override;
public:
  STDiscretizationT_DGDG_GLGL(const Meshes &meshes, DegreePair degree, int pDim,
                              bool singleEntries = false, bool blockDiagonal = false);

  STDiscretizationT_DGDG_GLGL(const STDiscretizationT_DGDG_GLGL &otherDisc,
                              bool singleEntries = false, bool blockDiagonal = false);

  std::shared_ptr<STDiscretizationT<>> create(const Meshes &meshes, DegreePair degree,
                                              int pDim) const override {
    return std::make_shared<STDiscretizationT_DGDG_GLGL>(meshes, degree, pDim);
  };

  void adaptQuadrature(LevelPair levels, int quadDeg, bool adaptCellQuad) override;

  const ShapeT<T, sDim, tDim> &GetTimeShape(int deg, bool testSpace = false) const override {
    return this->GetShape(INTERVAL, this->timeIdx(), deg);
  }

  const QuadratureT<T, sDim, tDim> &GetTimeQuad(int deg, bool testSpace = false) const override {
    return this->GetQuad(INTERVAL, this->timeIdx(), deg);
  }

  bool isDgInTime() const override { return true; }

  const SpaceTimeQuadratureT<T, sDim, tDim> &GetSpaceTimeQuad(DegreePair degs) const override {
    auto quadIter = this->st_quad.find(degs);
    if (quadIter == this->st_quad.end()) {
      SpaceTimeQuadratureT<T, sDim, tDim> stquad(this->GetCellQuad(degs.space),
                                                 GetTimeQuad(degs.time));
      quadIter = this->st_quad.insert({degs, stquad}).first;
    }
    return quadIter->second;
  }

  const SpaceTimeShapeT<T, sDim, tDim> &GetSpaceTimeShape(DegreePair degs) const override {
    auto shapeIter = this->st_shapes.find(degs);
    if (shapeIter == this->st_shapes.end()) {
      SpaceTimeShapeT<T, sDim, tDim> stshape(this->GetCellShape(degs.space),
                                             GetTimeShape(degs.time));
      auto st_quad = GetSpaceTimeQuad(degs);
      stshape.fillValues(st_quad.QPoints());
      shapeIter = this->st_shapes.insert({degs, stshape}).first;
    }
    return shapeIter->second;
  }
};

typedef STDiscretizationT_DGDG_GLGL<> STDiscretization_DGDG_GLGL;

template class STDiscretizationT_DGDG_GLGL<>;

#endif // SPACETIMEDISCRETIZATION_HPP