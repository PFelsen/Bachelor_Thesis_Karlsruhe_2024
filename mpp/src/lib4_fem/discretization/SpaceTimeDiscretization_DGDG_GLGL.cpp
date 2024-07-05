#include "SpaceTimeDiscretization_DGDG_GLGL.hpp"
#include "AdaptiveSpaceTimeDof_DGDG.hpp"
#include "SpaceTimeMatrixGraph.hpp"

#include "GaussLobattoShapes.hpp"

template<typename T, int sDim, int tDim>
void STDiscretizationT_DGDG_GLGL<T, sDim, tDim>::init() {
  this->timeFaceQuad = GetGLQuadratureT<T, sDim, tDim>(INTERVAL, 2 * this->defaultDegree.time);

  std::unique_ptr<ShapeQuadsOnCellT<T, sDim, tDim>> spaceShapeQuads = std::make_unique<
      ShapeQuadsOnCellT<T, sDim, tDim>>(QUADRILATERAL,
                                        GetGLQuadratureT<T, sDim, tDim>(INTERVAL,
                                                                        this->defaultDegree.space
                                                                            + 1));
  std::unique_ptr<ShapeQuadsOnCellT<T, sDim, tDim>> timeShapeQuads =
      std::make_unique<ShapeQuadsOnCellT<T, sDim, tDim>>(INTERVAL, this->timeFaceQuad);

  for (int degree = 0; degree <= 6; ++degree) {
    auto shape = [degree](CELLTYPE c) { return createGaussLobattoShape<T, sDim, tDim>(c, degree); };
    spaceShapeQuads->InsertShapeCellQuad(std::unique_ptr<ShapeT<T, sDim, tDim>>(
                                             shape(QUADRILATERAL)),
                                         GetGLQuadratureT<T, sDim, tDim>(QUADRILATERAL,
                                                                         degree + 1));
    timeShapeQuads->InsertShapeCellQuad(std::unique_ptr<ShapeT<T, sDim, tDim>>(shape(INTERVAL)),
                                        GetGLQuadratureT<T, sDim, tDim>(INTERVAL, degree + 1));
  }
  this->shapeQuadsMap[std::make_pair(QUADRILATERAL, this->cellIdx())] = std::move(spaceShapeQuads);
  this->shapeQuadsMap[std::make_pair(INTERVAL, this->timeIdx())] = std::move(timeShapeQuads);
}

template<typename T, int sDim, int tDim>
STDiscretizationT_DGDG_GLGL<T, sDim, tDim>::STDiscretizationT_DGDG_GLGL(const Meshes &M,
                                                                        DegreePair degree, int pDim,
                                                                        bool singleEntries,
                                                                        bool blockDiagonal) :
    STDiscretizationT<T, sDim, tDim>(M, degree, pDim, singleEntries, blockDiagonal) {
  init();
}

template<typename T, int sDim, int tDim>
STDiscretizationT_DGDG_GLGL<T, sDim, tDim>::STDiscretizationT_DGDG_GLGL(
    const STDiscretizationT_DGDG_GLGL &otherDisc, bool singleEntries, bool blockDiagonal) :
    STDiscretizationT<T, sDim, tDim>(otherDisc, singleEntries, blockDiagonal) {
  init();
}

template<typename T, int sDim, int tDim>
void STDiscretizationT_DGDG_GLGL<T, sDim, tDim>::adaptQuadrature(LevelPair levels, int quadDeg,
                                                                 bool adaptCellQuad) {
  DegreePair maxdeg = (*this)(levels).GetDoF().GetMaxDegree();
  int maxSpace = PPM->Max(max(quadDeg, int(maxdeg.space)), levels.commSplit);
  int maxTime = PPM->Max(max(quadDeg, int(maxdeg.time)), levels.commSplit);

  this->shapeQuadsMap.at(std::pair<CELLTYPE, int>(QUADRILATERAL, this->cellIdx()))
      ->UpdateFaceQuad(GetQuadratureT<T, sDim, tDim>(QUADRILATERAL, 2 * maxSpace));
  this->timeFaceQuad = GetGLQuadratureT<T, sDim, tDim>(INTERVAL, 2 * maxTime);
}

template<typename T, int sDim, int tDim>
std::unique_ptr<IDoF>
STDiscretizationT_DGDG_GLGL<T, sDim, tDim>::createDoF(LevelPair levels) const {
  GaussLobattoNodalPointProvider npp;
  return std::make_unique<AdaptiveSpaceTimeDof_DGDG>(this->meshes[levels], this->pDim, npp);
}
