#include "SpaceTimeDiscretization_PGDG.hpp"
#include "AdaptiveSpaceTimeDof_PGDG.hpp"
#include "SpaceTimeMatrixGraph.hpp"

template<typename T, int sDim, int tDim>
void STDiscretizationT_PGDG<T, sDim, tDim>::init() {
  std::unique_ptr<ShapeQuadsOnCellT<T, sDim, tDim>> spaceShapeQuads =
      std::make_unique<ShapeQuadsOnCellT<
          T, sDim, tDim>>(QUADRILATERAL,
                          GetQuadratureT<T, sDim, tDim>(INTERVAL, 2 * this->defaultDegree.space));

  std::unique_ptr<ShapeQuadsOnCellT<T, sDim, tDim>> timeShapeQuads =
      std::make_unique<ShapeQuadsOnCellT<
          T, sDim, tDim>>(INTERVAL,
                          GetQuadratureT<T, sDim, tDim>(POINT, 2 * this->defaultDegree.time));

  std::unique_ptr<ShapeQuadsOnCellT<T, sDim, tDim>> testShapeQuads =
      std::make_unique<ShapeQuadsOnCellT<
          T, sDim, tDim>>(INTERVAL,
                          GetQuadratureT<T, sDim, tDim>(POINT, 2 * this->defaultDegree.time));

  for (int degree = 0; degree <= 6; ++degree) {
    auto shape = [degree](CELLTYPE c) { return createLagrangeShape<T, sDim, tDim>(c, degree); };
    spaceShapeQuads->InsertShapeCellQuad(std::unique_ptr<ShapeT<T, sDim, tDim>>(
                                             shape(QUADRILATERAL)),
                                         GetQuadratureT<T, sDim, tDim>(QUADRILATERAL, 2 * degree));
    timeShapeQuads->InsertShapeCellQuad(std::unique_ptr<ShapeT<T, sDim, tDim>>(shape(INTERVAL)),
                                        GetQuadratureT<T, sDim, tDim>(INTERVAL,
                                                                      max(0, 2 * (degree - 1))));
    testShapeQuads->InsertShapeCellQuad(std::unique_ptr<ShapeT<T, sDim, tDim>>(shape(INTERVAL)),
                                        GetQuadratureT<T, sDim, tDim>(INTERVAL, 2 * degree));
  }
  this->shapeQuadsMap[std::make_pair(QUADRILATERAL, this->cellIdx())] = std::move(spaceShapeQuads);
  this->shapeQuadsMap[std::make_pair(INTERVAL, this->timeIdx())] = std::move(timeShapeQuads);
  this->shapeQuadsMap[std::make_pair(INTERVAL, this->testIdx())] = std::move(testShapeQuads);
  this->timeFaceQuad = GetQuadratureT<T, sDim, tDim>(INTERVAL, 2 * this->defaultDegree.time);
}

template<typename T, int sDim, int tDim>
STDiscretizationT_PGDG<T, sDim, tDim>::STDiscretizationT_PGDG(const Meshes &M, DegreePair degree,
                                                              int pDim, bool singleEntries,
                                                              bool blockDiagonal) :
    STDiscretizationT<T, sDim, tDim>(M, degree, pDim, singleEntries, blockDiagonal) {
  init();
}

template<typename T, int sDim, int tDim>
STDiscretizationT_PGDG<T, sDim, tDim>::STDiscretizationT_PGDG(
    const STDiscretizationT_PGDG &otherDisc, bool singleEntries, bool blockDiagonal) :
    STDiscretizationT<T, sDim, tDim>(otherDisc, singleEntries, blockDiagonal) {
  init();
}

template<typename T, int sDim, int tDim>
void STDiscretizationT_PGDG<T, sDim, tDim>::adaptQuadrature(LevelPair levels, int quadDeg,
                                                            bool adaptCellQuad) {
  DegreePair maxdeg = (*this)(levels).GetDoF().GetMaxDegree();
  int maxSpace = PPM->Max(max(quadDeg, int(maxdeg.space)), levels.commSplit);
  int maxTime = PPM->Max(max(quadDeg, int(maxdeg.time)), levels.commSplit) - 1;

  this->shapeQuadsMap.at(std::pair<CELLTYPE, int>(QUADRILATERAL, this->cellIdx()))
      ->UpdateFaceQuad(GetQuadratureT<T, sDim, tDim>(QUADRILATERAL, 2 * maxSpace));
  this->timeFaceQuad = GetQuadratureT<T, sDim, tDim>(INTERVAL, 2 * maxTime);
}

// Todo Think about time Degree for PGDG!!!!
template<typename T, int sDim, int tDim>
std::unique_ptr<IDoF> STDiscretizationT_PGDG<T, sDim, tDim>::createDoF(LevelPair levels) const {
  return std::make_unique<AdaptiveSpaceTimeDof_PGDG>(this->meshes[levels], this->pDim,
                                                     new EquidistantNodalPointProvider());
}