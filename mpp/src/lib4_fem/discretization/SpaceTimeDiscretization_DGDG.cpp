#include "SpaceTimeDiscretization_DGDG.hpp"

#include "AdaptiveSpaceTimeDof_DGDG.hpp"
#include "LagrangeShapes.hpp"
#include "SpaceTimeMatrixGraph.hpp"

template<typename T, int sDim, int tDim>
void STDiscretizationT_DGDG<T, sDim, tDim>::init() {
  CELLTYPE space_ct = SpaceCellType(this->cell_type);
  using ShapeMemory = ShapeQuadsOnCellT<T, sDim, tDim>;
  using Quad = const QuadratureT<T, sDim, tDim> &;

  Quad SQ = GetQuadratureT<T, sDim, tDim>(FaceCellType(space_ct), 2 * this->defaultDegree.space);
  auto spaceShapeQuads = std::make_unique<ShapeMemory>(space_ct, SQ);

  Quad TQ = GetQuadratureT<T, sDim, tDim>(POINT, 2 * this->defaultDegree.time);
  auto timeShapeQuads = std::make_unique<ShapeMemory>(INTERVAL, TQ);

  for (int degree = 0; degree <= 6; ++degree) {
    auto shape = [degree](CELLTYPE c) {
      return createLagrangeShapeUnique<T, sDim, tDim>(c, degree);
    };
    Quad SQdeg = GetQuadratureT<T, sDim, tDim>(space_ct, 2 * degree);
    spaceShapeQuads->InsertShapeCellQuad(shape(space_ct), SQdeg);
    Quad TQdeg = GetQuadratureT<T, sDim, tDim>(INTERVAL, 2 * degree);
    timeShapeQuads->InsertShapeCellQuad(shape(INTERVAL), TQdeg);
  }
  this->shapeQuadsMap[{space_ct, this->cellIdx()}] = std::move(spaceShapeQuads);
  this->shapeQuadsMap[{INTERVAL, this->timeIdx()}] = std::move(timeShapeQuads);
  this->timeFaceQuad = GetQuadratureT<T, sDim, tDim>(INTERVAL, 2 * this->defaultDegree.time);
}

template<typename T, int sDim, int tDim>
STDiscretizationT_DGDG<T, sDim, tDim>::STDiscretizationT_DGDG(const Meshes &M, DegreePair degree,
                                                              int pDim, bool singleEntries,
                                                              bool blockDiagonal) :
    STDiscretizationT<T, sDim, tDim>(M, degree, pDim, singleEntries, blockDiagonal) {
  init();
}

template<typename T, int sDim, int tDim>
STDiscretizationT_DGDG<T, sDim, tDim>::STDiscretizationT_DGDG(
    const STDiscretizationT_DGDG &otherDisc, bool singleEntries, bool blockDiagonal) :
    STDiscretizationT<T, sDim, tDim>(otherDisc, singleEntries, blockDiagonal) {
  init();
}

template<typename T, int sDim, int tDim>
void STDiscretizationT_DGDG<T, sDim, tDim>::adaptQuadrature(LevelPair levels, int quadDeg,
                                                            bool adaptCellQuad) {
  using ShapeMemory = ShapeQuadsOnCellT<T, sDim, tDim>;

  DegreePair maxdeg = (*this)(levels).GetDoF().GetMaxDegree();
  int maxSpace = PPM->Max(max(quadDeg, int(maxdeg.space)), levels.commSplit);
  int maxTime = PPM->Max(max(quadDeg, int(maxdeg.time)), levels.commSplit);

  CELLTYPE space_ct = SpaceCellType(this->cell_type);

  if (adaptCellQuad) {
    std::unique_ptr<ShapeMemory> spaceShapeQuads =
        std::make_unique<ShapeMemory>(space_ct,
                                      GetQuadratureT<T, sDim, tDim>(FaceCellType(space_ct),
                                                                    2 * maxSpace));

    std::unique_ptr<ShapeMemory> timeShapeQuads =
        std::make_unique<ShapeMemory>(INTERVAL, GetQuadratureT<T, sDim, tDim>(POINT, 2 * maxTime));

    for (int degree = 0; degree <= 6; ++degree) {
      auto shape = [degree](CELLTYPE c) { return createLagrangeShape<T, sDim, tDim>(c, degree); };
      spaceShapeQuads->InsertShapeCellQuad(std::unique_ptr<ShapeT<T, sDim, tDim>>(shape(space_ct)),
                                           GetQuadratureT<T, sDim, tDim>(space_ct, 2 * maxSpace));
      timeShapeQuads->InsertShapeCellQuad(std::unique_ptr<ShapeT<T, sDim, tDim>>(shape(INTERVAL)),
                                          GetQuadratureT<T, sDim, tDim>(INTERVAL, 2 * maxTime));
    }
    this->shapeQuadsMap[std::make_pair(space_ct, this->cellIdx())] = std::move(spaceShapeQuads);
    this->shapeQuadsMap[std::make_pair(INTERVAL, this->timeIdx())] = std::move(timeShapeQuads);
  }

  this->shapeQuadsMap.at(std::pair<CELLTYPE, int>(space_ct, this->cellIdx()))
      ->UpdateFaceQuad(GetQuadratureT<T, sDim, tDim>(FaceCellType(space_ct), 2 * maxSpace));
  this->timeFaceQuad = GetQuadratureT<T, sDim, tDim>(INTERVAL, 2 * maxTime);
  if (adaptCellQuad) {
    this->st_quad.clear();
    this->st_shapes.clear();
    for (short sdeg = 0; sdeg < 6; sdeg++) {
      for (short tdeg = 0; tdeg < 6; tdeg++) {
        DegreePair key{sdeg, tdeg};
        SpaceTimeShapeT<T, sDim, tDim> stshape(this->GetCellShape(sdeg), GetTimeShape(tdeg, false));
        auto st_quad = GetSpaceTimeQuad({short(maxSpace), short(maxTime)});
        stshape.fillValues(st_quad.QPoints());
        this->st_shapes.insert({key, stshape});
      }
    }
  }
}

template<typename T, int sDim, int tDim>
std::unique_ptr<IDoF> STDiscretizationT_DGDG<T, sDim, tDim>::createDoF(LevelPair levels) const {
  EquidistantNodalPointProvider npp;
  return std::make_unique<AdaptiveSpaceTimeDof_DGDG>(this->meshes[levels], this->pDim, npp);
}
