#ifndef SPACETIMEDISCRETIZATION_HPP
#define SPACETIMEDISCRETIZATION_HPP

#include "AdaptiveSpaceTimeDof.hpp"
#include "Config.hpp"
#include "Functools.hpp"
#include "IDiscretization.hpp"
#include "Mesh.hpp"
#include "Shapes.hpp"
#include "SpaceTimeMatrixGraph.hpp"
#include "SpaceTimeQuadrature.hpp"
#include "SpaceTimeShape.hpp"

#include <memory>


template<class Value>
using SpaceTimeMap = std::unordered_map<DegreePair, Value>;

template<typename T = double, int sDim = SpaceDimension, int tDim = TimeDimension>
class STDiscretizationT : public IDiscretizationT<T, sDim, tDim> {
protected:
  const CELLTYPE cell_type;
  QuadratureT<T, sDim, tDim> timeFaceQuad;

  bool singleEntries = false;
  bool blockDiagonal = false;
  DegreePair defaultDegree{0, 0};
  int pDim = 0;

  const STDiscretizationT *otherDisc = nullptr; // as reference

  mutable SpaceTimeMap<SpaceTimeShapeT<T, sDim, tDim>> st_shapes;
  mutable SpaceTimeMap<SpaceTimeQuadratureT<T, sDim, tDim>> st_quad;

  const ShapeQuadsOnCellT<T, sDim, tDim> &ShapeQuads(CELLTYPE type, int n = 0) const override {
    if (otherDisc) return otherDisc->ShapeQuads(type, n);
    return *this->shapeQuadsMap.at(std::make_pair(type, n));
  }

  int maxDegree(int n) const override { THROW("Not needed in space time") }

  int quadExactUpTo(int n, int degree) const override { THROW("Not needed in space time") }

  int faceQuadExactUpTo(int n) const override{THROW("Not needed in space time")}

  std::unique_ptr<ShapeT<T, sDim, tDim>> createShape(CELLTYPE type, int n,
                                                     int degree) const override {
    THROW("Not needed in space time")
  }

  using IDiscretization::createDoF;

  std::unique_ptr<IMatrixGraph> createMatrixGraph(LevelPair levels) const override {
    return std::make_unique<SpaceTimeMatrixGraph>(this->meshes[levels],
                                                  std::move(this->createDoF(levels)), defaultDegree,
                                                  isDgInTime(), singleEntries, levels,
                                                  std::unordered_map<Point, DegreePair>{},
                                                  blockDiagonal);
  }

  std::unique_ptr<IMatrixGraph>
  createMatrixGraph(LevelPair levels, const std::unordered_map<Point, DegreePair> &map) const {
    return std::make_unique<SpaceTimeMatrixGraph>(this->meshes[levels],
                                                  std::move(this->createDoF(levels)), defaultDegree,
                                                  isDgInTime(), singleEntries, levels, map,
                                                  blockDiagonal);
  }

  STDiscretizationT(const Meshes &meshes, DegreePair degree, int pDim, bool singleEntries = false,
                    bool blockDiagonal = false);

  STDiscretizationT(const STDiscretizationT &otherDisc, bool singleEntries = false,
                    bool blockDiagonal = false);
public:
  using IDiscretization ::GetShape;

  const ShapeT<T, sDim, tDim> &GetShape(const Cell &c, int n = 0, int degree = 0) const override {
    Exit("Not save.") return GetSpaceTimeShape(c);
  }

  using IDiscretization::GetQuad;

  const QuadratureT<T, sDim, tDim> &GetQuad(const Cell &c, int n = 0,
                                            int degree = 0) const override {
    Exit("Not save.") return GetSpaceTimeQuad(c);
  }

  virtual const SpaceTimeShapeT<T, sDim, tDim> &GetSpaceTimeShape(DegreePair degs) const = 0;

  const SpaceTimeShapeT<T, sDim, tDim> &GetSpaceTimeShape(const Cell &c) const {
    Exit("Not safe.") DegreePair deg = (*this)().GetDoF().GetDegree(c());
    return GetSpaceTimeShape(deg);
  }

  virtual const SpaceTimeQuadratureT<T, sDim, tDim> &GetSpaceTimeQuad(DegreePair degs) const = 0;

  virtual const SpaceTimeQuadratureT<T, sDim, tDim> &GetSpaceTimeQuad(const Cell &c) const {
    Exit("Not safe.") DegreePair deg = (*this)().GetDoF().GetDegree(c);
    return GetSpaceTimeQuad(deg);
  }

  virtual void adaptQuadrature(LevelPair levels, int quadDeg, bool adaptCellQuad) = 0;

  virtual const ShapeT<T, sDim, tDim> &GetTimeShape(int deg, bool testSpace = false) const = 0;

  virtual const QuadratureT<T, sDim, tDim> &GetTimeQuad(int deg, bool testSpace = false) const = 0;

  constexpr int cellIdx() const { return 0; }

  constexpr int timeIdx() const { return 1; }

  constexpr int testIdx() const { return 2; }

  const ShapeT<T, sDim, tDim> &GetCellShape(const Cell &c, int deg) const {
    return this->GetShape(c, cellIdx(), deg);
  }

  const ShapeT<T, sDim, tDim> &GetCellShape(int deg) const {
    return this->GetShape(SpaceCellType(cell_type), cellIdx(), deg);
  }

  const QuadratureT<T, sDim, tDim> &GetCellQuad(const Cell &c, int deg = 0) const {
    return this->GetQuad(c, cellIdx(), deg);
  }

  // TODO: remove face
  const QuadratureT<T, sDim, tDim> &GetCellFaceQuad(const Cell &c, int face, int deg = 0) const {
    return this->GetFaceQuad(c, cellIdx());
  }

  const QuadratureT<T, sDim, tDim> &GetCellQuad(int deg = 0) const {
    return this->GetQuad(SpaceCellType(cell_type), cellIdx(), deg);
  }

  const vector<vector<PointT<T, sDim, tDim>>> &FaceQPoints() const {
    return IDiscretization::FaceQPoints(SpaceCellType(cell_type), cellIdx());
  }

  const QuadratureT<T, sDim, tDim> &GetCellFaceQuad(int face, int deg = 0) const {
    return this->GetFaceQuad(SpaceCellType(cell_type), cellIdx());
  }

  const QuadratureT<T, sDim, tDim> &GetTimeFaceQuad(int face, int deg = 0) const {
    return timeFaceQuad;
  }

  virtual bool isDgInTime() const = 0;

  virtual std::shared_ptr<STDiscretizationT<>> create(const Meshes &meshes, DegreePair degree,
                                                      int pDim) const = 0;

  void communicate(LevelPair level) const {
    (*this)(level).GetDoF().communicate();

    (*this)(level).update();
  }

  void CreateAdaptivityLevel(const std::unordered_map<Point, DegreePair> &map, LevelPair pair);
};

typedef STDiscretizationT<> STDiscretization;

template class STDiscretizationT<>;


#endif // SPACETIMEDISCRETIZATION_HPP
