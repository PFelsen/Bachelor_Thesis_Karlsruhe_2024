#ifndef IDISCRETIZATION_HPP
#define IDISCRETIZATION_HPP

#include "Config.hpp"
#include "MatrixGraph.hpp"
#include "Meshes.hpp"
#include "Quadrature.hpp"
#include "Shape.hpp"
#include "Synchronization.hpp"

class Vector;

/**
 * Contains space level (first) and time level (second)
 * std::pair is used since it provides comparison operators for maps automatically
 */

/**
 * Contains shape and corresponding quadrature for a fixed degree/order
 */
template<typename T, int sDim, int tDim>
class ShapeCellQuadT {
  const QuadratureT<T, sDim, tDim> &cellQuad;
  std::unique_ptr<ShapeT<T, sDim, tDim>> shape = nullptr;
public:
  ShapeCellQuadT(std::unique_ptr<ShapeT<T, sDim, tDim>> &&shape,
                 const QuadratureT<T, sDim, tDim> &cellQuad,
                 const vector<vector<PointT<T, sDim, tDim>>> &localFaceQuadPoints) :
      shape(std::move(shape)), cellQuad(cellQuad) {
    this->shape->fillValues(cellQuad.QPoints());
    UpdateFaceValues(localFaceQuadPoints);
  }

  void UpdateFaceValues(const vector<vector<PointT<T, sDim, tDim>>> &localFaceQuadPoints) {
    shape->fillFaceValues(localFaceQuadPoints);
  }

  const QuadratureT<T, sDim, tDim> &CellQuad() const { return cellQuad; }

  const ShapeT<T, sDim, tDim> &GetShape() const { return *shape; }

  template<typename T0, int sDim0, int tDim0>
  friend class ShapeQuadsOnCellT;
};

/**
 * Contains ShapeQuad's of different degrees/orders, the face quadrature and local face quadrature
 * points Note that
 */
template<typename T, int sDim, int tDim>
class ShapeQuadsOnCellT {
  CELLTYPE type;
  /// face quadrature should be the same for all degrees/orders
  const QuadratureT<T, sDim, tDim> *faceQuad = nullptr; // as reference
  /// local face quadrature points corresponding to faceQuad (first index = face)
  mutable vector<vector<PointT<T, sDim, tDim>>> localFaceQuadPoints{};
  /// int corresponds to different degrees/orders of the shapes
  mutable std::vector<std::unique_ptr<ShapeCellQuadT<T, sDim, tDim>>> shapeCellQuads{};
public:
  ShapeQuadsOnCellT(ShapeQuadsOnCellT &&sq) :
      type(sq.type), faceQuad(sq.faceQuad), localFaceQuadPoints(std::move(sq.localFaceQuadPoints)),
      shapeCellQuads(std::move(sq.shapeCellQuads)) {}

  ShapeQuadsOnCellT(CELLTYPE type, const QuadratureT<T, sDim, tDim> &faceQuad) : type(type) {
    UpdateFaceQuad(faceQuad);
  }

  void UpdateFaceQuad(const QuadratureT<T, sDim, tDim> &faceQuad) {
    this->faceQuad = &faceQuad;
    localFaceQuadPoints.clear();
    const Cell *refCell = ReferenceCell(type);
    localFaceQuadPoints.resize(refCell->Faces());
    for (int face = 0; face < refCell->Faces(); ++face) {
      localFaceQuadPoints[face].resize(faceQuad.size());
      for (int q = 0; q < faceQuad.size(); ++q) {
        PointT<T, sDim, tDim> QP = faceQuad.QPoint(q);
        localFaceQuadPoints[face][q] = refCell->FaceLocalToGlobal(face, QP);
      }
    }
    for (int i = 0; i < shapeCellQuads.size(); ++i)
      shapeCellQuads[i]->UpdateFaceValues(localFaceQuadPoints);
  }

  /// Note that degree also can be the order (cf. RTElement)
  void InsertShapeCellQuad(std::unique_ptr<ShapeT<T, sDim, tDim>> &&shape,
                           const QuadratureT<T, sDim, tDim> &cellQuad) {
    shape->fillValues(cellQuad.QPoints());
    shape->fillFaceValues(localFaceQuadPoints);
    shapeCellQuads.push_back(std::make_unique<ShapeCellQuadT<T, sDim, tDim>>(std::move(shape),
                                                                             cellQuad,
                                                                             localFaceQuadPoints));
  }

  const QuadratureT<T, sDim, tDim> &FaceQuad() const { return *faceQuad; }

  const QuadratureT<T, sDim, tDim> &CellQuad(int degree = 0) const {
    return shapeCellQuads.at(degree)->CellQuad();
  }

  const ShapeT<T, sDim, tDim> &GetShape(int degree = 0) const {
    return shapeCellQuads.at(degree)->GetShape();
  }

  const vector<vector<PointT<T, sDim, tDim>>> &

  LocalFaceQPoints() const {
    return localFaceQuadPoints;
  }

  template<typename T0, int sDim0, int tDim0>
  friend class IDiscretizationT;
};

void transformGenericCellData(DataSynchronization &synchronization, const std::vector<int> &indices,
                              const Vector &data);

void transformGenericPointData(DataSynchronization &synchronization,
                               const std::vector<int> &indices, const Vector &data);

template<typename T = double, int sDim = SpaceDimension, int tDim = TimeDimension>
class IDiscretizationT {
protected:
  int verbose = 1;
  string discName;

  const Meshes &meshes;

  /**
   * Shapes and Quadratures
   * Note that the integer index numbers the different shape types (and not the degree)
   * Different degrees are contained in ShapeQuadsOnCellT
   */
  mutable std::map<std::pair<CELLTYPE, int>, std::unique_ptr<ShapeQuadsOnCellT<T, sDim, tDim>>>
      shapeQuadsMap{};

  virtual int quadExactUpTo(int n, int degree) const = 0;

  virtual const QuadratureT<T, sDim, tDim> &getQuad(CELLTYPE type, int n, int degree) const {
    return GetQuadratureT<T, sDim, tDim>(type, quadExactUpTo(n, degree));
  };

  virtual int faceQuadExactUpTo(int n) const = 0;

  virtual const QuadratureT<T, sDim, tDim> &getFaceQuad(CELLTYPE type, int n) const {
    return GetQuadratureT<T, sDim, tDim>(FaceCellType(type), faceQuadExactUpTo(n));
  }

  virtual std::unique_ptr<ShapeT<T, sDim, tDim>> createShape(CELLTYPE type, int n,
                                                             int degree) const = 0;

  virtual int maxDegree(int n) const = 0;

  virtual const ShapeQuadsOnCellT<T, sDim, tDim> &ShapeQuads(CELLTYPE type, int n = 0) const {
    if (const auto element = shapeQuadsMap.find(std::make_pair(type, n));
        element != shapeQuadsMap.end()) {
      return *element->second;
    }

    std::pair<CELLTYPE, int> key(type, n);
    auto shapeQuads =
        std::make_unique<ShapeQuadsOnCellT<T, sDim, tDim>>(type, getFaceQuad(type, n));

    for (int degree = 0; degree <= maxDegree(n); ++degree)
      shapeQuads->InsertShapeCellQuad(std::move(createShape(type, n, degree)),
                                      getQuad(type, n, degree));
    shapeQuadsMap[key] = std::move(shapeQuads);

    return *shapeQuadsMap.at(key);
  }

  /**
   * MatrixGraph
   * Note that the keys are now the same levels as for the corresponding mesh
   */
  mutable std::unordered_map<LevelPair, std::unique_ptr<IMatrixGraph>> graphs{};

  // TODO: use Levels in the following functions!
  virtual std::unique_ptr<IDoF> createDoF(LevelPair levels) const = 0;

  // override this if other MatrixGraph is needed
  virtual std::unique_ptr<IMatrixGraph> createMatrixGraph(LevelPair levels) const {
    return std::make_unique<MatrixGraph>(meshes[levels], std::move(createDoF(levels)));
  }

  explicit IDiscretizationT(const Meshes &meshes, const string name) :
      meshes(meshes), discName(name) {}
public:
  virtual ~IDiscretizationT() = default;

  void PrintInfo() const {
    mout.PrintInfo("Discretization", verbose, PrintInfoEntry("Discretization Name", discName));
  }

  string DiscName() const { return discName; }

  const Meshes &GetMeshes() const { return meshes; }

  virtual const ShapeT<T, sDim, tDim> &GetShape(const Cell &c, int n = 0, int degree = 0) const {
    return GetShape(c.ReferenceType(), n, degree);
  }

  const ShapeT<T, sDim, tDim> &GetShape(CELLTYPE type, int n = 0, int degree = 0) const {
    return ShapeQuads(type, n).GetShape(degree);
  }

  virtual const QuadratureT<T, sDim, tDim> &GetQuad(const Cell &c, int n = 0,
                                                    int degree = 0) const {
    return GetQuad(c.ReferenceType(), n, degree);
  }

  const QuadratureT<T, sDim, tDim> &GetQuad(CELLTYPE type, int n = 0, int degree = 0) const {
    return ShapeQuads(type, n).CellQuad(degree);
  }

  const QuadratureT<T, sDim, tDim> &GetFaceQuad(const Cell &c, int n = 0) const {
    return GetFaceQuad(c.ReferenceType(), n);
  }

  const QuadratureT<T, sDim, tDim> &GetFaceQuad(CELLTYPE type, int n = 0) const {
    return ShapeQuads(type, n).FaceQuad();
  }

  const vector<vector<PointT<T, sDim, tDim>>> &FaceQPoints(const Cell &c, int n = 0) const {
    return FaceQPoints(c.ReferenceType(), n);
  }

  const vector<vector<PointT<T, sDim, tDim>>> &FaceQPoints(CELLTYPE type, int n = 0) const {
    return ShapeQuads(type, n).LocalFaceQPoints();
  }

  IMatrixGraph &operator()(LevelPair levels = {-1, -1}) const {
    levels = meshes.AdaptLevels(levels);

    if (const auto element = graphs.find(levels); element != graphs.end()) {
      return *element->second;
    }

    const int coarseLevel = meshes.Settings().coarseLevel;
    const int fineLevel = meshes.Settings().fineLevel;
    if (coarseLevel <= levels.space && levels.space <= fineLevel) {
      graphs[levels] = std::move(createMatrixGraph(levels));
      return *graphs.at(levels);
    }

    Exit("No MatrixGraph for level at spaceLevel: " + std::to_string(levels.space)
          + " and timeLevel: " + std::to_string(levels.time))
  }

  /**
   * @brief Transforms data to plot into a container @see DataSynchronization
   * accepts
   *
   * @param synchronization the synchronization object which accepts the
   * transformed point data via @see DataSynchronization#AddPointData and/or the
   * cell data via @see DataSynchronization#AddCellData
   * @param indices the indices of @see Vector where the data to plot can be
   * found
   * @param data the data
   */
  virtual void transformData(DataSynchronization &synchronization, const std::vector<int> &indices,
                             const Vector &data) const {
    Warning("Plot type for Discretization not defined.");
  }

  /**
   * @brief Transforms deformations of a mesh into a common format @see
   * MeshSynchronization accepts
   *
   * @param synchronization the synchronization object which accepts the
   * transformed deformation
   * @param displacement the displacement
   */
  virtual void transformDisplacement(MeshSynchronization &synchronization,
                                     const Vector &displacement) const {
    Warning("Deformation type for Discretization not defined.");
  }
};

template<typename T = double, int sDim = SpaceDimension, int tDim = TimeDimension>
std::ostream &operator<<(std::ostream &s, const IDiscretizationT<T, sDim, tDim> &D) {
  return s << D.DiscName();
}

template<typename T = double, int sDim = SpaceDimension, int tDim = TimeDimension>
class NonAdaptiveIDiscretizationT : public IDiscretizationT<T, sDim, tDim> {
protected:
  int quadsExact{-1};

  int quadExactUpTo(int n, int degree) const override { return quadsExact; }

  int faceQuadExactUpTo(int n) const override { return quadsExact; }

  int maxDegree(int n) const override { return 0; }

  NonAdaptiveIDiscretizationT(const Meshes &meshes, std::string name, int quadsExactUpTo) :
      IDiscretizationT<T, sDim, tDim>(meshes, name), quadsExact(quadsExactUpTo) {}

  NonAdaptiveIDiscretizationT(const NonAdaptiveIDiscretizationT<T, sDim, tDim> &disc,
                              std::string name) :
      IDiscretizationT<T, sDim, tDim>(disc.GetMeshes(), name), quadsExact(disc.quadsExact) {}
};

using IDiscretization = IDiscretizationT<>;
using NonAdaptiveIDiscretization = NonAdaptiveIDiscretizationT<>;

#ifdef BUILD_IA

using IAIDiscretization = IDiscretizationT<IAInterval, SpaceDimension, TimeDimension>;
using IANonAdaptiveIDiscretization =
    NonAdaptiveIDiscretizationT<IAInterval, SpaceDimension, TimeDimension>;

#endif


#endif // IDISCRETIZATION_HPP
