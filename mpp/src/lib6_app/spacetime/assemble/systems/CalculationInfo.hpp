#ifndef __CALCULATION_INFO__
#define __CALCULATION_INFO__

#include "IterableSystem.hpp"

template<class Element>
struct ElementHolder {
  const Element &element;
  DGRowEntries &rowEntries;

  ElementHolder(const Element &element, DGRowEntries &rowEntries)
      : element(element), rowEntries(rowEntries) {}
};

template<class Element>
struct ElementWeightCacheable : public ElementHolder<Element> {
  const size_t quadratureSize;
  std::vector<double> weightCache;

  ElementWeightCacheable(const Element &element,
                         std::function<double(Point)> &weightFunction,
                         DGRowEntries &rowEntries)
      : ElementHolder<Element>(element, rowEntries),
        quadratureSize(element.GetQuadratureSize()),
        weightCache(quadratureSize) {
    for (size_t quadratureIndex = 0; quadratureIndex < quadratureSize;
         ++quadratureIndex) {
      const Point &quadraturePoint = element.qPoint[quadratureIndex];
      const double weight =
          element.qWeight[quadratureIndex] * weightFunction(quadraturePoint);
      weightCache[quadratureIndex] = weight;
    }
  }
};

struct InternalRowEntry {
  Scalar *internalAccess;
  size_t accessModifier;
  size_t multiplier;

  InternalRowEntry(Matrix &matrix, const Point &cellPoint) : multiplier(0) {
    const auto row = matrix.find_row(cellPoint);
    internalAccess = matrix(row, row);
    accessModifier = row.n();
  }

  inline Scalar *begin() const {
    return internalAccess + multiplier * accessModifier;
  }

  inline operator size_t() const { return multiplier; }
};

inline InternalRowEntry &operator++(InternalRowEntry &internalEntry) {
  internalEntry.multiplier++;
  return internalEntry;
}

inline InternalRowEntry &operator+=(InternalRowEntry &internalEntry,
                                    size_t amount) {
  internalEntry.multiplier += amount;
  return internalEntry;
}

struct CalculationInfo {
  const Point &point;
  Matrix &matrix;
  DGRowEntries rowEntries;
  SpaceTimeViscoAcousticDGTElement element;
  std::shared_ptr<AcousticProblem> problem;
  const double rho;
  Vector &vector;
  cell cellIterator;
  std::vector<double> kappaInverse;
  std::vector<double> kappaTauInverse;
  std::shared_ptr<STDiscretization> discretization;
  std::function<double(Point)> &weightFunction;
  ElementWeightCacheable<SpaceTimeViscoAcousticDGTElement> cache;
  const size_t quadratureSize;
  std::vector<double> &weightCache;

  CalculationInfo(Matrix &matrix, std::shared_ptr<AcousticProblem> problem,
                  Vector &vector, cell cellIterator,
                  std::function<double(Point)> &weightFunction,
                  std::shared_ptr<STDiscretization> discretization)
      : point(cellIterator()),
        matrix(matrix),
        rowEntries(matrix, *cellIterator, *cellIterator, false),
        element(matrix, cellIterator, problem->nL()),
        problem(problem),
        rho(problem->Rho(*cellIterator, point)),
        vector(vector),
        cellIterator(cellIterator),
        kappaInverse(2 + problem->nL()),
        kappaTauInverse(2 + problem->nL()),
        discretization(discretization),
        weightFunction(weightFunction),
        cache(element, weightFunction, rowEntries),
        quadratureSize(cache.quadratureSize),
        weightCache(cache.weightCache) {

    kappaInverse[0] = 0.0;
    kappaTauInverse[0] = kappaTauInverse[1] = 0.0;
    kappaInverse[1] = 1.0 / problem->Kappa_i(*cellIterator, point, 0);
    for (int j = 2; j < 2 + problem->nL(); j++) {
      kappaInverse[j] = 1.0 / problem->Kappa_i(*cellIterator, point, j - 1);
      kappaTauInverse[j] = kappaInverse[j] / problem->Tau_i(*cellIterator, point, j - 1);
    }
  }
};

struct SecondElementInfo {
  SpaceTimeViscoAcousticDGTFaceElement element;
  DGRowEntries &rowEntries;

  SecondElementInfo(const STDiscretization &disc, const VectorMatrixBase &M,
                    const cell &c, int f_id, int nL, DGRowEntries *rows)
      : element(disc, M, c, f_id, nL), rowEntries(*rows) {}
};

struct CalculateFaceInfo {
  bool isOnBoundary;
  Point facePoint;
  int boundaryID;
  const cell faceCellIterator;
  Point pointForFaceCell;
  int neighbourFaceID;
  SpaceTimeViscoAcousticDGTFaceElement element;
  DGRowEntries rowEntries;
  SecondElementInfo secondElementInfo;
  size_t quadratureSize;
  std::vector<double> weightCache;
  std::vector<double> fluxCache;
  double alpha[4];
  double z_K;

  CalculateFaceInfo(CalculationInfo &info, const size_t faceID,
                    const Matrix &matrix, const double z_K)
      : isOnBoundary(matrix.OnBoundary(*info.cellIterator, faceID)),
        facePoint(info.cellIterator.Face(faceID)),
        boundaryID(isOnBoundary ? info.problem->BndID(facePoint) : -1),
        faceCellIterator(isOnBoundary ? info.cellIterator
                                      : matrix.find_neighbour_cell(info.cellIterator, faceID)),
        pointForFaceCell(faceCellIterator()),
        neighbourFaceID(
            matrix.find_neighbour_face_id(facePoint, faceCellIterator)),
        element(*info.discretization, info.matrix, info.cellIterator, faceID,
                info.problem->nL()),
        rowEntries(info.matrix, *info.cellIterator, *faceCellIterator, false),
        secondElementInfo(*info.discretization, info.matrix, faceCellIterator,
                          neighbourFaceID, info.problem->nL(), &rowEntries),
        quadratureSize(element.nQ()),
        z_K(z_K) {
    weightCache.resize(info.quadratureSize);
    const auto fullDimension = info.element.GetFullDimensions();
    fluxCache.resize(fullDimension * quadratureSize);
    for (const auto &[index, weight, point]:
        quadratureIterableWithPoint(*this)) {
      const VectorField &quadratureNormal = element.qNormal[index];
      weightCache[index] = element.qWeight[index] * info.weightFunction(point);
      for (const auto [i, velocity]:
          velocityIterable(index, *this, size_t{0})) {
        fluxCache[i + fullDimension * index] = quadratureNormal * velocity;
      }
    }

    const double z_Kf = sqrt(info.problem->Rho(*faceCellIterator, pointForFaceCell) *
                             info.problem->Kappa(*faceCellIterator, pointForFaceCell));
    const auto alpha1 = 1.0 / (z_K + z_Kf);
    alpha[0] = alpha1;
    alpha[1] = z_Kf * z_K * alpha1;
    alpha[2] = z_Kf * alpha1;
    alpha[3] = z_K * alpha1;
  }
};

#endif