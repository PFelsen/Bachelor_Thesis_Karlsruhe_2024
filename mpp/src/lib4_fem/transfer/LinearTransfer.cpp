#include "LinearTransfer.hpp"
#include <LagrangeNodalPoints.hpp>
#include <Vector.hpp>

void addRowListEntry(int fineId, const std::vector<int> &coarseIds,
                     const std::vector<double> &coarseWeights, RowList &rowList) {
  for (int i = 0; i < coarseIds.size(); ++i) {
    if (coarseWeights[i] > 0.0) { rowList.AddRowWeight(fineId, coarseIds[i], coarseWeights[i]); }
  }
}

int findInterpolatedRow(const std::vector<Point> &corners, const std::vector<double> &weights,
                        const Vector &fine) {
  Point finePoint = std::inner_product(
      corners.begin(), corners.end(), weights.begin(), Origin,
      [](Point sum, Point next) { return sum + next; },
      [](Point p, double weight) { return weight * p; });

  row fineRow = fine.find_row(finePoint);
  if (fineRow == fine.rows_end()) {
    return -1;
    THROW("Couldn't find interpolated row")
  }
  return fineRow.Id();
}

void linearTransferInterval(const vector<Point> &corners, const vector<int> &coarseIds,
                            const Vector &fine, double lvl, RowList &rowList) {
  for (int n = 0; n <= lvl; ++n) {
    std::vector<double> lagrangeWeights{((double)(lvl - n)) / lvl, n / lvl};
    int fineRowId = findInterpolatedRow(corners, lagrangeWeights, fine);
    if (fineRowId < 0) continue;
    if (rowList.CoarseRows(fineRowId).empty()) {
      addRowListEntry(fineRowId, coarseIds, lagrangeWeights, rowList);
    }
  }
}

void linearTransferTriangle(const vector<Point> &corners, const vector<int> &coarseIds,
                            const Vector &fine, double lvl, RowList &rowList) {
  for (int n = 0; n <= lvl; ++n) {
    for (int m = 0; m <= lvl - n; ++m) {
      std::vector<double> lagrangeWeights{(lvl - m - n) / lvl, m / lvl, n / lvl};

      int fineRowId = findInterpolatedRow(corners, lagrangeWeights, fine);

      if (fineRowId < 0) continue;
      if (rowList.CoarseRows(fineRowId).empty()) {
        addRowListEntry(fineRowId, coarseIds, lagrangeWeights, rowList);
      }
    }
  }
}

void linearTransferQuad(const vector<Point> &corners, const vector<int> &coarseIds,
                        const Vector &fine, double lvl, RowList &rowList) {
  for (int i = 0; i <= lvl; i++) {
    for (int j = 0; j <= lvl; j++) {
      double x = i / lvl;
      double y = j / lvl;

      std::vector<double> lagrangeWeights{
          (1 - x) * (1 - y),
          (1 - x) * y,
          x * (1 - y),
          x * y,
      };

      int fineRowId = findInterpolatedRow(corners, lagrangeWeights, fine);
      if (fineRowId < 0) continue;

      if (rowList.CoarseRows(fineRowId).empty()) {
        addRowListEntry(fineRowId, coarseIds, lagrangeWeights, rowList);
      }
    }
  }
}

void linearTransferTetrahedron(const vector<Point> &corners, const vector<int> &coarseIds,
                               const Vector &fine, double lvl, RowList &rowList) {
  for (int n = 0; n <= lvl; ++n) {
    for (int m = 0; m <= lvl - n; ++m) {
      for (int k = 0; k <= lvl - n - m; ++k) {
        std::vector<double> lagrangeWeights{(lvl - m - n - k) / lvl, m / lvl, n / lvl, k / lvl};

        int fineRowId = findInterpolatedRow(corners, lagrangeWeights, fine);
        if (fineRowId < 0) continue;

        if (rowList.CoarseRows(fineRowId).empty()) {
          addRowListEntry(fineRowId, coarseIds, lagrangeWeights, rowList);
        }
      }
    }
  }
}

void linearTransferHexahedron(const vector<Point> &corners, const vector<int> &coarseIds,
                              const Vector &fine, double lvl, RowList &rowList) {

  for (int i = 0; i <= lvl; i++) {
    for (int j = 0; j <= lvl; j++) {
      for (int k = 0; k <= lvl; ++k) {
        double a = i / ((double)lvl);
        double b = j / ((double)lvl);
        double c = k / ((double)lvl);

        std::vector<double> lagrangeWeights{
            (1 - a) * (1 - b) * (1 - c),
            (1 - a) * (1 - b) * c,
            (1 - a) * b * (1 - c),
            (1 - a) * b * c,
            a * (1 - b) * (1 - c),
            a * (1 - b) * c,
            a * b * (1 - c),
            a * b * c,
        };

        int fineRowId = findInterpolatedRow(corners, lagrangeWeights, fine);
        if (fineRowId < 0) continue;

        if (rowList.CoarseRows(fineRowId).empty()) {
          addRowListEntry(fineRowId, coarseIds, lagrangeWeights, rowList);
        }
      }
    }
  }
}

void mpp::linearTransfer(CELLTYPE celltype, const vector<Point> &coarseNodalPoints,
                         const vector<int> &coarseIds, double lvl, const Vector &fine,
                         RowList &rowList) {
  switch (celltype) {
  case INTERVAL:
    linearTransferInterval(coarseNodalPoints, coarseIds, fine, lvl, rowList);
    break;
  case TRIANGLE:
    linearTransferTriangle(coarseNodalPoints, coarseIds, fine, lvl, rowList);
    break;
  case QUADRILATERAL:
    linearTransferQuad(coarseNodalPoints, coarseIds, fine, lvl, rowList);
    break;
  case TETRAHEDRON:
    linearTransferTetrahedron(coarseNodalPoints, coarseIds, fine, lvl, rowList);
    break;
  case HEXAHEDRON:
    linearTransferHexahedron(coarseNodalPoints, coarseIds, fine, lvl, rowList);
    break;
  default:
    THROW("Celltype not implemented")
  }
}

vector<int> mpp::nodalIds(const Vector &vec, const vector<Point> &nodalPoints) {
  std::vector<int> coarseIds(nodalPoints.size());
  for (int i = 0; i < nodalPoints.size(); ++i) {
    coarseIds[i] = vec.find_row(nodalPoints[i]).Id();
  }
  return coarseIds;
}

RowList mpp::constructLinearRowList(const Vector &coarse, const Vector &fine, int coarseDegree,
                                    int fineDegree) {
  RowList rowList(coarse, fine);

  int dlevel = (int)log2(coarseDegree);
  if (dlevel < log2(coarseDegree)) {
    THROW("degree " + std::to_string(coarseDegree) + " transfer not implemented")
  }
  if (((int)log2(coarseDegree)) + fine.SpaceLevel() < dlevel + coarse.SpaceLevel()) {
    THROW("Can't interpolate: Coarser mesh is finer than fine mesh.")
  }

  double lvl = pow(2, fine.SpaceLevel() - (coarse.SpaceLevel() + dlevel)) * fineDegree;

  std::vector<Point> childrenVerteces{};
  std::vector<Point> childVerteces{};
  std::vector<Point> coarseNodalPoints{};
  for (cell c = coarse.cells(); c != coarse.cells_end(); ++c) {
    std::vector<std::vector<std::unique_ptr<Cell>>> recursiveCells(dlevel + 1);
    recursiveCells[0].emplace_back(CreateCell(c.Type(), c.Subdomain(), c.AsVector()));
    for (int l = 1; l <= dlevel; ++l) {
      for (const auto &parent : recursiveCells[l - 1]) {
        const vector<Rule> &R = parent->Refine(childrenVerteces);
        for (int i = 0; i < R.size(); ++i) {
          R[i](childrenVerteces, childVerteces);
          recursiveCells[l].emplace_back(CreateCell(R[i].type(), c.Subdomain(), childVerteces));
        }
      }
    }
    for (const auto &child : recursiveCells[dlevel]) {
      LagrangeNodalPoints(*child, coarseNodalPoints, 1);
      vector<int> coarseIds = mpp::nodalIds(coarse, coarseNodalPoints);
      mpp::linearTransfer(child->ReferenceType(), coarseNodalPoints, coarseIds, lvl, fine, rowList);
    }
  }
  return rowList;
}