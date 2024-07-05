#ifndef LAGRANGENODALPOINTS_HPP
#define LAGRANGENODALPOINTS_HPP

#include "Cell.hpp"

/**
 * Numbering of nodal points on interval (example P3)
 *
 *   0 -- 1 -- 2 -- 3
 */

template<typename T, int sDim, int tDim>
void LagrangeNodalPointsInterval(const Cell &c, vector<PointT<T, sDim, tDim>> &z, int degree) {
  if (degree == 0) {
    z[0] = PointT<T, sDim, tDim>(c());
  } else {
    int index = 0;
    for (int i = 0; i < (degree + 1); i++) {
      z[index] = c[0] + ((c[1] - PointT<T, sDim, tDim>(c[0])) * (T(i) / degree));
      index++;
    }
  }
}

/**
 * Numbering of nodal points on triangle (example P3)
 *
 *   2
 *   |  -
 *   5    7
 *   |       -
 *   8    9    4
 *   |            -
 *   0 -- 3 -- 6 -- 1
 */

template<typename T, int sDim, int tDim>
void LagrangeNodalPointsTriangle(const Cell &c, vector<PointT<T, sDim, tDim>> &z, int index,
                                 int degree, int currentDegree, int numberOfRekCalls) {
  if (currentDegree == 0) { z[index++] = c(); }
  for (int n = 0; n < currentDegree; n++) {
    z[index++] = c[0] + ((c[1] - PointT<T, sDim, tDim>(c[0])) * (T(n) / degree))
                 + (c[1] + PointT<T, sDim, tDim>(c[2]) - 2 * PointT<T, sDim, tDim>(c[0]))
                       * (T(numberOfRekCalls) / degree);
    z[index++] = c[1] + ((c[2] - PointT<T, sDim, tDim>(c[1])) * (T(n) / degree))
                 + (c[2] + PointT<T, sDim, tDim>(c[0]) - 2 * PointT<T, sDim, tDim>(c[1]))
                       * (T(numberOfRekCalls) / degree);
    z[index++] = c[2] + ((c[0] - PointT<T, sDim, tDim>(c[2])) * (T(n) / degree))
                 + (c[1] + PointT<T, sDim, tDim>(c[0]) - 2 * PointT<T, sDim, tDim>(c[2]))
                       * (T(numberOfRekCalls) / degree);
  }
  if (index != ((degree + 1) * (degree + 2)) / 2) {
    LagrangeNodalPointsTriangle(c, z, index, degree, currentDegree - 3, ++numberOfRekCalls);
  }
}

template<typename T, int sDim, int tDim>
void LagrangeNodalPointsTriangle(const Cell &c, vector<PointT<T, sDim, tDim>> &z, int degree) {
  LagrangeNodalPointsTriangle(c, z, 0, degree, degree, 0);
}

template<typename T, int sDim, int tDim>
void LagrangeNodalPointsTrianglePrism(const Cell &c, vector<PointT<T, sDim, tDim>> &z, int degree) {
  int numberOfSpatialNodalPoints = z.size() / (degree + 1);
  if (z.size() % (degree + 1) != 0) { THROW("z.size must be divisible by deg+1") }
  vector<PointT<T, sDim, tDim>> spatialNodalPoints(numberOfSpatialNodalPoints);
  LagrangeNodalPointsTriangle(c, spatialNodalPoints, 0, degree, degree, 0);
  for (int d = 0; d < degree + 1; d++) {
    for (int s = 0; s < numberOfSpatialNodalPoints; s++) {
      z[d * numberOfSpatialNodalPoints + s] = spatialNodalPoints[s].WithT(((double)d) / degree);
    }
  }
}

/**
 * Numbering of nodal points on quadrilateral (example P3)
 *
 *  12 - 13 - 14 - 15
 *  |              |
 *  8    9    10   11
 *  |              |
 *  4    5    6    7
 *  |              |
 *  0 -- 1 -- 2 -- 3
 */

template<typename T, int sDim, int tDim>
void LagrangeNodalPointsQuadrilateral(const Cell &c, vector<PointT<T, sDim, tDim>> &z, int degree) {
  if (degree == 0) {
    z[0] = PointT<T, sDim, tDim>(c());
  } else {
    int index = 0;
    for (int i = 0; i < (degree + 1); i++) {
      for (int j = 0; j < (degree + 1); j++) {
        PointT<T, sDim, tDim> v = c[0] + j / T(degree) * (c[1] - c[0]);
        PointT<T, sDim, tDim> w = c[3] + j / T(degree) * (c[2] - c[3]);
        z[index] = v + i / T(degree) * (w - v);
        index++;
      }
    }
  }
}

/**
 * Numbering of nodal points on tetrahedron (example P3)
 *
 *  front:                                           rear:
 *
 *  19
 *  |  -
 *  16   17               18
 *  |       -             |  -
 *  10   11   12          13   14         15
 *  | /          -        | /     -       | /-
 *  0 -- 1 -- 2 -- 3      4 -- 5 -- 6     7 -- 8     9
 *                       /           \   /      \
 */

template<typename T, int sDim, int tDim>
void LagrangeNodalPointsTetrahedron(const Cell &c, vector<PointT<T, sDim, tDim>> &z, int degree) {
  int index = 0;
  if (degree == 0) {
    z[index++] = c();
    return;
  }
  for (int k = 0; k < degree + 1; ++k) {
    for (int n = 0; n < degree + 1 - k; ++n) {
      for (int m = 0; m < degree + 1 - k - n; ++m) {
        z[index++] = c[0] + ((c[1] - PointT<T, sDim, tDim>(c[0])) * (T(m) / degree))
                     + ((c[2] - PointT<T, sDim, tDim>(c[0])) * (T(n) / degree))
                     + ((c[3] - PointT<T, sDim, tDim>(c[0])) * (T(k) / degree));
      }
    }
  }
}

// Used to identify GetNodalPointsOnFace * Numbering of nodal points on tetrahedron (example P4)
// *
// *  front:                                           rear:
// *
// *  34
// *  |  -
// *  31  32                    33
// *  |       -                 |  -
// *  25   26   27              28   29              30
// *  |            -            |       -            |  -
// *  15   16   17   18         19   20   21         22   23          24
// *  | /               -       | /          -       | /      -       | / -
// *  0 -- 1 -- 2 -- 3 -- 4     5 -- 6 -- 7 -- 8     9 -- 10 -- 11    12 -- 13    14
// *                           /                \   /             \
// */
static std::vector<std::vector<std::vector<int>>>
    tetrahdronNodalMap{{{0, 2, 1}, {1, 2, 3}, {0, 3, 2}, {0, 1, 3}},
                       {{0, 3, 5, 1, 4, 2},
                        {2, 4, 5, 7, 8, 9},
                        {0, 6, 9, 3, 8, 5},
                        {0, 1, 2, 6, 7, 9}},
                       {{0, 4, 7, 9, 1, 5, 8, 2, 6, 3},
                        {3, 6, 8, 9, 12, 14, 15, 17, 18, 19},
                        {0, 10, 16, 19, 4, 13, 18, 7, 15, 9},
                        {0, 1, 2, 3, 10, 11, 12, 16, 17, 19}},
                       {
                           {0, 5, 9, 12, 14, 1, 6, 10, 13, 2, 7, 11, 3, 8, 4},         // OK
                           {4, 8, 11, 13, 14, 18, 21, 23, 24, 27, 29, 30, 32, 33, 34}, // OK
                           {0, 15, 25, 31, 34, 5, 19, 28, 33, 9, 22, 30, 12, 24, 14},  // OK
                           {0, 1, 2, 3, 4, 15, 16, 17, 18, 25, 26, 27, 31, 32, 34}     // OK
                       }};

static std::vector<std::vector<std::vector<int>>>
    tetrahdronEdgeNodalMap{{{0, 1}, {1, 2}, {0, 2}, {0, 3}, {1, 3}, {2, 3}},
                           {{0, 1, 2}, {2, 4, 5}, {0, 3, 5}, {0, 6, 9}, {2, 7, 9}, {5, 8, 9}},
                           {{0, 1, 2, 3},
                            {3, 6, 8, 9},
                            {0, 4, 7, 9},
                            {0, 10, 16, 19},
                            {3, 12, 17, 19},
                            {9, 15, 18, 19}},
                           {{0, 1, 2, 3, 4},
                            {4, 8, 11, 13, 14},
                            {0, 5, 9, 12, 14},
                            {0, 15, 25, 31, 34},
                            {4, 18, 27, 32, 34},
                            {14, 24, 30, 33, 34}}};

/**
 * Numbering of nodal points on hexahedron (example P3)
 *
 *  front:                                                            rear:
 *   /               /     /              /      /              /
 *  48 - 49 - 50 - 51     52 - 53 - 54 - 55     56 - 57 - 58 - 59     60 - 61 - 62 - 63
 *  |              |     /|             /|     /|             /|     /|             /|
 *  32   33   34   35     36   37   38   39     40   41   42   43     44   45   46   47
 *  |              |      |              |      |              |      |              |
 *  16   17   18   19     20   21   22   23     24   25   26   27     28   29   30   31
 *  | /            | /    | /            | /    | /            | /    |              |
 *  0 -- 1 -- 2 -- 3      4 -- 5 -- 6 -- 7      8 -- 9 -- 10 - 11     12 - 13 - 14 - 15
 *                       /              /      /              /      /              /
 */

template<typename T, int sDim, int tDim>
void LagrangeNodalPointsHexahedron(const Cell &c, vector<PointT<T, sDim, tDim>> &z, int degree) {
  if (degree == 0) {
    z[0] = PointT<T, sDim, tDim>(c());
  } else {
    int index = 0;
    for (int i = 0; i < (degree + 1); i++) {
      for (int j = 0; j < (degree + 1); j++) {
        for (int k = 0; k < (degree + 1); k++) {
          PointT<T, sDim, tDim> v = c[0] + k / T(degree) * (c[1] - c[0]);
          PointT<T, sDim, tDim> w = c[3] + k / T(degree) * (c[2] - c[3]);
          PointT<T, sDim, tDim> x = c[4] + k / T(degree) * (c[5] - c[4]);
          PointT<T, sDim, tDim> y = c[7] + k / T(degree) * (c[6] - c[7]);
          PointT<T, sDim, tDim> r = v + j / T(degree) * (w - v);
          PointT<T, sDim, tDim> s = x + j / T(degree) * (y - x);
          z[index] = r + i / T(degree) * (s - r);
          index++;
        }
      }
    }
  }
}

constexpr int LagrangeNodalPoints(CELLTYPE cellType, int degree = 1) {
  if (degree == 0) return 1;
  switch (cellType) {
  case INTERVAL:
    return degree + 1;
  case TRIANGLE:
    return ((degree + 1) * (degree + 2)) / 2;
  case QUADRILATERAL:
    return (degree + 1) * (degree + 1);
  case TETRAHEDRON:
    return ((degree + 1) * (degree + 2) * (degree + 3)) / 6;
  case HEXAHEDRON:
    return (degree + 1) * (degree + 1) * (degree + 1);
  default:
    THROW("Cell type not implemented in LagrangeDoF")
  }
}

template<typename cellT, typename T, int sDim, int tDim>
void LagrangeNodalPoints(const cellT &c, vector<PointT<T, sDim, tDim>> &z, int degree) {
  z.resize(LagrangeNodalPoints(c.ReferenceType(), degree));
  switch (c.ReferenceType()) {
  case INTERVAL:
    LagrangeNodalPointsInterval(c, z, degree);
    break;
  case TRIANGLE:
    LagrangeNodalPointsTriangle(c, z, degree);
    break;
  case QUADRILATERAL:
    LagrangeNodalPointsQuadrilateral(c, z, degree);
    break;
  case TETRAHEDRON:
    LagrangeNodalPointsTetrahedron(c, z, degree);
    break;
  case HEXAHEDRON:
    LagrangeNodalPointsHexahedron(c, z, degree);
    break;
  default:
    THROW("Unknown celltype: " + std::to_string(c.ReferenceType()))
  }
}

#endif // LAGRANGENODALPOINTS_HPP
