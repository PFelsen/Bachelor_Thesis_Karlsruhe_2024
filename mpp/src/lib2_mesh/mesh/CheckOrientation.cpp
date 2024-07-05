#include "CheckOrientation.hpp"
#include <algorithm>

// checks, if the entries corresponding to a larger dimension are set to zero
bool additionalComponentsZero(const std::vector<Point> &corners, int from) {
  for (int k = 0; k < corners.size(); ++k) {
    for (int i = from; i < corners[k].SpaceDim(); ++i) {
      if (corners[k][i] != 0.0) return false;
    }
  }
  return true;
}

void CheckOrientation::CheckInterval(std::vector<Point> &corners, std::vector<int> &cornerIndices) {
  if (corners.size() != 2) THROW("Number of corners does not fit for INTERVAL")

  if (!additionalComponentsZero(corners, 1)) THROW("Non zero entry in corners of INTERVAL")

  if (std::abs(corners[0][0] - corners[1][0]) < GeometricTolerance)
    THROW("Corners do not represent an INTERVAL")
  if (corners[0][0] > corners[1][0]) {
    std::swap(corners[0], corners[1]);
    std::swap(cornerIndices[0], cornerIndices[1]);
    Warning("Changed order of points for type INTERVAL")
  }
}

double det2D(const Point &p0, const Point &p1) {
  double d = p0[0] * p1[1] - p0[1] * p1[0];
  if (std::abs(d) < GeometricTolerance) THROW("Point are located at a line")
  return d;
}

void CheckOrientation::CheckTriangle(std::vector<Point> &corners, std::vector<int> &cornerIndices) {
  if (corners.size() != 3) THROW("Number of corners does not fit for TRIANGLE")

  if (!additionalComponentsZero(corners, 2)) THROW("Non zero entry in corners of TRIANGLE")

  if (det2D(corners[1] - corners[0], corners[2] - corners[0]) < 0) {
    std::swap(corners[1], corners[2]);
    std::swap(cornerIndices[1], cornerIndices[2]);
    Warning("Changed order of points for type TRIANGLE")
  }
}

void CheckOrientation::CheckQuadrilateral(std::vector<Point> &corners,
                                          std::vector<int> &cornerIndices) {
  if (corners.size() != 4) THROW("Number of corners does not fit for QUADRILATERAL")

  if (!additionalComponentsZero(corners, 2)) THROW("Non zero entry in corners of QUADRILATERAL")

  int cntAngle = 0;
  for (int i = 0; i < 4; ++i) {
    if (det2D(corners[i] - corners[(i + 1) % 4], corners[(i + 2) % 4] - corners[(i + 1) % 4]) > 0) {
      cntAngle++;
    }
  }

  switch (cntAngle) {
  case 0: // correctly oriented
    break;
  case 1: // non-convex quadrilateral
  case 3:
    THROW("Corners do not represent a convex QUADRILATERAL")
  case 2: { // edges cross each other
    double det1 = det2D(corners[1] - corners[0], corners[3] - corners[0]);
    double det2 = det2D(corners[0] - corners[1], corners[2] - corners[1]);
    if (det1 < 0 && det2 < 0) {
      std::swap(corners[0], corners[3]);
      std::swap(cornerIndices[0], cornerIndices[3]);
    } else if (det1 < 0 && det2 > 0) {
      std::swap(corners[0], corners[1]);
      std::swap(cornerIndices[0], cornerIndices[1]);
    } else if (det1 > 0 && det2 < 0) {
      std::swap(corners[2], corners[3]);
      std::swap(cornerIndices[2], cornerIndices[3]);
    } else {
      std::swap(corners[1], corners[2]);
      std::swap(cornerIndices[1], cornerIndices[2]);
    }
    Warning("Changed order of points for type QUADRILATERAL") break;
  }
  case 4: // conversely oriented
    std::reverse(corners.begin(), corners.end());
    std::reverse(cornerIndices.begin(), cornerIndices.end());
    Warning("Changed order of points for type QUADRILATERAL") break;
  default:
    THROW("This case should not occur")
  }
}

double det3D(const Point &p0, const Point &p1, const Point &p2) {
  double d = det(p0, p1, p2);
  if (std::abs(d) < GeometricTolerance) THROW("Point are located in a plane")
  return d;
}

void CheckOrientation::CheckTetrahedron(std::vector<Point> &corners,
                                        std::vector<int> &cornerIndices) {
  if (corners.size() != 4) THROW("Number of corners does not fit for TETRAHEDRON")

  if (det3D(corners[1] - corners[0], corners[2] - corners[0], corners[3] - corners[0]) < 0) {
    std::swap(corners[1], corners[2]);
    std::swap(cornerIndices[1], cornerIndices[2]);
    Warning("Changed order of points for type TETRAHEDRON")
  }
}

double det2D(const Point &p0, const Point &p1, const Point &normal) {
  double d = (p0 ^ p1) * normal;
  if (std::abs(d) < GeometricTolerance) THROW("Point are located in a plane")
  return d;
}

void checkOrientationSubQuad(const Point &normal, std::vector<std::pair<Point, int>> &corners) {
  int cntAngle = 0;
  for (int i = 0; i < 4; ++i) {
    if (det2D(corners[i].first - corners[(i + 1) % 4].first,
              corners[(i + 2) % 4].first - corners[(i + 1) % 4].first, normal)
        > 0) {
      cntAngle++;
    }
  }

  switch (cntAngle) {
  case 0: // correctly oriented
    break;
  case 1: // non-convex cell
  case 3:
    THROW("Corners do not represent a convex QUADRILATERAL in HEXAHEDRON")
  case 2: { // edges cross each other
    double det1 =
        det2D(corners[1].first - corners[0].first, corners[3].first - corners[0].first, normal);
    double det2 =
        det2D(corners[0].first - corners[1].first, corners[2].first - corners[1].first, normal);
    if (det1 < 0 && det2 < 0) {
      std::swap(corners[0], corners[3]);
    } else if (det1 < 0 && det2 > 0) {
      std::swap(corners[0], corners[1]);
    } else if (det1 > 0 && det2 < 0) {
      std::swap(corners[2], corners[3]);
    } else {
      std::swap(corners[1], corners[2]);
    }
    Warning("Changed order of points for type HEXAHEDRON") break;
  }
  case 4: // conversely oriented
    std::reverse(corners.begin(), corners.end());
    Warning("Changed order of points for type HEXAHEDRON") break;
  default:
    THROW("This case should not occur")
  }
}

bool checkOrientationSubQuad(const Point &c10, const Point &c20, const Point &c30) {
  // check: points are located in a plane
  if (std::abs(det(c10, c20, c30)) >= GeometricTolerance) return false;

  Point normal = c10 ^ c30;
  if ((c10 ^ c20) * normal >= GeometricTolerance && (c20 ^ c30) * normal >= GeometricTolerance)
    return true;
  return false;
}

bool correctOrientationTopToBottom(const std::vector<std::pair<Point, int>> &bottom,
                                   const std::vector<std::pair<Point, int>> &top) {
  return checkOrientationSubQuad(bottom[1].first - bottom[0].first, top[1].first - bottom[0].first,
                                 top[0].first - bottom[0].first);
}

void CheckOrientation::CheckHexahedron(std::vector<Point> &corners,
                                       std::vector<int> &cornerIndices) {
  if (corners.size() != 8) THROW("Number of corners does not fit for HEXAHEDRON")

  // first, find suitable points for bottom and top planes
  for (int i = 1; i < 8; ++i) {
    Point ci0 = corners[i] - corners[0];
    for (int j = i + 1; j < 8; ++j) {
      Point cj0 = corners[j] - corners[0];
      for (int k = j + 1; k < 8; ++k) {
        Point ck0 = corners[k] - corners[0];
        // check: corners[0], corners[i], corners[j], corners[k] are in a plane
        if (std::abs(det(ci0, cj0, ck0)) < 1e-8) {
          std::vector<int> others(7);
          std::iota(others.begin(), others.end(), 1);
          std::remove(others.begin(), others.end(), i);
          std::remove(others.begin(), others.end(), j);
          std::remove(others.begin(), others.end(), k);

          Point normal = ci0 ^ cj0;
          double tmp = normal * (corners[others[0]] - corners[0]);

          // check: remaining corners are located at one side of the plane
          bool isOnOneSide = true;
          for (int l = 1; l < 4; ++l) {
            isOnOneSide &= (tmp * (normal * (corners[others[l]] - corners[0])) > 0);
          }
          if (isOnOneSide) {
            if (tmp < 0) normal *= -1.0;
            std::vector<std::pair<Point, int>> bottom = {{corners[0], cornerIndices[0]},
                                                         {corners[i], cornerIndices[i]},
                                                         {corners[j], cornerIndices[j]},
                                                         {corners[k], cornerIndices[k]}};
            std::vector<std::pair<Point, int>> top =
                {{corners[others[0]], cornerIndices[others[0]]},
                 {corners[others[1]], cornerIndices[others[1]]},
                 {corners[others[2]], cornerIndices[others[2]]},
                 {corners[others[3]], cornerIndices[others[3]]}};

            // check: orientation of bottom and top quadrilateral
            checkOrientationSubQuad(normal, bottom);
            checkOrientationSubQuad(normal, top);

            // turn top quadrilateral if necessary
            for (int l = 0; l < 4; ++l) {
              if (correctOrientationTopToBottom(bottom, top)) {
                // check: faces are contained in a plane (note: bottom, top and front are already
                // checked)
                for (int m = 0; m < 3; ++m) {
                  Point c10 = bottom[(2 + m) % 4].first - bottom[1 + m].first;
                  Point c20 = top[(2 + m) % 4].first - bottom[1 + m].first;
                  Point c30 = top[1 + m].first - bottom[1 + m].first;

                  if (!checkOrientationSubQuad(c10, c20, c30))
                    THROW("Corners do not represent a HEXAHEDRON")
                }

                // finally, combine corners and set cornerIndices
                for (int i = 0; i < 4; ++i) {
                  corners[i] = bottom[i].first;
                  corners[i + 4] = top[i].first;
                  cornerIndices[i] = bottom[i].second;
                  cornerIndices[i + 4] = top[i].second;
                }
                return;
              }
              std::rotate(top.begin(), top.begin() + 1, top.end());
            }
            THROW("Corners do not represent a HEXAHEDRON")
          }
        }
      }
    }
  }
  THROW("Corners do not represent a HEXAHEDRON")
}

void CheckOrientation::Check(CELLTYPE type, std::vector<Point> &corners,
                             std::vector<int> &cornerIndices) {
  switch (type) {
  case INTERVAL:
    CheckInterval(corners, cornerIndices);
    break;
  case TRIANGLE:
    CheckTriangle(corners, cornerIndices);
    break;
  case QUADRILATERAL:
    CheckQuadrilateral(corners, cornerIndices);
    break;
  case TETRAHEDRON:
    CheckTetrahedron(corners, cornerIndices);
    break;
  case HEXAHEDRON:
    CheckHexahedron(corners, cornerIndices);
    break;
  default:
    THROW("checkOrientation not implemented for celltype " + std::to_string(type))
  }
}