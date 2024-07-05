#include "RTDoF.hpp"

short RTDoF::NumberOfNodalPoints(const Cell &c) const {
  switch (c.ReferenceType()) {
  case TRIANGLE:
    return 3 * (order + 1) + (order * (order + 1)) / 2;
  default:
    if (order == 0) return c.Faces();
    THROW("Cell type not implemented in RTDoF")
  }
}

std::vector<Point> RTDoF::GetNodalPoints(const Cell &c) const {
  switch (c.ReferenceType()) {
  case TRIANGLE:
    return RTNodalPointsTriangle<Point>(c, order);
  default:
    if (order == 0) {
      std::vector<Point> z(c.Faces());
      for (int i = 0; i < c.Faces(); ++i)
        z[i] = c.Face(i);
      return z;
    } else {
      THROW("Cell type not implemented in RTDoF")
    }
  }
}

std::vector<short> RTDoF::DoFSizesAtNodalPoints(const Cell &c) const {
  switch (c.ReferenceType()) {
  case TRIANGLE: {
    std::vector<short> d(3 * (order + 1), m);
    for (int i = 3 * (order + 1); i < NumberOfNodalPoints(c); ++i)
      d.emplace_back(short(2 * m));
    return d;
  }
  default:
    return std::vector<short>(c.Faces(), m);
  }
}

short RTDoF::NumberOfNodalPointsOnFace(const Cell &c, int faceId) const {
  switch (c.ReferenceType()) {
  case TRIANGLE:
    return order + 1;
  default:
    return 1;
  }
}

short RTDoF::IdOfNodalPointOnFace(const Cell &c, int faceId, int k) const {
  switch (c.ReferenceType()) {
  case TRIANGLE:
    switch (faceId) {
    case 0:
      return 3 * k;
    case 1:
      return 3 * k + 1;
    case 2:
      return 3 * k + 2;
    }
    break;
  default:
    return faceId;
  }
  return faceId;
}

short RTDoF::NumberOfNodalPointsOnEdge(const Cell &c, int edgeId) const {
  switch (c.ReferenceType()) {
  case TRIANGLE:
    return order + 1;
  default:
    return 1;
  }
}

short RTDoF::IdOfNodalPointOnEdge(const Cell &c, int edgeId, int k) const {
  switch (c.ReferenceType()) {
  case TRIANGLE:
    switch (edgeId) {
    case 0:
      return 3 * k;
    case 1:
      return 3 * k + 1;
    case 2:
      return 3 * k + 2;
    }
    break;
  default:
    return edgeId;
  }
  return edgeId;
}