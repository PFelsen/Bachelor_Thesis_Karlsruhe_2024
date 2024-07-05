#include "DGDoF.hpp"
#include "LagrangeNodalPoints.hpp"

#include <ext/numeric>

short DGDoF::NumberOfNodalPoints(const Cell &c) const {
  return EquidistantNodalPointCount(c.ReferenceType(), cellStandardDegree);
}

std::vector<Point> DGDoF::GetNodalPoints(const Cell &c) const {
  return EquidistantNodalPoints(c, cellStandardDegree);
}

std::vector<short> DGDoF::DoFSizesAtNodalPoints(const Cell &c) const {
  return std::vector<short>(NumberOfNodalPoints(c), m);
}

short DGDoF::NumberOfNodalPointsOnFace(const Cell &c, int faceId) const {
  return EquidistantNodalPointCountOnFace(c.ReferenceType(), cellStandardDegree, faceId);
}

short DGDoF::IdOfNodalPointOnFace(const Cell &c, int faceId, int k) const {
  return EquidistantNodalPointIdOnFace(c.ReferenceType(), cellStandardDegree, faceId, k);
}

short DGDoF::NumberOfNodalPointsOnEdge(const Cell &c, int edgeId) const {
  return EquidistantNodalPointCountOnEdge(c.ReferenceType(), cellStandardDegree, edgeId);
}

short DGDoF::IdOfNodalPointOnEdge(const Cell &c, int edgeId, int k) const {
  return EquidistantNodalPointIdOnEdge(c.ReferenceType(), cellStandardDegree, edgeId, k);
}

short DGDoF::NumberOfStoragePoints(const Cell &c) const { return 1; }

std::vector<Point> DGDoF::GetStoragePoints(const Cell &c) const { return std::vector<Point>{c()}; }

std::vector<short> DGDoF::AllocationSizesAtStoragePoints(const Cell &c) const {
  return std::vector<short>{
      short(m * EquidistantNodalPointCount(c.ReferenceType(), cellStandardDegree))};
}

short DGDoF::NumberOfStoragePointsOnFace(const Cell &c, int faceId) const { return 0; }

short DGDoF::IdOfStoragePointOnFace(const Cell &c, int faceId, int k) const {
  THROW("Should never reach this function")
}

short DGDoF::NumberOfStoragePointsOnEdge(const Cell &c, int edgeId) const { return 0; }

short DGDoF::IdOfStoragePointOnEdge(const Cell &c, int edgeId, int k) const {
  THROW("Should never reach this function")
}

void DGDoF::set_cell_deg(const Cell &c, int deg) {
  vector<short> *c_deg;
  Point z = c();
  if (cellDegree.find(z) != cellDegree.end()) {
    c_deg = &cellDegree[z];
  } else {
    c_deg = &cellDegree[z];
    c_deg->resize(1, -1);
  }
  (*c_deg)[0] = deg; // Todo vorher [ud]
}