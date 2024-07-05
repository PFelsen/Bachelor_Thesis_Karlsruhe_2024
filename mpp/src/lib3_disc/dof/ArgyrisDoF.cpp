#include "ArgyrisDoF.hpp"

short ArgyrisDoF::NumberOfNodalPoints(const Cell &c) const { return c.Corners() + c.Edges(); }

std::vector<Point> ArgyrisDoF::GetNodalPoints(const Cell &c) const {
  std::vector<Point> sp(NumberOfNodalPoints(c));
  int n = 0;
  for (int i = 0; i < c.Corners(); ++i)
    sp[n++] = c[i];
  for (int i = 0; i < c.Edges(); ++i)
    sp[n++] = c.Edge(i);
  return sp;
}

std::vector<short> ArgyrisDoF::DoFSizesAtNodalPoints(const Cell &c) const {
  return std::vector<short>{short(6 * m), short(6 * m), short(6 * m), m, m, m};
}

std::vector<short> ArgyrisDoF::AccumulatedAllocationSizes(const Cell &c) const {
  return std::vector<short>{short(6 * m),  short(12 * m), short(18 * m),
                            short(19 * m), short(20 * m), short(21 * m)};
}

short ArgyrisDoF::NumberOfNodalPointsOnFace(const Cell &c, int faceId) const { return 3; }

short ArgyrisDoF::IdOfNodalPointOnFace(const Cell &c, int faceId, int k) const {
  int fc = c.FaceCorners(faceId);
  if (k < fc) return c.facecorner(faceId, k);
  return c.Corners() + c.faceedge(faceId, k - fc);
}

short ArgyrisDoF::NumberOfNodalPointsOnEdge(const Cell &c, int edgeId) const {
  return NumberOfNodalPointsOnFace(c, edgeId);
}

short ArgyrisDoF::IdOfNodalPointOnEdge(const Cell &c, int edgeId, int k) const {
  return IdOfNodalPointOnFace(c, edgeId, k);
}
