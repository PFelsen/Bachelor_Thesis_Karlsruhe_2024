#include "SerendipityDoF.hpp"

short SerendipityDoF::NumberOfNodalPoints(const Cell &c) const {
  switch (c.ReferenceType()) {
  case QUADRILATERAL:
    return 8;
  case HEXAHEDRON:
    return 20;
  default:
    return EquidistantNodalPointCount(c.ReferenceType(), 2);
  }
}

std::vector<Point> SerendipityDoF::GetNodalPoints(const Cell &c) const {
  switch (c.ReferenceType()) {
  case QUADRILATERAL:
  case HEXAHEDRON: {
    std::vector<Point> z(NumberOfNodalPoints(c));
    int n = 0;
    for (int i = 0; i < c.Corners(); ++i)
      z[n++] = c[i];
    for (int i = 0; i < c.Edges(); ++i)
      z[n++] = c.Edge(i);
    return z;
  }
  default:
    return EquidistantNodalPoints(c, 2);
  }
}

std::vector<short> SerendipityDoF::DoFSizesAtNodalPoints(const Cell &c) const {
  return std::vector<short>(NumberOfNodalPoints(c), m);
}

short SerendipityDoF::NumberOfNodalPointsOnFace(const Cell &c, int faceId) const {
  switch (c.ReferenceType()) {
  case QUADRILATERAL:
  case HEXAHEDRON:
    return c.FaceCorners(faceId) + c.FaceEdges(faceId);
  default:
    return EquidistantNodalPointCountOnFace(c.ReferenceType(), 2, faceId);
  }
}

short SerendipityDoF::IdOfNodalPointOnFace(const Cell &c, int faceId, int k) const {
  switch (c.ReferenceType()) {
  case QUADRILATERAL:
  case HEXAHEDRON: {
    int fc = c.FaceCorners(faceId);
    if (k < fc) return c.facecorner(faceId, k);
    return c.Corners() + c.faceedge(faceId, k - fc);
  }
  default:
    return EquidistantNodalPointIdOnFace(c.ReferenceType(), 2, faceId, k);
  }
}

short SerendipityDoF::NumberOfNodalPointsOnEdge(const Cell &c, int edgeId) const {
  switch (c.ReferenceType()) {
  case QUADRILATERAL:
  case HEXAHEDRON:
    return c.EdgeCorners(edgeId);
  default:
    return EquidistantNodalPointCountOnEdge(c.ReferenceType(), 2, edgeId);
  }
}

short SerendipityDoF::IdOfNodalPointOnEdge(const Cell &c, int edgeId, int k) const {
  switch (c.ReferenceType()) {
  case QUADRILATERAL:
  case HEXAHEDRON:
    return c.edgecorner(edgeId, k);
  default:
    return EquidistantNodalPointIdOnEdge(c.ReferenceType(), 2, edgeId, k);
  }
}