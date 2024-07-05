#include "CurlDoF.hpp"

short CurlDoF::NumberOfNodalPoints(const Cell &c) const { return c.Edges(); }

std::vector<Point> CurlDoF::GetNodalPoints(const Cell &c) const {
  std::vector<Point> sp(NumberOfNodalPoints(c));
  for (int i = 0; i < c.Edges(); ++i)
    sp[i] = c.Edge(i);
  return sp;
}

std::vector<short> CurlDoF::DoFSizesAtNodalPoints(const Cell &c) const {
  return std::vector<short>(NumberOfNodalPoints(c), m);
}

short CurlDoF::NumberOfNodalPointsOnFace(const Cell &c, int faceId) const {
  return c.FaceEdges(faceId);
}

short CurlDoF::IdOfNodalPointOnFace(const Cell &c, int faceId, int k) const {
  return c.faceedge(faceId, k);
}

short CurlDoF::NumberOfNodalPointsOnEdge(const Cell &c, int edgeId) const { THROW("TODO") }

short CurlDoF::IdOfNodalPointOnEdge(const Cell &c, int edgeId, int k) const { THROW("TODO") }

short Curl2TetDoF::NumberOfNodalPoints(const Cell &c) const { return 20; }

std::vector<Point> Curl2TetDoF::GetNodalPoints(const Cell &c) const {
  std::vector<Point> sp(NumberOfNodalPoints(c));
  for (int i = 0; i < c.Edges(); ++i) {
    sp[2 * i] = c.Edge(i);
    sp[2 * i + 1] = 0.5 * (c.Edge(i) + c.EdgeCorner(i, (c.EdgeCorner(i, 0) < c.EdgeCorner(i, 1))));
  }
  for (int i = 0; i < c.Faces(); ++i) {
    sp[2 * (i + c.Edges())] = c.Face(i);
    sp[2 * (i + c.Edges()) + 1] = 0.5 * (c.Face(i) + c.FaceCorner(i, SmallestFaceCorner(c, i)));
  }
  return sp;
}

std::vector<short> Curl2TetDoF::DoFSizesAtNodalPoints(const Cell &c) const {
  return std::vector<short>(NumberOfNodalPoints(c), m);
}

short Curl2TetDoF::NumberOfNodalPointsOnFace(const Cell &c, int faceId) const { return 8; }

short Curl2TetDoF::IdOfNodalPointOnFace(const Cell &c, int faceId, int k) const {
  if (k < 2) return IdOfNodalPointOnEdge(c, c.faceedge(faceId, 0), k);
  else if (k < 4) return IdOfNodalPointOnEdge(c, c.faceedge(faceId, 1), k - 2);
  else if (k < 6) return IdOfNodalPointOnEdge(c, c.faceedge(faceId, 2), k - 4);
  else return 2 * (faceId + c.Edges()) + k - 6;
}

short Curl2TetDoF::NumberOfNodalPointsOnEdge(const Cell &c, int edgeId) const { return 2; }

short Curl2TetDoF::IdOfNodalPointOnEdge(const Cell &c, int edgeId, int k) const {
  return 2 * edgeId + k;
}

int Curl2TetDoF::SmallestFaceCorner(const Cell &c, int i) const {
  if (c.FaceCorner(i, 0) < c.FaceCorner(i, 1)) {
    if (c.FaceCorner(i, 0) < c.FaceCorner(i, 2)) return 0;
    return 2;
  } else if (c.FaceCorner(i, 1) < c.FaceCorner(i, 2)) return 1;
  return 2;
}

Point Curl2DoF::FaceCor(const Cell &c, int i, int j) const {
  if (i == 0) {
    switch (j) {
    case 0:
      return c[1];
    case 1:
      return c[0];
    case 2:
      return c[3];
    case 3:
      return c[2];
    }
    THROW("Not Implemented!")
    return -1;
  } else {
    return c.FaceCorner(i, j);
  }
}

short Curl2DoF::NumberOfNodalPoints(const Cell &c) const { return 54; }

std::vector<Point> Curl2DoF::GetNodalPoints(const Cell &c) const {
  std::vector<Point> sp(NumberOfNodalPoints(c));
  int n = 0;
  for (int i = 0; i < c.Edges(); ++i)
    sp[n++] = c.Edge(i);
  for (int i = 0; i < c.Edges(); ++i)
    sp[n++] = 0.5 * (c.Edge(i) + c.EdgeCorner(i, (c.EdgeCorner(i, 0) < c.EdgeCorner(i, 1))));
  for (int i = 0; i < c.Faces(); ++i) {
    Point fm = c.Face(i);
    Point f0 = FaceCor(c, i, 0);
    Point f1 = FaceCor(c, i, 1);
    sp[n++] = fm + 0.25 * (2 * (f1 > f0) - 1) * (f1 - f0);
    sp[n++] = fm + 0.25 * (FaceCor(c, i, 3) - f0);
    sp[n++] = fm - 0.25 * (FaceCor(c, i, 3) - f0);
    sp[n++] = fm;
  }
  for (int i = 0; i < 5; ++i)
    sp[n++] = 0.5 * (c() + c.Face(i));
  sp[n++] = c();
  return sp;
}

std::vector<short> Curl2DoF::DoFSizesAtNodalPoints(const Cell &c) const {
  return std::vector<short>(NumberOfNodalPoints(c), m);
}

short Curl2DoF::NumberOfNodalPointsOnFace(const Cell &c, int faceId) const { return 12; }

short Curl2DoF::IdOfNodalPointOnFace(const Cell &c, int faceId, int k) const {
  int nfe = c.FaceEdges(faceId);
  if (k < nfe) return c.faceedge(faceId, k);
  k -= nfe;
  if (k < nfe) return c.Edges() + c.faceedge(faceId, k);
  return 2 * c.Edges() + 4 * faceId + k - nfe;
}

short Curl2DoF::NumberOfNodalPointsOnEdge(const Cell &c, int edgeId) const { return 2; }

short Curl2DoF::IdOfNodalPointOnEdge(const Cell &c, int edgeId, int k) const {
  return edgeId + 12 * k;
}