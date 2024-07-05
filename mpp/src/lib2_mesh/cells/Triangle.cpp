#include "Triangle.hpp"

Triangle::Triangle(const vector<Point> &z, short sd) : Cell(z, sd) {}

CELLTYPE Triangle::Type() const { return TRIANGLE; }

CELLTYPE Triangle::ReferenceType() const { return TRIANGLE; }

const Point &Triangle::LocalCorner(int i) const { return LocalCornersTri()[i]; }

const Point &Triangle::LocalEdge(int i) const { return LocalEdgesTri()[i]; }

const Point &Triangle::LocalFace(int i) const { return LocalEdgesTri()[i]; }

const Point &Triangle::LocalFaceNormal(int i) const { return LocalFaceNormalsTri()[i]; }

const Point &Triangle::LocalCenter() const { return *LocalCenterTri(); }

Point Triangle::FaceLocalToGlobal(int face, const Point &localOnFace) const {
  return faceLocalToGlobal(face, localOnFace);
}

Point Triangle::FaceLocalToLocal(int faceId, const Point &localOnFace) const {
  return faceLocalToLocal(faceId, localOnFace);
}

Point Triangle::FaceLocalToLocalOriented(int faceId, const Point &localOnFace) const {
  return faceLocalToLocalOriented(faceId, localOnFace);
}

Point Triangle::LocalToGlobal(const Point &z) const { return localToGlobal(z); }

Transformation Triangle::GetTransformation(const Point &z) const { return getTransformation(z); }

#ifdef BUILD_IA

IAPoint Triangle::FaceLocalToGlobal(int face, const IAPoint &local) const {
  return faceLocalToGlobal(face, local);
}

IAPoint Triangle::FaceLocalToLocal(int face, const IAPoint &local) const {
  return faceLocalToLocal(face, local);
}

IAPoint Triangle::FaceLocalToLocalOriented(int faceId, const IAPoint &localOnFace) const {
  return faceLocalToLocalOriented(faceId, localOnFace);
}

IAPoint Triangle::LocalToGlobal(const IAPoint &z) const { return localToGlobal(z); }

IAPoint Triangle::GlobalToLocal(const IAPoint &z) const { return globalToLocal(z); }

IATransformation Triangle::GetTransformation(const IAPoint &z) const {
  return getTransformation(z);
}

IAInterval Triangle::LocalFaceAreaIA(int faceId) const {
  switch (faceId) {
  case 0:
  case 2:
    return IAInterval(1.0);
  case 1:
    return IAInterval::Sqrt2();
  default:
    THROW("faceId " + std::to_string(faceId) + " not possible in Triangle")
  }
}

const IAPoint *Triangle::LocalFaceNormalsIATri() {
  static const IAPoint *localFaceNormalsIA = ReferenceFaceNormals<IAPoint>(*ReferenceTriangle());
  return localFaceNormalsIA;
}

const IAPoint &Triangle::LocalFaceNormalIA(int faceId) const {
  return LocalFaceNormalsIATri()[faceId];
}

#endif

int Triangle::Edges() const { return 3; }

Point Triangle::Edge(int i) const {
  switch (i) {
  case 0:
    return 0.5 * (Corner(0) + Corner(1));
  case 1:
    return 0.5 * (Corner(1) + Corner(2));
  case 2:
    return 0.5 * (Corner(2) + Corner(0));
  }
  THROW("Not Implemented")
}

short Triangle::edgecorner(unsigned short i, unsigned short j) const {
#if DebugLevel > 0
  if (!(i < 3 && j < 2))
    THROW("edgecorner (" + to_string(i) + ", " + to_string(j) + ") does not exist")
#endif
  return (i + j) % 3;
}

int Triangle::Faces() const { return 3; }

CELLTYPE Triangle::FaceType(int face) const { return INTERVAL; }

short Triangle::FaceCorners(int i) const { return 2; }

short Triangle::facecorner(unsigned short i, unsigned short j) const {
#if DebugLevel > 0
  if (!(i < 3 && j < 2))
    THROW("facecorner (" + to_string(i) + ", " + to_string(j) + ") does not exist")
#endif
  return (i + j) % 3;
}

short Triangle::FaceEdges(int i) const { return 1; }

short Triangle::faceedge(int i, int j) const { return short(i); }

double Triangle::FaceArea(int i) const { return dist(FaceCorner(i, 0), FaceCorner(i, 1)); }

double Triangle::LocalFaceArea(int i) const {
  switch (i) {
  case 0:
    return 1;
  case 1:
    return sqrt(2.0);
  case 2:
    return 1;
  }
  THROW("Not Implemented")
}

void Triangle::GetRefinePoints(std::vector<Point> &z) const {
  z.resize(6);
  int n = 0;
  for (int i = 0; i < Corners(); ++i)
    z[n++] = Corner(i);
  for (int i = 0; i < Edges(); ++i)
    z[n++] = Edge(i);
}

const vector<Rule> &Triangle::Refine(vector<Point> &z) const {
  GetRefinePoints(z);
  return *Rtri;
}

void Triangle::GetRefineBarycentricPoints(std::vector<Point> &z, const Point &shiftCenter) const {
  z.resize(4);
  for (int i = 0; i < 3; ++i)
    z[i] = Corner(i);
  z[3] = Center() + shiftCenter;
}

const std::vector<Rule> &Triangle::RefineBarycentric(std::vector<Point> &z,
                                                     const Point &shiftCenter) const {
  GetRefineBarycentricPoints(z, shiftCenter);
  return *Rtri_barycentric;
}

int Triangle::dim() const { return 2; }

bool Triangle::plane() const { return true; }

const Cell *Triangle::ReferenceCell() const { return ReferenceTriangle(); }

int Triangle::Children() const { return 4; }

Point Triangle::Child(int i) const {
  switch (i) {
  case 0:
    return (1.0 / 6.0) * (4 * Corner(0) + Corner(1) + Corner(2));
  case 1:
    return (1.0 / 6.0) * (Corner(0) + 4 * Corner(1) + Corner(2));
  case 2:
    return (1.0 / 6.0) * (Corner(0) + Corner(1) + 4 * Corner(2));
  case 3:
    return (1.0 / 3.0) * (Corner(0) + Corner(1) + Corner(2));
  default:
    THROW("Not Implemented")
  }
}

const Triangle *Triangle::ReferenceTriangle() {
  static const Triangle *refTriangle = new Triangle(Points(3, Triangle::LocalCornersTri()), -1);
  return refTriangle;
}

const Point *Triangle::LocalCornersTri() {
  static const Point localCorners[] = {Point(0.0, 0.0, 0.0), Point(1.0, 0.0, 0.0),
                                       Point(0.0, 1.0, 0.0)};
  return localCorners;
}

const Point *Triangle::LocalEdgesTri() {
  static const Point *localEdges = ReferenceEdges(*ReferenceTriangle());
  return localEdges;
}

const Point *Triangle::LocalFaceNormalsTri() {
  static const Point *localFaceNormals = ReferenceFaceNormals<Point>(*ReferenceTriangle());
  return localFaceNormals;
}

const Point *Triangle::LocalCenterTri() {
  static const Point *localCenter = ReferenceCenter(*ReferenceTriangle());
  return localCenter;
}

Point Triangle::GlobalToLocal(const Point &z) const { return globalToLocal(z); }
