#include "Intval.hpp"

Intval::Intval(const vector<Point> &x, int sd) : Cell(x, sd) {}

CELLTYPE Intval::Type() const { return INTERVAL; }

CELLTYPE Intval::ReferenceType() const { return INTERVAL; }

int Intval::Corners() const { return 2; }

const Point &Intval::LocalCorner(int i) const { return LocalCornersInt()[i]; }

const Point &Intval::LocalEdge(int i) const { return LocalEdgesInt()[i]; }

const Point &Intval::LocalFace(int i) const { return LocalEdgesInt()[i]; }

const Point &Intval::LocalFaceNormal(int i) const { return LocalFaceNormalsInt()[i]; }

const Point &Intval::LocalCenter() const { return *LocalCenterInt(); }

Point Intval::FaceLocalToGlobal(int face, const Point &localOnFace) const {
  return faceLocalToGlobal(face, localOnFace);
}

Point Intval::FaceLocalToLocal(int faceId, const Point &localOnFace) const {
  return faceLocalToLocal(faceId, localOnFace);
}

Point Intval::FaceLocalToLocalOriented(int faceId, const Point &localOnFace) const {
  return faceLocalToLocalOriented(faceId, localOnFace);
}

Point Intval::LocalToGlobal(const Point &z) const { return localToGlobal(z); }

Transformation Intval::GetTransformation(const Point &z) const { return getTransformation(z); }

#ifdef BUILD_IA

IAPoint Intval::FaceLocalToGlobal(int face, const IAPoint &local) const {
  return faceLocalToGlobal(face, local);
}

IAPoint Intval::FaceLocalToLocal(int face, const IAPoint &local) const {
  return faceLocalToLocal(face, local);
}

IAPoint Intval::FaceLocalToLocalOriented(int faceId, const IAPoint &localOnFace) const {
  return faceLocalToLocalOriented(faceId, localOnFace);
}

IAPoint Intval::LocalToGlobal(const IAPoint &z) const { return localToGlobal(z); }

IAPoint Intval::GlobalToLocal(const IAPoint &z) const { return globalToLocal(z); }

IATransformation Intval::GetTransformation(const IAPoint &z) const { return getTransformation(z); }

IAInterval Intval::LocalFaceAreaIA(int i) const { return IAInterval(1.0); }

const IAPoint *Intval::LocalFaceNormalsIAInt() {
  static const IAPoint *localFaceNormalsIA = ReferenceFaceNormals<IAPoint>(*ReferenceIntval());
  return localFaceNormalsIA;
}

const IAPoint &Intval::LocalFaceNormalIA(int faceId) const {
  return LocalFaceNormalsIAInt()[faceId];
}

#endif

int Intval::Edges() const { return 2; }

Point Intval::Edge(int i) const { return Corner(i); }

short Intval::edgecorner(unsigned short edge, unsigned short corner) const {
#if DebugLevel > 0
  if (edge > 1)
    THROW("edgecorner(" + to_string(edge) + ", " + to_string(corner) + ") does not exist")
#endif
  return edge;
}

int Intval::Faces() const { return 2; }

Point Intval::Face(int i) const { return Corner(i); }

CELLTYPE Intval::FaceType(int face) const { return POINT; }

short Intval::FaceCorners(int i) const { return 1; }

short Intval::facecorner(unsigned short face, unsigned short corner) const { return face; }

short Intval::FaceEdges(int i) const { return 1; }

short Intval::faceedge(int i, int j) const { return i; }

double Intval::FaceArea(int i) const { return 1.0; }

double Intval::LocalFaceArea(int i) const { return 1.0; }

void Intval::GetRefinePoints(std::vector<Point> &a) const {
  a.resize(3);
  a[0] = Corner(0);
  a[1] = Corner(1);
  a[2] = 0.5 * (Corner(0) + Corner(1));
}

const vector<Rule> &Intval::Refine(vector<Point> &a) const {
  GetRefinePoints(a);
  return *Rint;
}

int Intval::dim() const { return 1; }

bool Intval::plane() const { return true; }

const Cell *Intval::ReferenceCell() const { return ReferenceIntval(); }

int Intval::Children() const { return 2; }

Point Intval::Child(int i) const {
  switch (i) {
  case 0:
    return (3.0 * Corner(0) + Corner(1)) / 4.0;
  case 1:
    return (Corner(0) + 3.0 * Corner(1)) / 4.0;
  }
  THROW("Not Implemented")
}

const Intval *Intval::ReferenceIntval() {
  static const Intval *refIntval = new Intval(Points(2, Intval::LocalCornersInt()), -1);
  return refIntval;
}

const Point *Intval::LocalCornersInt() {
  static const Point localCorners[] = {Point(0.0, 0.0, 0.0), Point(1.0, 0.0, 0.0)};
  return localCorners;
}

const Point *Intval::LocalEdgesInt() {
  static const Point *localEdges = LocalCornersInt();
  return localEdges;
}

const Point *Intval::LocalFaceNormalsInt() {

  static const Point localFaceNormals[] = {Point(-1.0, 0.0, 0.0), Point(1.0, 0.0, 0.0)};
  return localFaceNormals;
}

const Point *Intval::LocalCenterInt() {
  static const Point *localCenter = ReferenceCenter(*Intval::ReferenceIntval());
  return localCenter;
}

Point Intval::GlobalToLocal(const Point &z) const { return globalToLocal(z); }

bool Intval::PointInCell(const Point &z) const {
  return (abs(z[0] - Center()[0]) < 0.5 * abs(Corner(0)[0] - Corner(1)[0]) + GeometricTolerance);
}
