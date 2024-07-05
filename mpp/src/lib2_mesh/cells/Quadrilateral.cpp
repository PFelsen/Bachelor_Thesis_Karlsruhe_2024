#include "Quadrilateral.hpp"

Quadrilateral::Quadrilateral(const vector<Point> &z, short sd) : Cell(z, sd) {}

CELLTYPE Quadrilateral::Type() const { return QUADRILATERAL; }

CELLTYPE Quadrilateral::ReferenceType() const { return QUADRILATERAL; }

int Quadrilateral::Corners() const { return 4; }

const Point &Quadrilateral::LocalCorner(int i) const { return LocalCornersQuad()[i]; }

const Point &Quadrilateral::LocalEdge(int i) const { return LocalEdgesQuad()[i]; }

const Point &Quadrilateral::LocalFace(int i) const { return LocalEdgesQuad()[i]; }

const Point &Quadrilateral::LocalFaceNormal(int i) const { return LocalFaceNormalsQuad()[i]; }

const Point &Quadrilateral::LocalCenter() const { return *LocalCenterQuad(); }

Point Quadrilateral::FaceLocalToGlobal(int face, const Point &localOnFace) const {
  return faceLocalToGlobal(face, localOnFace);
}

Point Quadrilateral::FaceLocalToLocal(int faceId, const Point &localOnFace) const {
  return faceLocalToLocal(faceId, localOnFace);
}

Point Quadrilateral::FaceLocalToLocalOriented(int faceId, const Point &localOnFace) const {
  return faceLocalToLocalOriented(faceId, localOnFace);
}

Point Quadrilateral::LocalToGlobal(const Point &z) const { return localToGlobal(z); }

Point Quadrilateral::GlobalToLocal(const Point &z) const { return globalToLocal(z); }

Transformation Quadrilateral::GetTransformation(const Point &z) const {
  return getTransformation(z);
}

int Quadrilateral::Edges() const { return 4; }

Point Quadrilateral::Edge(int i) const {
  switch (i) {
  case 0:
    return 0.5 * (Corner(0) + Corner(1));
  case 1:
    return 0.5 * (Corner(1) + Corner(2));
  case 2:
    return 0.5 * (Corner(2) + Corner(3));
  case 3:
    return 0.5 * (Corner(3) + Corner(0));
  }
  THROW("Quadrilateral::Edge was called with i = " + std::to_string(i));
}

short Quadrilateral::edgecorner(unsigned short i, unsigned short j) const {
#if DebugLevel > 0
  if (!(i < 4 && j < 2))
    THROW("edgecorner (" + to_string(i) + ", " + to_string(j) + ") does not exist")
#endif
  return (i + j) % 4;
}

int Quadrilateral::Faces() const { return 4; }

Point Quadrilateral::Face(int i) const {
  switch (i) {
  case 0:
    return 0.5 * (Corner(0) + Corner(1));
  case 1:
    return 0.5 * (Corner(1) + Corner(2));
  case 2:
    return 0.5 * (Corner(2) + Corner(3));
  case 3:
    return 0.5 * (Corner(3) + Corner(0));
  }
  THROW("Not Implemented")
}

CELLTYPE Quadrilateral::FaceType(int face) const { return INTERVAL; }

short Quadrilateral::FaceCorners(int i) const { return 2; }

short Quadrilateral::facecorner(unsigned short i, unsigned short j) const {
#if DebugLevel > 0
  if (!(i < 4 && j < 2))
    THROW("facecorner (" + to_string(i) + ", " + to_string(j) + ") does not exist")
#endif
  return (i + j) % 4;
}

short Quadrilateral::FaceEdges(int i) const { return 1; }

short Quadrilateral::faceedge(int i, int j) const { return short(i); }

double Quadrilateral::LocalFaceArea(int i) const { return 1; }

double Quadrilateral::FaceArea(int i) const { return dist(FaceCorner(i, 0), FaceCorner(i, 1)); }

void Quadrilateral::GetRefinePoints(std::vector<Point> &z) const {
  z.resize(9);
  int n = 0;
  for (int i = 0; i < Corners(); ++i)
    z[n++] = Corner(i);
  for (int i = 0; i < Edges(); ++i)
    z[n++] = Edge(i);
  z[n] = Center();
}

const vector<Rule> &Quadrilateral::Refine(vector<Point> &z) const {
  GetRefinePoints(z);
  return *Rquad;
}

void Quadrilateral::GetRefineBarycentricPoints(std::vector<Point> &z,
                                               const Point &shiftCenter) const {
  z.resize(5);
  for (int i = 0; i < 4; ++i)
    z[i] = Corner(i);
  z[4] = Center() + shiftCenter;
}

const std::vector<Rule> &Quadrilateral::RefineBarycentric(std::vector<Point> &z,
                                                          const Point &shiftCenter) const {
  GetRefineBarycentricPoints(z, shiftCenter);
  return *Rquad_barycentric;
}

int Quadrilateral::dim() const { return 2; }

bool Quadrilateral::plane() const { return true; }

const Cell *Quadrilateral::ReferenceCell() const { return ReferenceQuadrilateral(); }

int Quadrilateral::Children() const { return 4; }

Point Quadrilateral::Child(int i) const {
  switch (i) {
  case 0:
    return 0.0625 * (9 * Corner(0) + 3 * Corner(1) + Corner(2) + 3 * Corner(3));
  case 1:
    return 0.0625 * (3 * Corner(0) + 9 * Corner(1) + 3 * Corner(2) + Corner(3));
  case 2:
    return 0.0625 * (Corner(0) + 3 * Corner(1) + 9 * Corner(2) + 3 * Corner(3));
  case 3:
    return 0.0625 * (3 * Corner(0) + Corner(1) + 3 * Corner(2) + 9 * Corner(3));
  }
  THROW("Not Implemented")
}

const Quadrilateral *Quadrilateral::ReferenceQuadrilateral() {
  static const Quadrilateral *refQuad =
      new Quadrilateral(Points(4, Quadrilateral::LocalCornersQuad()), -1);
  return refQuad;
}

const Point *Quadrilateral::LocalCornersQuad() {

  static const Point localCornersQ[] = {Point(0.0, 0.0), Point(1.0, 0.0), Point(1.0, 1.0),
                                        Point(0.0, 1.0)};
  return localCornersQ;
}

const Point *Quadrilateral::LocalEdgesQuad() {

  static const Point *localEdges = ReferenceEdges(*ReferenceQuadrilateral());
  return localEdges;
}

const Point *Quadrilateral::LocalFaceNormalsQuad() {

  static const Point *localFaceN = ReferenceFaceNormals<Point>(*ReferenceQuadrilateral());
  return localFaceN;
}

const Point *Quadrilateral::LocalCenterQuad() {
  static const Point *localCenter = ReferenceCenter(*ReferenceQuadrilateral());
  return localCenter;
}

bool Quadrilateral::PointInCell(const Point &z) const {
  Point loc = GlobalToLocal(z);
  if (loc[0] >= 0 && loc[0] <= 1 && loc[1] >= 0 && loc[1] <= 1) return true;
  else return false;
}

#ifdef BUILD_IA

IAPoint Quadrilateral::FaceLocalToGlobal(int face, const IAPoint &local) const {
  return faceLocalToGlobal(face, local);
}

IAPoint Quadrilateral::FaceLocalToLocal(int face, const IAPoint &local) const {
  return faceLocalToLocal(face, local);
}

IAPoint Quadrilateral::FaceLocalToLocalOriented(int faceId, const IAPoint &localOnFace) const {
  return faceLocalToLocalOriented(faceId, localOnFace);
}

IAPoint Quadrilateral::LocalToGlobal(const IAPoint &z) const { return localToGlobal(z); }

IAPoint Quadrilateral::GlobalToLocal(const IAPoint &z) const { return globalToLocal(z); }

IATransformation Quadrilateral::GetTransformation(const IAPoint &z) const {
  return getTransformation(z);
}

IAInterval Quadrilateral::LocalFaceAreaIA(int i) const { return IAInterval(1.0); }

const IAPoint *Quadrilateral::LocalFaceNormalsIAQuad() {
  static const IAPoint *localFaceNormalsIA =
      ReferenceFaceNormals<IAPoint>(*ReferenceQuadrilateral());
  return localFaceNormalsIA;
}

const IAPoint &Quadrilateral::LocalFaceNormalIA(int faceId) const {
  return LocalFaceNormalsIAQuad()[faceId];
}

#endif


Quadrilateral2::Quadrilateral2(const vector<Point> &z, short sd) : Quadrilateral(z, sd) {
  int n = 4;
  for (int i = 0; i < Edges(); ++i)
    corners[n++] = Edge(i);
}

Point Quadrilateral2::LocalToGlobal(const Point &z) const {
  return ((1 - z[0]) * (1 - z[1]) * (1 - 2 * z[0] - 2 * z[1])) * corners[0]
         + (-z[0] * (1 - z[1]) * (1 - 2 * z[0] + 2 * z[1])) * corners[1]
         + (-z[0] * z[1] * (3 - 2 * z[0] - 2 * z[1])) * corners[2]
         + (-z[1] * (1 - z[0]) * (1 + 2 * z[0] - 2 * z[1])) * corners[3]
         + (4 * z[0] * (1 - z[0]) * (1 - z[1])) * corners[4]
         + (4 * z[0] * z[1] * (1 - z[1])) * corners[5] + (4 * z[0] * z[1] * (1 - z[0])) * corners[6]
         + (4 * z[1] * (1 - z[0]) * (1 - z[1])) * corners[7];
}

Transformation Quadrilateral2::GetTransformation(const Point &z) const {
  vector<Point> y(2);
  for (int k = 0; k < y.size(); ++k) {
    y[k] = Point((1 - z[1]) * (-3 + 4 * z[0] + 2 * z[1]), (1 - z[0]) * (-3 + 4 * z[1] + 2 * z[0]))
               * corners[0][k]
           + Point((1 - z[1]) * (-1 + 4 * z[0] - 2 * z[1]), -z[0] * (1 + 2 * z[0] - 4 * z[1]))
                 * corners[1][k]
           + Point(-z[1] * (3 - 4 * z[0] - 2 * z[1]), -z[0] * (3 - 4 * z[1] - 2 * z[0]))
                 * corners[2][k]
           + Point(-z[1] * (1 - 4 * z[0] + 2 * z[1]), (z[0] - 1) * (1 + 2 * z[0] - 4 * z[1]))
                 * corners[3][k]
           + Point(4 * (1 - 2 * z[0]) * (1 - z[1]), -4 * z[0] * (1 - z[0])) * corners[4][k]
           + Point(4 * z[1] * (1 - z[1]), 4 * z[0] * (1 - 2 * z[1])) * corners[5][k]
           + Point(4 * z[1] * (1 - 2 * z[0]), 4 * z[0] * (1 - z[0])) * corners[6][k]
           + Point(-4 * z[1] * (1 - z[1]), 4 * (1 - z[0]) * (1 - 2 * z[1])) * corners[7][k];
  }
  return Transformation(y);
}

CELLTYPE Quadrilateral2::Type() const { return QUADRILATERAL2; }

void Quadrilateral2::GetRefinePoints(std::vector<Point> &z) const {
  z.resize(21);
  int n = 0;
  for (int i = 0; i < Corners(); ++i)
    z[n++] = Corner(i);
  for (int i = 0; i < Edges(); ++i)
    z[n++] = LocalToGlobal(LocalEdge(i));
  z[n++] = LocalToGlobal(LocalCenter());
  for (int i = 0; i < Edges(); ++i) {
    z[n++] = LocalToGlobal(0.5 * (LocalCorner(edgecorner(i, 0)) + LocalEdge(i)));
    z[n++] = LocalToGlobal(0.5 * (LocalCorner(edgecorner(i, 1)) + LocalEdge(i)));
    z[n++] = LocalToGlobal(0.5 * (LocalCenter() + LocalEdge(i)));
  }
}

const vector<Rule> &Quadrilateral2::Refine(vector<Point> &z) const {
  GetRefinePoints(z);
  return *R2quad;
}

const Quadrilateral2 *Quadrilateral2::ReferenceQuadrilateral2() {
  static const Quadrilateral2 *refQuad2 =
      new Quadrilateral2(Points(4, Quadrilateral::LocalCornersQuad(), 4,
                                Quadrilateral::LocalEdgesQuad()),
                         -1);
  return refQuad2;
}
