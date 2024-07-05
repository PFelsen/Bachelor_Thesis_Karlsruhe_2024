#include "Tetrahedron.hpp"

Tetrahedron::Tetrahedron(const vector<Point> &z, int sd) : Cell(z, sd) {}

CELLTYPE Tetrahedron::Type() const { return TETRAHEDRON; }

CELLTYPE Tetrahedron::ReferenceType() const { return TETRAHEDRON; }

int Tetrahedron::Corners() const { return 4; }

const Point &Tetrahedron::LocalCorner(int i) const { return LocalCornersTet()[i]; }

const Point &Tetrahedron::LocalEdge(int i) const { return LocalEdgesTet()[i]; }

const Point &Tetrahedron::LocalFace(int i) const { return LocalFacesTet()[i]; }

const Point &Tetrahedron::LocalFaceNormal(int i) const { return LocalFaceNormalsTet()[i]; }

const Point &Tetrahedron::LocalCenter() const { return *LocalCenterTet(); }

Point Tetrahedron::FaceLocalToGlobal(int face, const Point &localOnFace) const {
  return faceLocalToGlobal(face, localOnFace);
}

Point Tetrahedron::FaceLocalToLocal(int faceId, const Point &localOnFace) const {
  return faceLocalToLocal(faceId, localOnFace);
}

Point Tetrahedron::FaceLocalToLocalOriented(int faceId, const Point &localOnFace) const {
  return faceLocalToLocalOriented(faceId, localOnFace);
}

Point Tetrahedron::LocalToGlobal(const Point &z) const { return localToGlobal(z); }

Transformation Tetrahedron::GetTransformation(const Point &z) const { return getTransformation(z); }

int Tetrahedron::Edges() const { return 6; }

Point Tetrahedron::Edge(int i) const {
  switch (i) {
  case 0:
    return 0.5 * (Corner(0) + Corner(1));
  case 1:
    return 0.5 * (Corner(1) + Corner(2));
  case 2:
    return 0.5 * (Corner(0) + Corner(2));
  case 3:
    return 0.5 * (Corner(0) + Corner(3));
  case 4:
    return 0.5 * (Corner(1) + Corner(3));
  case 5:
    return 0.5 * (Corner(2) + Corner(3));
  }
  THROW("Not Implemented")
}

short Tetrahedron::edgecorner(unsigned short i, unsigned short j) const {
  switch (i) {
  case 0:
    switch (j) {
    case 0:
      return 0;
    case 1:
      return 1;
    }
  case 1:
    switch (j) {
    case 0:
      return 1;
    case 1:
      return 2;
    }
  case 2:
    switch (j) {
    case 0:
      return 0;
    case 1:
      return 2;
    }
  case 3:
    switch (j) {
    case 0:
      return 0;
    case 1:
      return 3;
    }
  case 4:
    switch (j) {
    case 0:
      return 1;
    case 1:
      return 3;
    }
  case 5:
    switch (j) {
    case 0:
      return 2;
    case 1:
      return 3;
    }
  }
  THROW("edgecorner (" + to_string(i) + ", " + to_string(j) + ") does not exist")
}

int Tetrahedron::Faces() const { return 4; }

Point Tetrahedron::Face(int i) const {
  switch (i) {
  case 0:
    return (1.0 / 3.0) * (Corner(0) + Corner(2) + Corner(1));
  case 1:
    return (1.0 / 3.0) * (Corner(1) + Corner(2) + Corner(3));
  case 2:
    return (1.0 / 3.0) * (Corner(0) + Corner(3) + Corner(2));
  case 3:
    return (1.0 / 3.0) * (Corner(0) + Corner(1) + Corner(3));
  }
  THROW("Not Implemented")
}

CELLTYPE Tetrahedron::FaceType(int face) const { return TRIANGLE; }

short Tetrahedron::FaceCorners(int i) const { return 3; }

short Tetrahedron::facecorner(unsigned short face, unsigned short corner) const {
  switch (face) {
  case 0:
    switch (corner) {
    case 0:
      return 0;
    case 1:
      return 2;
    case 2:
      return 1;
    }
  case 1:
    switch (corner) {
    case 0:
      return 1;
    case 1:
      return 2;
    case 2:
      return 3;
    }
  case 2:
    switch (corner) {
    case 0:
      return 0;
    case 1:
      return 3;
    case 2:
      return 2;
    }
  case 3:
    switch (corner) {
    case 0:
      return 0;
    case 1:
      return 1;
    case 2:
      return 3;
    }
  }
  THROW("facecorner (" + to_string(face) + ", " + to_string(corner) + ") does not exist")
}

short Tetrahedron::FaceEdges(int i) const { return 3; }

short Tetrahedron::faceedge(int i, int j) const {
  switch (i) {
  case 0:
    switch (j) {
    case 0:
      return 2;
    case 1:
      return 1;
    case 2:
      return 0;
    }
  case 1:
    switch (j) {
    case 0:
      return 1;
    case 1:
      return 5;
    case 2:
      return 4;
    }
  case 2:
    switch (j) {
    case 0:
      return 3;
    case 1:
      return 5;
    case 2:
      return 2;
    }
  case 3:
    switch (j) {
    case 0:
      return 0;
    case 1:
      return 4;
    case 2:
      return 3;
    }
  }
  THROW("faceedge (" + to_string(i) + ", " + to_string(j) + ") does not exist")
}

double Tetrahedron::LocalFaceArea(int i) const {
  switch (i) {
  case 0:
    return 0.5;
  case 1:
    return sqrt(0.75);
  case 2:
    return 0.5;
  case 3:
    return 0.5;
  }
  THROW("Not Implemented")
}

double Tetrahedron::GetVolume() const {
  double volume(0.0);
  /*const double p = (5 - sqrt(5.0)) / 20;
  const double q = (5 + 3 * sqrt(5.0)) / 20;
  Point x(p,p,p);
  Point y(q,p,p);
  Point Corner(p,q,p);
  Point b(p,p,q);
  volume += 1.0*det(y,z,b);
  volume += (-1.0)*det(x,z,b);
  volume += 1.0*det(x,y,b);
  volume += (-1.0)*det(x,y,z);*/
  volume += 1.0 * det(Corner(1), Corner(2), Corner(3));
  volume += (-1.0) * det(Corner(0), Corner(2), Corner(3));
  volume += 1.0 * det(Corner(0), Corner(1), Corner(3));
  volume += (-1.0) * det(Corner(0), Corner(1), Corner(2));

  return (1.0 / 6.0) * volume; //
}

VectorField Tetrahedron::GetNablaNi(int i) const {
  Scalar a2, a3, a4 = 0.0;
  Scalar b2, b3, b4 = 0.0;
  Scalar c2, c3, c4 = 0.0;
  Scalar v2, v3, v4 = 0.0;

  Point P1(1.0, Corner(0)[1], Corner(0)[2]);
  Point P2(1.0, Corner(1)[1], Corner(1)[2]);
  Point P3(1.0, Corner(2)[1], Corner(2)[2]);
  Point P4(1.0, Corner(3)[1], Corner(3)[2]);
  v2 = (-1) * det(P2, P3, P4);
  a2 = 1 * det(P1, P3, P4);
  b2 = (-1) * det(P1, P2, P4);
  c2 = 1 * det(P1, P2, P3);

  Point p1(1.0, Corner(0)[0], Corner(0)[2]);
  Point p2(1.0, Corner(1)[0], Corner(1)[2]);
  Point p3(1.0, Corner(2)[0], Corner(2)[2]);
  Point p4(1.0, Corner(3)[0], Corner(3)[2]);
  v3 = (1) * det(p2, p3, p4);
  a3 = (-1) * det(p1, p3, p4);
  b3 = 1 * det(p1, p2, p4);
  c3 = (-1) * det(p1, p2, p3);

  Point PP1(1.0, Corner(0)[0], Corner(0)[1]);
  Point PP2(1.0, Corner(1)[0], Corner(1)[1]);
  Point PP3(1.0, Corner(2)[0], Corner(2)[1]);
  Point PP4(1.0, Corner(3)[0], Corner(3)[1]);
  v4 = (-1) * det(PP2, PP3, PP4);
  a4 = 1 * det(PP1, PP3, PP4);
  b4 = (-1) * det(PP1, PP2, PP4);
  c4 = 1 * det(PP1, PP2, PP3);

  if (i == 0) {
    VectorField NablaNi(v2, v3, v4);
    return (1.0 / (6.0 * GetVolume())) * NablaNi;
  } else if (i == 1) {
    VectorField NablaNi(a2, a3, a4);
    return (1.0 / (6.0 * GetVolume())) * NablaNi;
  } else if (i == 2) {
    VectorField NablaNi(b2, b3, b4);
    return (1.0 / (6.0 * GetVolume())) * NablaNi;
  } else if (i == 3) {
    VectorField NablaNi(c2, c3, c4);
    return (1.0 / (6.0 * GetVolume())) * NablaNi;
  } else {
    THROW("Number out of range!!")
  }
}

double Tetrahedron::FaceArea(int i) const {
  return 0.5 * norm((FaceCorner(i, 1) - FaceCorner(i, 0)) ^ (FaceCorner(i, 2) - FaceCorner(i, 0)));
}

void Tetrahedron::GetRefinePoints(std::vector<Point> &z) const {
  z.resize(10);
  int n = 0;
  for (int i = 0; i < Corners(); ++i)
    z[n++] = Corner(i);
  for (int i = 0; i < Edges(); ++i)
    z[n++] = Edge(i);
}

const vector<Rule> &Tetrahedron::Refine(vector<Point> &z) const {
  GetRefinePoints(z);
  double d0 = norm(z[7] - z[5]);
  double d1 = norm(z[6] - z[8]);
  double d2 = norm(z[9] - z[4]);
  if (d0 < d1) {
    if (d0 < d2) return *Rtet;
    else return *Rtet_2;
  } else {
    if (d1 < d2) return *Rtet_1;
    else return *Rtet_2;
  }
}

void Tetrahedron::GetRefineBarycentricPoints(std::vector<Point> &z,
                                             const Point &shiftCenter) const {
  z.resize(5);
  for (int i = 0; i < 4; ++i)
    z[i] = Corner(i);
  z[4] = Center() + shiftCenter;
}

const std::vector<Rule> &Tetrahedron::RefineBarycentric(std::vector<Point> &z,
                                                        const Point &shiftCenter) const {
  GetRefineBarycentricPoints(z, shiftCenter);
  return *Rtet_barycentric;
}

int Tetrahedron::dim() const { return 3; }

bool Tetrahedron::plane() const { return false; }

const Cell *Tetrahedron::ReferenceCell() const { return ReferenceTetrahedron(); }

int Tetrahedron::Children() const { return 8; }

Point Tetrahedron::Child(int i) const {
  switch (i) {
  case 0:
    return 0.125 * (5 * Corner(0) + Corner(1) + Corner(2) + Corner(3));
  case 1:
    return 0.125 * (Corner(0) + 5 * Corner(1) + Corner(2) + Corner(3));
  case 2:
    return 0.125 * (Corner(0) + Corner(1) + 5 * Corner(2) + Corner(3));
  case 3:
    return 0.125 * (Corner(0) + Corner(1) + Corner(2) + 5 * Corner(3));
  case 4:
  case 5:
  case 6:
  case 7:
    vector<Point> y(10);
    const vector<Rule> &R = Refine(y);
    return 0.25 * (y[R[i][0]] + y[R[i][1]] + y[R[i][2]] + y[R[i][3]]);
  }
  THROW("Not Implemented")
}

const Tetrahedron *Tetrahedron::ReferenceTetrahedron() {
  static const Tetrahedron *refTetra =
      new Tetrahedron(Points(4, Tetrahedron::LocalCornersTet()), -1);
  return refTetra;
}

const Point *Tetrahedron::LocalCornersTet() {
  static const Point localCorners[] = {Point(0.0, 0.0, 0.0), Point(1.0, 0.0, 0.0),
                                       Point(0.0, 1.0, 0.0), Point(0.0, 0.0, 1.0)};
  return localCorners;
}

const Point *Tetrahedron::LocalEdgesTet() {
  static const Point *localEdges = ReferenceEdges(*ReferenceTetrahedron());
  return localEdges;
}

const Point *Tetrahedron::LocalFacesTet() {
  static const Point *localFaces = ReferenceFaces(*ReferenceTetrahedron());
  return localFaces;
}

const Point *Tetrahedron::LocalFaceNormalsTet() {
  static const Point *localFaceNormals = ReferenceFaceNormals<Point>(*ReferenceTetrahedron());
  return localFaceNormals;
}

const Point *Tetrahedron::LocalCenterTet() {
  static const Point *localCenter = ReferenceCenter(*ReferenceTetrahedron());
  return localCenter;
}

Point Tetrahedron::GlobalToLocal(const Point &z) const { return globalToLocal(z); }

#ifdef BUILD_IA

IAPoint Tetrahedron::FaceLocalToGlobal(int face, const IAPoint &local) const {
  return faceLocalToGlobal(face, local);
}

IAPoint Tetrahedron::FaceLocalToLocal(int face, const IAPoint &local) const {
  return faceLocalToLocal(face, local);
}

IAPoint Tetrahedron::FaceLocalToLocalOriented(int faceId, const IAPoint &localOnFace) const {
  return faceLocalToLocalOriented(faceId, localOnFace);
}

IAPoint Tetrahedron::LocalToGlobal(const IAPoint &z) const { return localToGlobal(z); }

IAPoint Tetrahedron::GlobalToLocal(const IAPoint &z) const { return globalToLocal(z); }

IATransformation Tetrahedron::GetTransformation(const IAPoint &z) const {
  return getTransformation(z);
}

const IAPoint *Tetrahedron::LocalFaceNormalsIATet() {
  static const IAPoint *localFaceNormalsIA = ReferenceFaceNormals<IAPoint>(*ReferenceTetrahedron());
  return localFaceNormalsIA;
}

IAInterval Tetrahedron::LocalFaceAreaIA(int faceId) const {
  switch (faceId) {
  case 0:
  case 2:
  case 3:
    return IAInterval::c2r();
  case 1:
    return IAInterval::Sqrt3d2();
  default:
    THROW("faceId " + std::to_string(faceId) + " not possible in Tetrahedron")
  }
}

const IAPoint &Tetrahedron::LocalFaceNormalIA(int faceId) const {
  return LocalFaceNormalsIATet()[faceId];
}

#endif
