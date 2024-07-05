#include "Hexahedron.hpp"

// ----------------------------------------------------
//    HEXAHEDRON
// ----------------------------------------------------
Hexahedron::Hexahedron(const vector<Point> &z, int sd) : Cell(z, sd){};

CELLTYPE Hexahedron::Type() const { return HEXAHEDRON; }

CELLTYPE Hexahedron::ReferenceType() const { return HEXAHEDRON; }

int Hexahedron::Corners() const { return 8; }

const Point &Hexahedron::LocalCorner(int i) const { return LocalCornersHex()[i]; }

const Point &Hexahedron::LocalEdge(int i) const { return LocalEdgesHex()[i]; }

const Point &Hexahedron::LocalFace(int i) const { return LocalFacesHex()[i]; }

const Point &Hexahedron::LocalFaceNormal(int i) const { return LocalFaceNormalsHex()[i]; }

double Hexahedron::LocalFaceArea(int i) const { return 1; }

double Hexahedron::FaceArea(int i) const {
  return norm((FaceCorner(i, 1) - FaceCorner(i, 0)) ^ (FaceCorner(i, 3) - FaceCorner(i, 0)));
}

const Point &Hexahedron::LocalCenter() const { return *LocalCenterHex(); }

Point Hexahedron::FaceLocalToGlobal(int face, const Point &local) const {
  return faceLocalToGlobal(face, local);
}

Point Hexahedron::FaceLocalToLocal(int face, const Point &local) const {
  return faceLocalToLocal(face, local);
}

Point Hexahedron::FaceLocalToLocalOriented(int faceId, const Point &localOnFace) const {
  return faceLocalToLocalOriented(faceId, localOnFace);
}

Point Hexahedron::LocalToGlobal(const Point &z) const { return localToGlobal(z); }

Transformation Hexahedron::GetTransformation(const Point &z) const { return getTransformation(z); }

#ifdef BUILD_IA

IAPoint Hexahedron::FaceLocalToGlobal(int face, const IAPoint &local) const {
  return faceLocalToGlobal(face, local);
}

IAPoint Hexahedron::FaceLocalToLocal(int face, const IAPoint &local) const {
  return faceLocalToLocal(face, local);
}

IAPoint Hexahedron::FaceLocalToLocalOriented(int faceId, const IAPoint &localOnFace) const {
  return faceLocalToLocalOriented(faceId, localOnFace);
}

IAPoint Hexahedron::LocalToGlobal(const IAPoint &z) const { return localToGlobal(z); }

IAPoint Hexahedron::GlobalToLocal(const IAPoint &z) const { return globalToLocal(z); }

IATransformation Hexahedron::GetTransformation(const IAPoint &z) const {
  return getTransformation(z);
}

IAInterval Hexahedron::LocalFaceAreaIA(int i) const { return IAInterval(1.0); }

const IAPoint *Hexahedron::LocalFaceNormalsIAHex() {
  static const IAPoint *localFaceNormalsIA = ReferenceFaceNormals<IAPoint>(*ReferenceHexahedron());
  return localFaceNormalsIA;
}

const IAPoint &Hexahedron::LocalFaceNormalIA(int faceId) const {
  return LocalFaceNormalsIAHex()[faceId];
}

#endif

int Hexahedron::Edges() const { return 12; }

Point Hexahedron::Edge(int i) const {
  switch (i) {
  case 0:
    return 0.5 * (Corner(0) + Corner(1));
  case 1:
    return 0.5 * (Corner(1) + Corner(2));
  case 2:
    return 0.5 * (Corner(2) + Corner(3));
  case 3:
    return 0.5 * (Corner(3) + Corner(0));
  case 4:
    return 0.5 * (Corner(0) + Corner(4));
  case 5:
    return 0.5 * (Corner(1) + Corner(5));
  case 6:
    return 0.5 * (Corner(2) + Corner(6));
  case 7:
    return 0.5 * (Corner(3) + Corner(7));
  case 8:
    return 0.5 * (Corner(4) + Corner(5));
  case 9:
    return 0.5 * (Corner(5) + Corner(6));
  case 10:
    return 0.5 * (Corner(6) + Corner(7));
  case 11:
    return 0.5 * (Corner(7) + Corner(4));
  }
  THROW("Not Implemented for i > 11")
}

short Hexahedron::edgecorner(unsigned short i, unsigned short j) const {
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
      return 2;
    case 1:
      return 3;
    }
  case 3:
    switch (j) {
    case 0:
      return 3;
    case 1:
      return 0;
    }
  case 4:
    switch (j) {
    case 0:
      return 0;
    case 1:
      return 4;
    }
  case 5:
    switch (j) {
    case 0:
      return 1;
    case 1:
      return 5;
    }
  case 6:
    switch (j) {
    case 0:
      return 2;
    case 1:
      return 6;
    }
  case 7:
    switch (j) {
    case 0:
      return 3;
    case 1:
      return 7;
    }
  case 8:
    switch (j) {
    case 0:
      return 4;
    case 1:
      return 5;
    }
  case 9:
    switch (j) {
    case 0:
      return 5;
    case 1:
      return 6;
    }
  case 10:
    switch (j) {
    case 0:
      return 6;
    case 1:
      return 7;
    }
  case 11:
    switch (j) {
    case 0:
      return 7;
    case 1:
      return 4;
    }
  }
  THROW("edgecorner (" + to_string(i) + ", " + to_string(j) + ") does not exist")
}

int Hexahedron::Faces() const { return 6; }

Point Hexahedron::Face(int i) const {
  switch (i) {
  case 0:
    return 0.25 * (Corner(0) + Corner(1) + Corner(2) + Corner(3));
  case 1:
    return 0.25 * (Corner(0) + Corner(1) + Corner(5) + Corner(4));
  case 2:
    return 0.25 * (Corner(1) + Corner(2) + Corner(5) + Corner(6));
  case 3:
    return 0.25 * (Corner(2) + Corner(3) + Corner(6) + Corner(7));
  case 4:
    return 0.25 * (Corner(0) + Corner(3) + Corner(4) + Corner(7));
  case 5:
    return 0.25 * (Corner(4) + Corner(5) + Corner(6) + Corner(7));
  }
  THROW("Not Implemented")
}

CELLTYPE Hexahedron::FaceType(int face) const { return QUADRILATERAL; }

short Hexahedron::FarFace(int i) const {
  if (i < Faces()) return 5 - i + (1 - 2 * i % 2);
  else THROW("Not Implemented")
  switch (i) {
  case 0:
    return 5;
  case 1:
    return 3;
  case 2:
    return 4;
  case 3:
    return 1;
  case 4:
    return 2;
  case 5:
    return 0;
  }
}

short Hexahedron::MiddleEdge(int i, int j) const {
  if (i == 5) i = 0;
  else if (i == 3) i = 1;
  else if (i == 4) i = 2;
  switch (i) {
  case 0:
    switch (j) {
    case 0:
      return 4;
    case 1:
      return 5;
    case 2:
      return 6;
    case 3:
      return 7;
    }
  case 1:
    switch (j) {
    case 0:
      return 1;
    case 1:
      return 3;
    case 2:
      return 9;
    case 3:
      return 11;
    }
  case 2:
    switch (j) {
    case 0:
      return 0;
    case 1:
      return 2;
    case 2:
      return 8;
    case 3:
      return 10;
    }
  }
  THROW("Not Implemented")
}

short Hexahedron::FaceCorners(int i) const { return 4; }

short Hexahedron::facecorner(unsigned short i, unsigned short j) const {
  switch (i) {
  case 0:
    switch (j) {
    case 0:
      return 0;
    case 1:
      return 3;
    case 2:
      return 2;
    case 3:
      return 1;
    }
  case 1:
    switch (j) {
    case 0:
      return 0;
    case 1:
      return 1;
    case 2:
      return 5;
    case 3:
      return 4;
    }
  case 2:
    switch (j) {
    case 0:
      return 1;
    case 1:
      return 2;
    case 2:
      return 6;
    case 3:
      return 5;
    }
  case 3:
    switch (j) {
    case 0:
      return 2;
    case 1:
      return 3;
    case 2:
      return 7;
    case 3:
      return 6;
    }
  case 4:
    switch (j) {
    case 0:
      return 3;
    case 1:
      return 0;
    case 2:
      return 4;
    case 3:
      return 7;
    }
  case 5:
    switch (j) {
    case 0:
      return 4;
    case 1:
      return 5;
    case 2:
      return 6;
    case 3:
      return 7;
    }
  }
  THROW("facecorner (" + to_string(i) + ", " + to_string(j) + ") does not exist")
}

short Hexahedron::FaceEdges(int i) const { return 4; }

bool Hexahedron::faceedgecircuit(int i) const {
  switch (i) {
  case 0:
    return true;
  case 1:
    return false;
  case 2:
    return false;
  case 3:
    return false;
  case 4:
    return false;
  case 5:
    return true;
  }
  THROW("Not Implemented")
}

short Hexahedron::faceedge(int i, int j) const {
  switch (i) {
  case 0:
    switch (j) {
    case 0:
      return 0;
    case 1:
      return 1;
    case 2:
      return 2;
    case 3:
      return 3;
    }
  case 1:
    switch (j) {
    case 0:
      return 0;
    case 1:
      return 4;
    case 2:
      return 5;
    case 3:
      return 8;
    }
  case 2:
    switch (j) {
    case 0:
      return 1;
    case 1:
      return 5;
    case 2:
      return 6;
    case 3:
      return 9;
    }
  case 3:
    switch (j) {
    case 0:
      return 2;
    case 1:
      return 6;
    case 2:
      return 7;
    case 3:
      return 10;
    }
  case 4:
    switch (j) {
    case 0:
      return 3;
    case 1:
      return 7;
    case 2:
      return 4;
    case 3:
      return 11;
    }
  case 5:
    switch (j) {
    case 0:
      return 8;
    case 1:
      return 9;
    case 2:
      return 10;
    case 3:
      return 11;
    }
  }
  THROW("faceedge (" + to_string(i) + ", " + to_string(j) + ") does not exist")
}

void Hexahedron::GetRefinePoints(std::vector<Point> &z) const {
  z.resize(27);
  int n = 0;
  for (int i = 0; i < Corners(); ++i)
    z[n++] = Corner(i);
  for (int i = 0; i < Edges(); ++i)
    z[n++] = Edge(i);
  for (int i = 0; i < Faces(); ++i)
    z[n++] = Face(i);
  z[n++] = Center();
}

const vector<Rule> &Hexahedron::Refine(vector<Point> &z) const {
  GetRefinePoints(z);
  return *Rhex;
}

void Hexahedron::GetRefineBarycentricPoints(std::vector<Point> &z, const Point &shiftCenter) const {
  z.resize(15);
  for (int i = 0; i < 8; ++i)
    z[i] = Corner(i);
  for (int i = 0, n = 8; i < 6; ++i, ++n)
    z[n] = Face(i);
  z[14] = Center() + shiftCenter;
}

const std::vector<Rule> &Hexahedron::RefineBarycentric(std::vector<Point> &z,
                                                       const Point &shiftCenter) const {
  GetRefineBarycentricPoints(z, shiftCenter);
  return *Rhex_barycentric;
}

int Hexahedron::dim() const { return 3; }

bool Hexahedron::plane() const { return false; }

const Hexahedron *Hexahedron::ReferenceHexahedron() {
  static const Hexahedron *refHexa = new Hexahedron(Points(8, Hexahedron::LocalCornersHex()), -1);
  return refHexa;
}

const Cell *Hexahedron::ReferenceCell() const { return ReferenceHexahedron(); }

int Hexahedron::Children() const { return 8; }

Point Hexahedron::Child(int i) const {
  switch (i) {
  case 0:
    return 0.015625
           * (Corner(0) + 3 * Corner(1) + 9 * Corner(2) + 3 * Corner(3) + 3 * Corner(4)
              + 9 * Corner(5) + 27 * Corner(6) + 9 * Corner(7));
  case 1:
    return 0.015625
           * (3 * Corner(0) + Corner(1) + 3 * Corner(2) + 9 * Corner(3) + 9 * Corner(4)
              + 3 * Corner(5) + 9 * Corner(6) + 27 * Corner(7));
  case 2:
    return 0.015625
           * (9 * Corner(0) + 3 * Corner(1) + Corner(2) + 3 * Corner(3) + 27 * Corner(4)
              + 9 * Corner(5) + 3 * Corner(6) + 9 * Corner(7));
  case 3:
    return 0.015625
           * (3 * Corner(0) + 9 * Corner(1) + 3 * Corner(2) + Corner(3) + 9 * Corner(4)
              + 27 * Corner(5) + 9 * Corner(6) + 3 * Corner(7));
  case 4:
    return 0.015625
           * (3 * Corner(0) + 9 * Corner(1) + 27 * Corner(2) + 9 * Corner(3) + Corner(4)
              + 3 * Corner(5) + 9 * Corner(6) + 3 * Corner(7));
  case 5:
    return 0.015625
           * (9 * Corner(0) + 3 * Corner(1) + 9 * Corner(2) + 27 * Corner(3) + 3 * Corner(4)
              + Corner(5) + 3 * Corner(6) + 9 * Corner(7));
  case 6:
    return 0.015625
           * (27 * Corner(0) + 9 * Corner(1) + 3 * Corner(2) + 9 * Corner(3) + 9 * Corner(4)
              + 3 * Corner(5) + Corner(6) + 3 * Corner(7));
  case 7:
    return 0.015625
           * (9 * Corner(0) + 27 * Corner(1) + 9 * Corner(2) + 3 * Corner(3) + 3 * Corner(4)
              + 9 * Corner(5) + 3 * Corner(6) + Corner(7));
  }
  Exit("Not implemented") return Origin;
}

const Point *Hexahedron::LocalCornersHex() {
  static const Point localCorners[] = {Point(0.0, 0.0, 0.0), Point(1.0, 0.0, 0.0),
                                       Point(1.0, 1.0, 0.0), Point(0.0, 1.0, 0.0),
                                       Point(0.0, 0.0, 1.0), Point(1.0, 0.0, 1.0),
                                       Point(1.0, 1.0, 1.0), Point(0.0, 1.0, 1.0)};

  return localCorners;
}

const Point *Hexahedron::LocalEdgesHex() {
  static const Point *localEdges = ReferenceEdges(*ReferenceHexahedron());
  return localEdges;
}

const Point *Hexahedron::LocalFacesHex() {
  static const Point *localFaces = ReferenceFaces(*ReferenceHexahedron());
  return localFaces;
}

const Point *Hexahedron::LocalFaceNormalsHex() {
  static const Point *localFaceNormals = ReferenceFaceNormals<Point>(*ReferenceHexahedron());
  return localFaceNormals;
}

const Point *Hexahedron::LocalCenterHex() {
  static const Point *localCenter = ReferenceCenter(*ReferenceHexahedron());
  return localCenter;
}

Point Hexahedron::GlobalToLocal(const Point &z) const { return globalToLocal(z); }

Hexahedron20::Hexahedron20(const vector<Point> &z, int sd) : Hexahedron(z, sd) {
  int n = 8;
  for (int i = 0; i < Edges(); ++i)
    corners[n++] = Edge(i);
}

Point Hexahedron20::LocalToGlobal(const Point &z) const {
  return ((1.0 - z[0]) * (1.0 - z[1]) * (1.0 - z[2]) * (1.0 - 2.0 * z[0] - 2.0 * z[1] - 2.0 * z[2]))
             * corners[0]
         + ((z[0]) * (1.0 - z[1]) * (1.0 - z[2])
            * (2.0 * (z[0]) + 2.0 * (1.0 - z[1]) + 2.0 * (1.0 - z[2]) - 5.0))
               * corners[1]
         + ((z[0]) * (z[1]) * (1.0 - z[2])
            * (2.0 * (z[0]) + 2.0 * (z[1]) + 2.0 * (1.0 - z[2]) - 5.0))
               * corners[2]
         + ((1.0 - z[0]) * (z[1]) * (1.0 - z[2])
            * (2.0 * (1.0 - z[0]) + 2.0 * (z[1]) + 2.0 * (1.0 - z[2]) - 5.0))
               * corners[3]
         + ((1.0 - z[0]) * (1.0 - z[1]) * (z[2])
            * (2.0 * (1.0 - z[0]) + 2.0 * (1.0 - z[1]) + 2.0 * (z[2]) - 5.0))
               * corners[4]
         + ((z[0]) * (1.0 - z[1]) * (z[2])
            * (2.0 * (z[0]) + 2.0 * (1.0 - z[1]) + 2.0 * (z[2]) - 5.0))
               * corners[5]
         + ((z[0]) * (z[1]) * (z[2]) * (2.0 * (z[0]) + 2.0 * (z[1]) + 2.0 * (z[2]) - 5.0))
               * corners[6]
         + ((1.0 - z[0]) * (z[1]) * (z[2])
            * (2.0 * (1.0 - z[0]) + 2.0 * (z[1]) + 2.0 * (z[2]) - 5.0))
               * corners[7]
         + (4.0 * (1.0 - z[0]) * (1.0 - z[1]) * (1.0 - z[2]) * z[0]) * corners[8]
         + (4.0 * (z[0]) * (1.0 - z[1]) * (1.0 - z[2]) * z[1]) * corners[9]
         + (4.0 * (1.0 - z[0]) * (z[1]) * (1.0 - z[2]) * z[0]) * corners[10]
         + (4.0 * (1.0 - z[0]) * (1.0 - z[1]) * (1.0 - z[2]) * z[1]) * corners[11]
         + (4.0 * (1.0 - z[0]) * (1.0 - z[1]) * (1.0 - z[2]) * z[2]) * corners[12]
         + (4.0 * (z[0]) * (1.0 - z[1]) * (1.0 - z[2]) * z[2]) * corners[13]
         + (4.0 * (z[0]) * (z[1]) * (1.0 - z[2]) * z[2]) * corners[14]
         + (4.0 * (1.0 - z[0]) * (z[1]) * (1.0 - z[2]) * z[2]) * corners[15]
         + (4.0 * (1.0 - z[0]) * (1.0 - z[1]) * (z[2]) * z[0]) * corners[16]
         + (4.0 * (z[0]) * (1.0 - z[1]) * (z[2]) * z[1]) * corners[17]
         + (4.0 * (1.0 - z[0]) * (z[1]) * (z[2]) * z[0]) * corners[18]
         + (4.0 * (1.0 - z[0]) * (1.0 - z[1]) * (z[2]) * z[1]) * corners[19];
}

Transformation Hexahedron20::GetTransformation(const Point &z) const {
  vector<Point> y(3);
  for (int k = 0; k < y.size(); ++k) {
    y[k] =
        Point(-(1.0 - z[1]) * (1.0 - z[2]) * (1.0 - 2.0 * z[0] - 2.0 * z[1] - 2.0 * z[2])
                  - 2.0 * (1.0 - z[0]) * (1.0 - z[1]) * (1.0 - z[2]),
              -(1.0 - z[0]) * (1.0 - z[2]) * (1.0 - 2.0 * z[0] - 2.0 * z[1] - 2.0 * z[2])
                  - 2.0 * (1.0 - z[0]) * (1.0 - z[1]) * (1.0 - z[2]),
              -(1.0 - z[0]) * (1.0 - z[1]) * (1.0 - 2.0 * z[0] - 2.0 * z[1] - 2.0 * z[2])
                  - 2.0 * (1.0 - z[0]) * (1.0 - z[1]) * (1.0 - z[2]))
            * corners[0][k]
        + Point((1.0 - z[1]) * (1.0 - z[2])
                        * (2.0 * (z[0]) + 2.0 * (1.0 - z[1]) + 2.0 * (1.0 - z[2]) - 5.0)
                    + 2.0 * (z[0]) * (1.0 - z[1]) * (1.0 - z[2]),
                -z[0] * (1.0 - z[2])
                        * (2.0 * (z[0]) + 2.0 * (1.0 - z[1]) + 2.0 * (1.0 - z[2]) - 5.0)
                    - 2.0 * (z[0]) * (1.0 - z[1]) * (1.0 - z[2]),
                -z[0] * (1.0 - z[1])
                        * (2.0 * (z[0]) + 2.0 * (1.0 - z[1]) + 2.0 * (1.0 - z[2]) - 5.0)
                    - 2.0 * (z[0]) * (1.0 - z[1]) * (1.0 - z[2]))
              * corners[1][k]
        + Point(z[1] * (1.0 - z[2]) * (2.0 * (z[0]) + 2.0 * (z[1]) + 2.0 * (1.0 - z[2]) - 5.0)
                    + 2.0 * z[0] * z[1] * (1.0 - z[2]),
                z[0] * (1.0 - z[2]) * (2.0 * (z[0]) + 2.0 * (z[1]) + 2.0 * (1.0 - z[2]) - 5.0)
                    + 2.0 * z[0] * z[1] * (1.0 - z[2]),
                -z[0] * z[1] * (2.0 * (z[0]) + 2.0 * (z[1]) + 2.0 * (1.0 - z[2]) - 5.0)
                    - 2.0 * z[0] * z[1] * (1.0 - z[2]))
              * corners[2][k]
        + Point(-z[1] * (1.0 - z[2])
                        * (2.0 * (1.0 - z[0]) + 2.0 * (z[1]) + 2.0 * (1.0 - z[2]) - 5.0)
                    - 2.0 * (1.0 - z[0]) * z[1] * (1.0 - z[2]),
                (1.0 - z[0]) * (1.0 - z[2])
                        * (2.0 * (1.0 - z[0]) + 2.0 * (z[1]) + 2.0 * (1.0 - z[2]) - 5.0)
                    + 2.0 * (1.0 - z[0]) * z[1] * (1.0 - z[2]),
                -(1.0 - z[0]) * z[1]
                        * (2.0 * (1.0 - z[0]) + 2.0 * (z[1]) + 2.0 * (1.0 - z[2]) - 5.0)
                    - 2.0 * (1.0 - z[0]) * z[1] * (1.0 - z[2]))
              * corners[3][k]
        + Point(-(1.0 - z[1]) * z[2] * (2.0 * (1.0 - z[0]) + 2.0 * (1.0 - z[1]) + 2.0 * z[2] - 5.0)
                    - 2.0 * (1.0 - z[0]) * (1.0 - z[1]) * z[2],
                -(1.0 - z[0]) * z[2] * (2.0 * (1.0 - z[0]) + 2.0 * (1.0 - z[1]) + 2.0 * z[2] - 5.0)
                    - 2.0 * (1.0 - z[0]) * (1.0 - z[1]) * z[2],
                (1.0 - z[0]) * (1.0 - z[1])
                        * (2.0 * (1.0 - z[0]) + 2.0 * (1.0 - z[1]) + 2.0 * z[2] - 5.0)
                    + 2.0 * (1.0 - z[0]) * (1.0 - z[1]) * z[2])
              * corners[4][k]
        + Point((1.0 - z[1]) * z[2] * (2.0 * z[0] + 2.0 * (1.0 - z[1]) + 2.0 * z[2] - 5.0)
                    + 2.0 * z[0] * (1.0 - z[1]) * z[2],
                -z[0] * z[2] * (2.0 * z[0] + 2.0 * (1.0 - z[1]) + 2.0 * z[2] - 5.0)
                    - 2.0 * z[0] * (1.0 - z[1]) * z[2],
                z[0] * (1.0 - z[1]) * (2.0 * z[0] + 2.0 * (1.0 - z[1]) + 2.0 * z[2] - 5.0)
                    + 2.0 * z[0] * (1.0 - z[1]) * z[2])
              * corners[5][k]
        + Point(z[1] * z[2] * (2.0 * z[0] + 2.0 * z[1] + 2.0 * z[2] - 5.0)
                    + 2.0 * z[0] * z[1] * z[2],
                z[0] * z[2] * (2.0 * z[0] + 2.0 * z[1] + 2.0 * z[2] - 5.0)
                    + 2.0 * z[0] * z[1] * z[2],
                z[0] * z[1] * (2.0 * (z[0]) + 2.0 * (z[1]) + 2.0 * z[2] - 5.0)
                    + 2.0 * z[0] * z[1] * z[2])
              * corners[6][k]
        + Point(-z[1] * z[2] * (2.0 * (1.0 - z[0]) + 2.0 * z[1] + 2.0 * z[2] - 5.0)
                    - 2.0 * (1.0 - z[0]) * z[1] * z[2],
                (1.0 - z[0]) * z[2] * (2.0 * (1.0 - z[0]) + 2.0 * z[1] + 2.0 * z[2] - 5.0)
                    + 2.0 * (1.0 - z[0]) * z[1] * z[2],
                (1.0 - z[0]) * z[1] * (2.0 * (1.0 - z[0]) + 2.0 * z[1] + 2.0 * z[2] - 5.0)
                    + 2.0 * (1.0 - z[0]) * z[1] * z[2])
              * corners[7][k]
        + Point(4.0 * (1.0 - 2.0 * z[0]) * (1.0 - z[1]) * (1.0 - z[2]),
                -4.0 * (1.0 - z[0]) * (1.0 - z[2]) * z[0],
                -4.0 * (1.0 - z[0]) * (1.0 - z[1]) * z[0])
              * corners[8][k]
        + Point(4.0 * (1.0 - z[1]) * (1.0 - z[2]) * z[1],
                4.0 * (z[0]) * (1.0 - 2.0 * z[1]) * (1.0 - z[2]),
                -4.0 * (z[0]) * (1.0 - z[1]) * z[1])
              * corners[9][k]
        + Point(4.0 * (1.0 - 2.0 * z[0]) * z[1] * (1.0 - z[2]),
                4.0 * (1.0 - z[0]) * (1.0 - z[2]) * z[0], -4.0 * (1.0 - z[0]) * z[1] * z[0])
              * corners[10][k]
        + Point(-4.0 * (1.0 - z[1]) * (1.0 - z[2]) * z[1],
                4.0 * (1.0 - z[0]) * (1.0 - 2.0 * z[1]) * (1.0 - z[2]),
                -4.0 * (1.0 - z[0]) * (1.0 - z[1]) * z[1])
              * corners[11][k]
        + Point(-4.0 * (1.0 - z[1]) * (1.0 - z[2]) * z[2],
                -4.0 * (1.0 - z[0]) * (1.0 - z[2]) * z[2],
                4.0 * (1.0 - z[0]) * (1.0 - z[1]) * (1.0 - 2.0 * z[2]))
              * corners[12][k]
        + Point(4.0 * (1.0 - z[1]) * (1.0 - z[2]) * z[2], -4.0 * (z[0]) * (1.0 - z[2]) * z[2],
                4.0 * z[0] * (1.0 - z[1]) * (1.0 - 2.0 * z[2]))
              * corners[13][k]
        + Point(4.0 * z[1] * (1.0 - z[2]) * z[2], 4.0 * z[0] * (1.0 - z[2]) * z[2],
                4.0 * z[0] * z[1] * (1.0 - 2.0 * z[2]))
              * corners[14][k]
        + Point(-4.0 * z[1] * (1.0 - z[2]) * z[2], 4.0 * (1.0 - z[0]) * (1.0 - z[2]) * z[2],
                4.0 * (1.0 - z[0]) * z[1] * (1.0 - 2.0 * z[2]))
              * corners[15][k]
        + Point(4.0 * (1.0 - 2.0 * z[0]) * (1.0 - z[1]) * z[2], -4.0 * (1.0 - z[0]) * z[2] * z[0],
                4.0 * (1.0 - z[0]) * (1.0 - z[1]) * z[0])
              * corners[16][k]
        + Point(4.0 * (1.0 - z[1]) * z[2] * z[1], 4.0 * z[0] * (1.0 - 2.0 * z[1]) * z[2],
                4.0 * z[0] * (1.0 - z[1]) * z[1])
              * corners[17][k]
        + Point(4.0 * (1.0 - 2.0 * z[0]) * z[1] * z[2], 4.0 * (1.0 - z[0]) * z[2] * z[0],
                4.0 * (1.0 - z[0]) * z[1] * z[0])
              * corners[18][k]
        + Point(-4.0 * (1.0 - z[1]) * z[2] * z[1], 4.0 * (1.0 - z[0]) * (1.0 - 2.0 * z[1]) * z[2],
                4.0 * (1.0 - z[0]) * (1.0 - z[1]) * z[1])
              * corners[19][k];
  }
  return Transformation(y);
}

CELLTYPE Hexahedron20::Type() const { return HEXAHEDRON20; }

void Hexahedron20::GetRefinePoints(std::vector<Point> &z) const {
  z.resize(81);
  int n = 0;
  for (int i = 0; i < Corners(); ++i) {
    z[n++] = Corner(i);
  }
  for (int i = 0; i < Edges(); ++i) {
    z[n++] = LocalToGlobal(LocalEdge(i));
  }
  for (int i = 0; i < Faces(); ++i) {
    z[n++] = LocalToGlobal(LocalFace(i));
  }
  z[n++] = LocalToGlobal(LocalCenter());
  for (int i = 0; i < Edges(); ++i) {
    z[n++] = LocalToGlobal(0.5 * (LocalCorner(edgecorner(i, 0)) + LocalEdge(i)));
    z[n++] = LocalToGlobal(0.5 * (LocalCorner(edgecorner(i, 1)) + LocalEdge(i)));
  }
  for (int i = 0; i < Faces(); ++i) {
    z[n++] = LocalToGlobal(0.5 * (LocalEdge(faceedge(i, 0)) + LocalFace(i)));
    z[n++] = LocalToGlobal(0.5 * (LocalEdge(faceedge(i, 1)) + LocalFace(i)));
    z[n++] = LocalToGlobal(0.5 * (LocalEdge(faceedge(i, 2)) + LocalFace(i)));
    z[n++] = LocalToGlobal(0.5 * (LocalEdge(faceedge(i, 3)) + LocalFace(i)));
    z[n++] = LocalToGlobal(0.5 * (LocalCenter() + LocalFace(i)));
  }
}

const vector<Rule> &Hexahedron20::Refine(vector<Point> &z) const {
  GetRefinePoints(z);
  return *R20hex;
}

const Hexahedron20 *Hexahedron20::ReferenceHexahedron20() {
  static const Hexahedron20 *refHexa20 =
      new Hexahedron20(Points(8, Hexahedron::LocalCornersHex(), 12, Hexahedron::LocalEdgesHex()),
                       -1);
  return refHexa20;
}

static const double c_27[27][3] = {{0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {1.0, 1.0, 0.0},
                                   {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}, {1.0, 0.0, 1.0},
                                   {1.0, 1.0, 1.0}, {0.0, 1.0, 1.0}, {0.5, 0.0, 0.0},
                                   {1.0, 0.5, 0.0}, {0.5, 1.0, 0.0}, {0.0, 0.5, 0.0},
                                   {0.0, 0.0, 0.5}, {1.0, 0.0, 0.5}, {1.0, 1.0, 0.5},
                                   {0.0, 1.0, 0.5}, {0.5, 0.0, 1.0}, {1.0, 0.5, 1.0},
                                   {0.5, 1.0, 1.0}, {0.0, 0.5, 1.0}, {0.5, 0.5, 0.0},
                                   {0.5, 0.0, 0.5}, {1.0, 0.5, 0.5}, {0.5, 1.0, 0.5},
                                   {0.0, 0.5, 0.5}, {0.5, 0.5, 1.0}, {0.5, 0.5, 0.5}};

inline double q01(const double x) { return x * (2 * x - 1); }

inline double q00(const double x) { return (x - 1) * (2 * x - 1); }

inline double q11(const double x) { return 4 * x * (1 - x); }

inline double Dq01(const double x) { return 4 * x - 1; }

inline double Dq00(const double x) { return 4 * x - 3; }

inline double Dq11(const double x) { return -8 * x + 4; }

inline double q1(const double x) { return x * (2 * x - 1); }

inline double q2(const double x) { return (x - 1) * (2 * x - 1); }

inline double q3(const double x) { return 4 * x * (1 - x); }

inline double Dq1(const double x) { return 4 * x - 1; }

inline double Dq2(const double x) { return 4 * x - 3; }

inline double Dq3(const double x) { return -8 * x + 4; }

inline double shapefct(int i, int k, const double x) {
  double s = c_27[i][k];
  if (s == 1.0) return q1(x);
  if (s == 0.0) return q2(x);
  if (s == 0.5) return q3(x);
  THROW("Not Implemented")
}

inline double Dshapefct(int i, int k, const double x) {
  double s = c_27[i][k];
  if (s == 1.0) return Dq1(x);
  if (s == 0.0) return Dq2(x);
  if (s == 0.5) return Dq3(x);
  THROW("Not Implemented")
}

Hexahedron27::Hexahedron27(const vector<Point> &z, int sd) : Hexahedron(z, sd) {
  int n = 8;
  for (int i = 0; i < Edges(); ++i)
    corners[n++] = Edge(i);
  for (int i = 0; i < Faces(); ++i)
    corners[n++] = Face(i);
  corners[n++] = Center();
}

Point Hexahedron27::LocalToGlobal(const Point &z) const {
  return q00(z[0]) * q00(z[1]) * q00(z[2]) * corners[0]
         + q01(z[0]) * q00(z[1]) * q00(z[2]) * corners[1]
         + q01(z[0]) * q01(z[1]) * q00(z[2]) * corners[2]
         + q00(z[0]) * q01(z[1]) * q00(z[2]) * corners[3]
         + q00(z[0]) * q00(z[1]) * q01(z[2]) * corners[4]
         + q01(z[0]) * q00(z[1]) * q01(z[2]) * corners[5]
         + q01(z[0]) * q01(z[1]) * q01(z[2]) * corners[6]
         + q00(z[0]) * q01(z[1]) * q01(z[2]) * corners[7]
         + q11(z[0]) * q00(z[1]) * q00(z[2]) * corners[8]
         + q01(z[0]) * q11(z[1]) * q00(z[2]) * corners[9]
         + q11(z[0]) * q01(z[1]) * q00(z[2]) * corners[10]
         + q00(z[0]) * q11(z[1]) * q00(z[2]) * corners[11]
         + q00(z[0]) * q00(z[1]) * q11(z[2]) * corners[12]
         + q01(z[0]) * q00(z[1]) * q11(z[2]) * corners[13]
         + q01(z[0]) * q01(z[1]) * q11(z[2]) * corners[14]
         + q00(z[0]) * q01(z[1]) * q11(z[2]) * corners[15]
         + q11(z[0]) * q00(z[1]) * q01(z[2]) * corners[16]
         + q01(z[0]) * q11(z[1]) * q01(z[2]) * corners[17]
         + q11(z[0]) * q01(z[1]) * q01(z[2]) * corners[18]
         + q00(z[0]) * q11(z[1]) * q01(z[2]) * corners[19]
         + q11(z[0]) * q11(z[1]) * q00(z[2]) * corners[20]
         + q11(z[0]) * q00(z[1]) * q11(z[2]) * corners[21]
         + q01(z[0]) * q11(z[1]) * q11(z[2]) * corners[22]
         + q11(z[0]) * q01(z[1]) * q11(z[2]) * corners[23]
         + q00(z[0]) * q11(z[1]) * q11(z[2]) * corners[24]
         + q11(z[0]) * q11(z[1]) * q01(z[2]) * corners[25]
         + q11(z[0]) * q11(z[1]) * q11(z[2]) * corners[26];
}

Transformation Hexahedron27::GetTransformation(const Point &z) const {
  vector<Point> y(3);
  for (int k = 0; k < y.size(); ++k) {
    y[k] = Point();
    for (int i = 0; i < 27; ++i) {
      y[k] += Point(Dshapefct(i, 0, z[0]) * shapefct(i, 1, z[1]) * shapefct(i, 2, z[2]),
                    shapefct(i, 0, z[0]) * Dshapefct(i, 1, z[1]) * shapefct(i, 2, z[2]),
                    shapefct(i, 0, z[0]) * shapefct(i, 1, z[1]) * Dshapefct(i, 2, z[2]))
              * corners[i][k];
    }
  }
  return Transformation(y);
}

CELLTYPE Hexahedron27::Type() const { return HEXAHEDRON27; }

void Hexahedron27::GetRefinePoints(std::vector<Point> &z) const {
  THROW("Refinement not correct. See also end of Cell.C")
  z.resize(125);
  int n = 0;
  for (int i = 0; i < Corners(); ++i) {
    z[n++] = Corner(i);
  }
  for (int i = 0; i < Edges(); ++i) {
    z[n++] = LocalToGlobal(LocalEdge(i));
  }
  for (int i = 0; i < Faces(); ++i) {
    z[n++] = LocalToGlobal(LocalFace(i));
  }
  z[n++] = LocalToGlobal(LocalCenter());
  for (int i = 0; i < Edges(); ++i) {
    z[n++] = LocalToGlobal(0.5 * (LocalCorner(edgecorner(i, 0)) + LocalEdge(i)));
    z[n++] = LocalToGlobal(0.5 * (LocalCorner(edgecorner(i, 1)) + LocalEdge(i)));
  }
  for (int i = 0; i < Faces(); ++i) {
    z[n++] = LocalToGlobal(0.5 * (LocalEdge(faceedge(i, 0)) + LocalFace(i)));
    z[n++] = LocalToGlobal(0.5 * (LocalEdge(faceedge(i, 1)) + LocalFace(i)));
    z[n++] = LocalToGlobal(0.5 * (LocalEdge(faceedge(i, 2)) + LocalFace(i)));
    z[n++] = LocalToGlobal(0.5 * (LocalEdge(faceedge(i, 3)) + LocalFace(i)));
    z[n++] = LocalToGlobal(0.5 * (LocalCenter() + LocalFace(i)));
  }
  for (int i = 0; i < Faces(); ++i) {
    z[n++] = LocalToGlobal(0.5 * (LocalCorner(facecorner(i, 0)) + LocalFace(i)));
    z[n++] = LocalToGlobal(0.5 * (LocalCorner(facecorner(i, 1)) + LocalFace(i)));
    z[n++] = LocalToGlobal(0.5 * (LocalCorner(facecorner(i, 2)) + LocalFace(i)));
    z[n++] = LocalToGlobal(0.5 * (LocalCorner(facecorner(i, 3)) + LocalFace(i)));
  }
  for (int i = 0; i < Edges(); ++i)
    z[n++] = LocalToGlobal(0.5 * (LocalEdge(i) + LocalCenter()));
  for (int i = 0; i < Corners(); ++i)
    z[n++] = LocalToGlobal(0.5 * (LocalCorner(i) + LocalCenter()));
}

const vector<Rule> &Hexahedron27::Refine(vector<Point> &z) const {
  GetRefinePoints(z);
  return *R27hex;
}

const Hexahedron27 *Hexahedron27::ReferenceHexahedron27() {
  static const Hexahedron27 *refHexa27 =
      new Hexahedron27(Points(8, Hexahedron::LocalCornersHex(), 12, Hexahedron::LocalEdgesHex(), 6,
                              Hexahedron::LocalFacesHex(), 1, Hexahedron::LocalCenterHex()),
                       -1);
  return refHexa27;
}
