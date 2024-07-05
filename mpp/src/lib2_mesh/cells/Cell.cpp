#include "Cell.hpp"
#include "Face.hpp"
#include "Mesh.hpp"

#include "Intval.hpp"
#include "Quadrilateral.hpp"
#include "Triangle.hpp"

#include "Buffer.hpp"


#if SpaceDimension >= 3
#include "Hexahedron.hpp"
#include "Tetrahedron.hpp"
#endif

#ifdef USE_SPACETIME
#include "TCell.hpp"
#endif


bool checkQuadType(const Cell &C) {
  for (int i = 0; i < C.Edges(); ++i)
    if (C[4 + i] != C.Edge(i)) return true;
  return false;
}

Cell *CreateCell(CELLTYPE tp, int subdomain, const vector<Point> &z, double a, double b) {
  if (tp == INTERVAL) return new Intval(z, subdomain);
#if SpaceDimension >= 2
  if (tp == TRIANGLE) return new Triangle(z, subdomain);
  if (tp == QUADRILATERAL) return new Quadrilateral(z, subdomain);
  if (tp == QUADRILATERAL2) {
    Cell *c = new Quadrilateral2(z, subdomain);
    if (checkQuadType(*c)) return c;
    // If Quad is not actually a Quadrilateral2
    delete c;
    vector<Point> zNew(z.begin(), z.begin() + 4);
    c = new Quadrilateral(zNew, subdomain);
    return c;
  }
#endif
#if SpaceDimension >= 3
  if (tp == TETRAHEDRON) return new Tetrahedron(z, subdomain);
  if (tp == HEXAHEDRON) return new Hexahedron(z, subdomain);
  if (tp == HEXAHEDRON20) return new Hexahedron20(z, subdomain);
  if (tp == HEXAHEDRON27) return new Hexahedron27(z, subdomain);
#endif
#ifdef USE_SPACETIME

  if (tp == SPACETIME_INTERVAL) {
    Cell *c = CreateCell(INTERVAL, subdomain, z);
    TimeInterval I(a, b);
    return new TCell(c, I);
  }
#if SpaceDimension >= 2
  if (tp == SPACETIME_QUADRILATERAL) {
    Cell *c = CreateCell(QUADRILATERAL, subdomain, z);
    TimeInterval I(a, b);
    return new TCell(c, I);
  }
  if (tp == SPACETIME_TRIANGLE) {
    Cell *c = CreateCell(TRIANGLE, subdomain, z);
    TimeInterval I(a, b);
    return new TCell(c, I);
  }
#endif
#if SpaceDimension >= 3
  if (tp == SPACETIME_HEXAHEDRON) {
    Cell *c = CreateCell(HEXAHEDRON, subdomain, z);
    TimeInterval I(a, b);
    return new TCell(c, I);
  }
#endif
#endif
  Exit(std::to_string(tp) + " not implemented; CreateCell in Cell.cpp")
}

const Cell *ReferenceCell(CELLTYPE tp) {
  switch (tp) {
  case INTERVAL:
    return Intval::ReferenceIntval();
#if SpaceDimension >= 2
  case TRIANGLE:
    return Triangle::ReferenceTriangle();
  case QUADRILATERAL:
    return Quadrilateral::ReferenceQuadrilateral();
#endif
#if SpaceDimension >= 3
  case TETRAHEDRON:
    return Tetrahedron::ReferenceTetrahedron();
  case HEXAHEDRON:
    return Hexahedron::ReferenceHexahedron();
#endif
#ifdef USE_SPACETIME
  case SPACETIME_INTERVAL:
    return ReferenceCell(SpaceCellType(SPACETIME_INTERVAL));
#if SpaceDimension >= 2
  case SPACETIME_QUADRILATERAL:
    return ReferenceCell(SpaceCellType(SPACETIME_QUADRILATERAL));
#endif
#endif
  default:
    THROW("Not implemented")
  }
}

Buffer &operator<<(Buffer &b, const Cell &c) {
  b << int(c.Type()) << c.Subdomain() << short(c.SpaceCell().Size());
  if (c.IsSpaceTime()) { b << c.min() << c.max(); }
  for (int i = 0; i < c.SpaceCell().Size(); ++i)
    b << c.SpaceCell().Corner(i);
#ifdef USE_DATAMESH
  b << c.GetData();
#endif
  return b;
}

Buffer &operator<<(Buffer &b, const cell &c) { return b << *c; }

Buffer &operator>>(Buffer &b, Cell *&c) {
  Assert(c == nullptr);
  int tp;
  short sd, n;
  b >> tp >> sd >> n;
  auto type = CELLTYPE(tp);
  double min = -1, max = -1;
#ifdef USE_SPACETIME
  if (isSpaceTimeCellType(type)) { b >> min >> max; }
#endif
  vector<Point> x(n);
  for (int i = 0; i < n; ++i)
    b >> x[i];
  c = CreateCell(type, sd, x, min, max);

#ifdef USE_DATAMESH
  DataContainer data;
  b >> data;
  c->SetData(data);
#endif
  return b;
}

bool LessCell_x(const cell &c0, const cell &c1) { return LessPoint_x(c0(), c1()); }

bool LessCell_y(const cell &c0, const cell &c1) { return LessPoint_y(c0(), c1()); }

bool LessCell_z(const cell &c0, const cell &c1) { return LessPoint_z(c0(), c1()); }

bool onFace(const Point &z, const Cell &D, int j) {
  for (int k = 0; k < D.FaceCorners(j); ++k)
    if (z == D.FaceCorner(j, k)) return true;
  return false;
}

bool onFace(const Cell &C, int i, const Cell &D, int j) {
  for (int k = 0; k < C.FaceCorners(i); ++k)
    for (int l = 0; l < D.FaceCorners(j); ++l)
      if (C.FaceCorner(i, k) == D.FaceCorner(j, l)) return true;
  return false;
}

bool onTriangleFace(const Cell &C, int i, const Cell &D, int j) {
  int n = 0;
  for (int k = 0; k < C.FaceEdges(i); ++k)
    for (int l = 0; l < D.FaceCorners(j); ++l)
      if (C.FaceEdge(i, k) == D.FaceCorner(j, l)) ++n;
  if (n == 3) return true;
  if (n < 2) return false;
  for (int k = 0; k < C.FaceCorners(i); ++k)
    for (int l = 0; l < D.FaceCorners(j); ++l)
      if (C.FaceCorner(i, k) == D.FaceCorner(j, l)) return true;
  return false;
}

void RefineFaces(const Cell &C, const Cell &c, vector<short> &f) {
  f.resize(c.Faces());
  for (int j = 0; j < c.Faces(); ++j)
    f[j] = -1;
  for (int i = 0; i < C.Faces(); ++i) {
    if (C.FaceCorners(i) == 3) {
      for (int j = 0; j < c.Faces(); ++j)
        if (onTriangleFace(C, i, c, j)) f[j] = short(i);
      continue;
    }
    Point z = C.Face(i);
    for (int j = 0; j < c.Faces(); ++j)
      if (onFace(z, c, j))
        if (onFace(C, i, c, j)) f[j] = short(i);
  }
}

const vector<Rule> *RefineCell(const Cell &C, vector<Rule> &rules) {
  vector<Point> z;
  vector<Point> x;
  C.GetRefinePoints(z);

  for (Rule &rule : rules) {
    rule(z, x);
    Cell *c = CreateCell(rule.type(), -1, x);
    RefineFaces(C, *c, rule.face);
    delete c;
  }
  return new vector<Rule>(rules);
}

const vector<Rule> *RefineBarycentricCell(const Cell &C, vector<Rule> &rules) {
  vector<Point> z;
  vector<Point> x;
  C.GetRefineBarycentricPoints(z, Origin);

  for (Rule &rule : rules) {
    rule(z, x);
    Cell *c = CreateCell(rule.type(), -1, x);
    rule.face = std::vector<short>(c->Faces(), -1);
    for (int f = 0; f < c->Faces(); ++f) {
      Point face = c->Face(f);
      for (int ff = 0; ff < C.Faces(); ++ff) {
        if (face == C.Face(ff)) {
          rule.face[f] = short(ff);
          break;
        }
      }
    }
    delete c;
  }
  return new vector<Rule>(rules);
}

vector<Rule> rulesInt = {Rule(INTERVAL, {0, 2}), Rule(INTERVAL, {2, 1})};

const vector<Rule> *Intval::Rint = RefineCell(*Intval::ReferenceIntval(), rulesInt);

#if SpaceDimension >= 2
vector<Rule> R_tri = {Rule(TRIANGLE, {0, 3, 5}), Rule(TRIANGLE, {3, 1, 4}),
                      Rule(TRIANGLE, {5, 3, 4}), Rule(TRIANGLE, {5, 4, 2})};

const vector<Rule> *Triangle::Rtri = RefineCell(*Triangle::ReferenceTriangle(), R_tri);

vector<Rule> R_tri_barycentric = {Rule(TRIANGLE, {0, 1, 3}), Rule(TRIANGLE, {1, 2, 3}),
                                  Rule(TRIANGLE, {2, 0, 3})};

const vector<Rule> *Triangle::Rtri_barycentric =
    RefineBarycentricCell(*Triangle::ReferenceTriangle(), R_tri_barycentric);

vector<Rule> R_quad = {Rule(QUADRILATERAL, {0, 4, 8, 7}), Rule(QUADRILATERAL, {1, 5, 8, 4}),
                       Rule(QUADRILATERAL, {2, 6, 8, 5}), Rule(QUADRILATERAL, {3, 7, 8, 6})};

const vector<Rule> *Quadrilateral::Rquad =
    RefineCell(*Quadrilateral::ReferenceQuadrilateral(), R_quad);

vector<Rule> R_quad_barycentric = {Rule(TRIANGLE, {0, 1, 4}), Rule(TRIANGLE, {1, 2, 4}),
                                   Rule(TRIANGLE, {2, 3, 4}), Rule(TRIANGLE, {3, 0, 4})};

const vector<Rule> *Quadrilateral::Rquad_barycentric =
    RefineBarycentricCell(*Quadrilateral::ReferenceQuadrilateral(), R_quad_barycentric);

vector<Rule> R2_quad = {Rule(QUADRILATERAL2, {0, 4, 8, 7, 9, 11, 20, 19}),
                        Rule(QUADRILATERAL2, {1, 5, 8, 4, 12, 14, 11, 10}),
                        Rule(QUADRILATERAL2, {2, 6, 8, 5, 15, 17, 14, 13}),
                        Rule(QUADRILATERAL2, {3, 7, 8, 6, 18, 20, 17, 16})};

const vector<Rule> *Quadrilateral2::R2quad =
    RefineCell(*Quadrilateral2::ReferenceQuadrilateral2(), R2_quad);
#endif
#if SpaceDimension >= 3

vector<Rule> R_tet = {Rule(TETRAHEDRON, {0, 4, 6, 7}), Rule(TETRAHEDRON, {1, 5, 4, 8}),
                      Rule(TETRAHEDRON, {2, 6, 5, 9}), Rule(TETRAHEDRON, {7, 8, 9, 3}),
                      Rule(TETRAHEDRON, {4, 5, 6, 7}), Rule(TETRAHEDRON, {4, 7, 8, 5}),
                      Rule(TETRAHEDRON, {5, 8, 9, 7}), Rule(TETRAHEDRON, {6, 9, 7, 5})};

vector<Rule> R_tet1 = {Rule(TETRAHEDRON, {0, 4, 6, 7}), Rule(TETRAHEDRON, {1, 5, 4, 8}),
                       Rule(TETRAHEDRON, {2, 6, 5, 9}), Rule(TETRAHEDRON, {7, 8, 9, 3}),
                       Rule(TETRAHEDRON, {7, 8, 4, 6}), Rule(TETRAHEDRON, {7, 6, 9, 8}),
                       Rule(TETRAHEDRON, {8, 9, 5, 6}), Rule(TETRAHEDRON, {4, 5, 6, 8})};

vector<Rule> R_tet2 = {Rule(TETRAHEDRON, {0, 4, 6, 7}), Rule(TETRAHEDRON, {1, 5, 4, 8}),
                       Rule(TETRAHEDRON, {2, 6, 5, 9}), Rule(TETRAHEDRON, {7, 8, 9, 3}),
                       Rule(TETRAHEDRON, {6, 9, 7, 4}), Rule(TETRAHEDRON, {6, 4, 5, 9}),
                       Rule(TETRAHEDRON, {9, 5, 8, 4}), Rule(TETRAHEDRON, {7, 8, 4, 9})};

vector<Rule> R_tet_barycentric = {Rule(TETRAHEDRON, {0, 1, 3, 4}), Rule(TETRAHEDRON, {1, 2, 3, 4}),
                                  Rule(TETRAHEDRON, {2, 0, 3, 4}), Rule(TETRAHEDRON, {0, 2, 1, 4})};

const vector<Rule> *Tetrahedron::Rtet = RefineCell(*Tetrahedron::ReferenceTetrahedron(), R_tet);

const vector<Rule> *Tetrahedron::Rtet_1 = RefineCell(*Tetrahedron::ReferenceTetrahedron(), R_tet1);

const vector<Rule> *Tetrahedron::Rtet_2 = RefineCell(*Tetrahedron::ReferenceTetrahedron(), R_tet2);

const vector<Rule> *Tetrahedron::Rtet_barycentric =
    RefineBarycentricCell(*Tetrahedron::ReferenceTetrahedron(), R_tet_barycentric);

vector<Rule> R_hex = {Rule(HEXAHEDRON, {0, 8, 20, 11, 12, 21, 26, 24}),
                      Rule(HEXAHEDRON, {8, 1, 9, 20, 21, 13, 22, 26}),
                      Rule(HEXAHEDRON, {11, 20, 10, 3, 24, 26, 23, 15}),
                      Rule(HEXAHEDRON, {20, 9, 2, 10, 26, 22, 14, 23}),
                      Rule(HEXAHEDRON, {12, 21, 26, 24, 4, 16, 25, 19}),
                      Rule(HEXAHEDRON, {21, 13, 22, 26, 16, 5, 17, 25}),
                      Rule(HEXAHEDRON, {24, 26, 23, 15, 19, 25, 18, 7}),
                      Rule(HEXAHEDRON, {26, 22, 14, 23, 25, 17, 6, 18})};

const vector<Rule> *Hexahedron::Rhex = RefineCell(*Hexahedron::ReferenceHexahedron(), R_hex);

vector<Rule> R_hex_barycentric =
    {Rule(TETRAHEDRON, {14, 8, 0, 1}),  Rule(TETRAHEDRON, {14, 8, 1, 2}),
     Rule(TETRAHEDRON, {14, 8, 2, 3}),  Rule(TETRAHEDRON, {14, 8, 3, 0}),
     Rule(TETRAHEDRON, {14, 9, 0, 1}),  Rule(TETRAHEDRON, {14, 9, 1, 5}),
     Rule(TETRAHEDRON, {14, 9, 5, 4}),  Rule(TETRAHEDRON, {14, 9, 4, 0}),
     Rule(TETRAHEDRON, {14, 10, 1, 2}), Rule(TETRAHEDRON, {14, 10, 2, 6}),
     Rule(TETRAHEDRON, {14, 10, 6, 5}), Rule(TETRAHEDRON, {14, 10, 5, 1}),
     Rule(TETRAHEDRON, {14, 11, 2, 3}), Rule(TETRAHEDRON, {14, 11, 3, 7}),
     Rule(TETRAHEDRON, {14, 11, 7, 6}), Rule(TETRAHEDRON, {14, 11, 6, 2}),
     Rule(TETRAHEDRON, {14, 12, 3, 0}), Rule(TETRAHEDRON, {14, 12, 0, 4}),
     Rule(TETRAHEDRON, {14, 12, 4, 7}), Rule(TETRAHEDRON, {14, 12, 7, 3}),
     Rule(TETRAHEDRON, {14, 13, 4, 5}), Rule(TETRAHEDRON, {14, 13, 5, 6}),
     Rule(TETRAHEDRON, {14, 13, 6, 7}), Rule(TETRAHEDRON, {14, 13, 7, 4})};

const vector<Rule> *Hexahedron::Rhex_barycentric =
    RefineBarycentricCell(*Hexahedron::ReferenceHexahedron(), R_hex_barycentric);

vector<Rule> R20_hex = {Rule(HEXAHEDRON20, {0,  8,  20, 11, 12, 21, 26, 24, 27, 51,
                                            54, 34, 35, 56, 55, 71, 57, 60, 75, 73}),
                        Rule(HEXAHEDRON20, {8,  1,  9,  20, 21, 13, 22, 26, 28, 29,
                                            52, 51, 56, 37, 61, 55, 58, 62, 65, 60}),
                        Rule(HEXAHEDRON20, {11, 20, 10, 3,  24, 26, 23, 15, 54, 53,
                                            32, 33, 71, 55, 66, 41, 75, 70, 68, 72}),
                        Rule(HEXAHEDRON20, {20, 9,  2,  10, 26, 22, 14, 23, 52, 30,
                                            31, 53, 55, 61, 39, 66, 65, 63, 67, 70}),
                        Rule(HEXAHEDRON20, {12, 21, 26, 24, 4,  16, 25, 19, 57, 60,
                                            75, 73, 36, 59, 80, 74, 43, 76, 79, 50}),
                        Rule(HEXAHEDRON20, {21, 13, 22, 26, 16, 5,  17, 25, 58, 62,
                                            65, 60, 59, 38, 64, 80, 44, 45, 77, 76}),
                        Rule(HEXAHEDRON20, {24, 26, 23, 15, 19, 25, 18, 7,  75, 70,
                                            68, 72, 74, 80, 69, 42, 79, 78, 48, 49}),
                        Rule(HEXAHEDRON20, {26, 22, 14, 23, 25, 17, 6,  18, 65, 63,
                                            67, 70, 80, 64, 40, 69, 77, 46, 47, 78})};

const vector<Rule> *Hexahedron20::R20hex =
    RefineCell(*Hexahedron20::ReferenceHexahedron20(), R20_hex);

// ReferenceHexahedron27
const vector<Rule> R27_hex =
    {Rule(HEXAHEDRON27, {0,  8,  20, 11, 12, 21, 26, 24, 27,  51,  54, 34,  35, 56,
                         55, 71, 57, 60, 75, 73, 81, 85, 105, 108, 98, 109, 117}),
     Rule(HEXAHEDRON27, {8,  1,  9,  20, 21, 13, 22, 26, 28, 29,  52,  51,  56, 37,
                         61, 55, 58, 62, 65, 60, 84, 86, 89, 106, 105, 110, 118}),
     Rule(HEXAHEDRON27, {11, 20, 10, 3,  24, 26, 23, 15,  54,  53, 32, 33,  71, 55,
                         66, 41, 75, 70, 68, 72, 82, 108, 107, 94, 97, 112, 120}),
     Rule(HEXAHEDRON27, {20, 9,  2,  10, 26, 22, 14, 23,  52, 30, 31,  53,  55, 61,
                         39, 66, 65, 63, 67, 70, 83, 106, 90, 93, 107, 111, 119}),
     Rule(HEXAHEDRON27, {12, 21, 26, 24, 4,  16, 25,  19, 57,  60,  75, 73,  36, 59,
                         80, 74, 43, 76, 79, 50, 109, 88, 113, 116, 99, 101, 121}),
     Rule(HEXAHEDRON27, {21, 13, 22, 26, 16, 5,  17,  25, 58, 62,  65,  60,  59, 38,
                         64, 80, 44, 45, 77, 76, 110, 87, 92, 114, 113, 102, 122}),
     Rule(HEXAHEDRON27, {24, 26, 23, 15, 19, 25, 18,  7,   75,  70, 68,  72,  74, 80,
                         69, 42, 79, 78, 48, 49, 112, 116, 115, 95, 100, 104, 124}),
     Rule(HEXAHEDRON27, {26, 22, 14, 23, 25, 17, 6,   18,  65, 63, 67,  70,  80, 64,
                         40, 69, 77, 46, 47, 78, 111, 114, 91, 96, 115, 103, 123})};

// const vector<Rule> *Hexahedron27::R27hex =  RefineCell(*Hexahedron27::ReferenceHexahedron27(),
// R27_hex);
const vector<Rule> *Hexahedron27::R27hex = nullptr;

#endif
