#include "CurlElement.hpp"


using namespace std;

int FaceOrientation(const Cell &c, int f) {
  return 2 * ((f == 0) ? (c[0] > c[1]) : (c.FaceCorner(f, 1) > c.FaceCorner(f, 0))) - 1;
}

int EdgeOrientation(const Cell &c, int e) {
  if (c.EdgeCorner(e, 0) < c.EdgeCorner(e, 1)) return 1;
  return -1;
}

int EdgeCornerIndex(const int edge, const int corner) {
  switch (edge) {
  case 0: {
    if (corner == 0) return 0;
    else return 1;
  }
  case 1: {
    if (corner == 0) return 0;
    else return 2;
  }
  case 2: {
    if (corner == 0) return 0;
    else return 3;
  }
  case 3: {
    if (corner == 0) return 1;
    else return 2;
  }
  case 4: {
    if (corner == 0) return 1;
    else return 3;
  }
  case 5: {
    if (corner == 0) return 2;
    else return 3;
  }
  }
  THROW("Not implemented!")
  return 0;
}

int FaceCornerIndex(const int face, const int corner) {
  switch (face) {
  case 0:
    switch (corner) {
    case 0:
      return 0;
    case 1:
      return 1;
    case 2:
      return 2;
    }
  case 1:
    switch (corner) {
    case 0:
      return 0;
    case 1:
      return 1;
    case 2:
      return 3;
    }
  case 2:
    switch (corner) {
    case 0:
      return 0;
    case 1:
      return 2;
    case 2:
      return 3;
    }
  case 3:
    switch (corner) {
    case 0:
      return 1;
    case 1:
      return 2;
    case 2:
      return 3;
    }
  }
  THROW("Not implemented!")
  return 0;
}

CurlElement::CurlElement(const VectorMatrixBase &base, const Cell &c) :
    Element(base, c), S(base.GetDisc().GetShape(c)), vectorfield(this->nQ(), this->size()),
    curlfield(this->nQ(), this->size()) {
  if (base.DoFsName() == "Curl2TetDoF") {
    vector<Point> vertices(4);

    for (int i = 0; i < 4; ++i)
      vertices[i] = c.Corner(i);

    for (int i = 0; i < 4; ++i)
      for (int j = i + 1; j < 4; ++j)
        if (vertices[i] > vertices[j]) {
          Point tmp = vertices[i];

          vertices[i] = vertices[j];
          vertices[j] = tmp;
        }

    Tensor F;

    for (int i = 0; i < F.Dim(); ++i)
      for (int j = 0; j < F.Dim(); ++j)
        F[i][j] = vertices[j + 1][i] - vertices[0][i];

    Scalar detF = F[0][0] * F[1][1] * F[2][2] + F[0][1] * F[1][2] * F[2][0]
                  + F[0][2] * F[1][0] * F[2][1] - F[0][1] * F[1][0] * F[2][2]
                  - F[0][2] * F[1][1] * F[2][0] - F[0][0] * F[1][2] * F[2][1];
    Tensor T;

    T[0][0] = (F[1][1] * F[2][2] - F[1][2] * F[2][1]) / detF;
    T[0][1] = (F[1][2] * F[2][0] - F[1][0] * F[2][2]) / detF;
    T[0][2] = (F[1][0] * F[2][1] - F[1][1] * F[2][0]) / detF;
    T[1][0] = (F[0][2] * F[2][1] - F[0][1] * F[2][2]) / detF;
    T[1][1] = (F[0][0] * F[2][2] - F[0][2] * F[2][0]) / detF;
    T[1][2] = (F[0][1] * F[2][0] - F[0][0] * F[2][1]) / detF;
    T[2][0] = (F[0][1] * F[1][2] - F[0][2] * F[1][1]) / detF;
    T[2][1] = (F[0][2] * F[1][0] - F[0][0] * F[1][2]) / detF;
    T[2][2] = (F[0][0] * F[1][1] - F[0][1] * F[1][0]) / detF;

    int n_edges = c.Edges();
    int sz = int(size());
    RMatrix matrix(sz);

    for (int edge = 0; edge < n_edges; ++edge) {
      Point e_0 = c.LocalCorner(EdgeCornerIndex(edge, 0));
      Point e_1 = c.LocalCorner(EdgeCornerIndex(edge, 1));

      /*
      Point e0 = c.EdgeCorner (edge, 0);
      Point e1 = c.EdgeCorner (edge, 1);
      if (e1 > e0)
          tangent[2 * edge] = e1 - e0;
      else
          tangent[2 * edge] = e0 - e1;
      */

      tangent[2 * edge] = vertices[EdgeCornerIndex(edge, 1)] - vertices[EdgeCornerIndex(edge, 0)];
      tangent[2 * edge + 1] = tangent[2 * edge];

      Point z_0 = (0.5 + 0.5 / sqrt(3.0)) * e_0 + (0.5 - 0.5 / sqrt(3.0)) * e_1;
      Point z_1 = (0.5 - 0.5 / sqrt(3.0)) * e_0 + (0.5 + 0.5 / sqrt(3.0)) * e_1;

      for (int i = 0; i < size(); ++i) {
        VectorField shape_value = S.LocalVector(z_0, i);
        VectorField Tshape_value;

        for (int j = 0; j < Tshape_value.Dim(); ++j) {
          Tshape_value[j] = 0.0;

          for (int k = 0; k < Tshape_value.Dim(); ++k)
            Tshape_value[j] += T[j][k] * shape_value[k];
        }

        matrix[2 * edge][i] = Tshape_value * tangent[2 * edge];
        shape_value = S.LocalVector(z_1, i);

        for (int j = 0; j < Tshape_value.Dim(); ++j) {
          Tshape_value[j] = 0.0;

          for (int k = 0; k < T.Dim(); ++k)
            Tshape_value[j] += T[j][k] * shape_value[k];
        }

        matrix[2 * edge + 1][i] = Tshape_value * tangent[2 * edge];
      }
    }

    for (int face = 0; face < c.Faces(); ++face) {
      Point f_0 = vertices[FaceCornerIndex(face, 0)];
      Point f_1 = vertices[FaceCornerIndex(face, 1)];
      Point f_2 = vertices[FaceCornerIndex(face, 2)];

      tangent[2 * (n_edges + face)] = f_1 - f_0;
      tangent[2 * (n_edges + face) + 1] = f_2 - f_0;

      Point f_0_local = c.LocalCorner(FaceCornerIndex(face, 0));
      Point f_1_local = c.LocalCorner(FaceCornerIndex(face, 1));
      Point f_2_local = c.LocalCorner(FaceCornerIndex(face, 2));
      Point z_01 = 0.5 * (f_0_local + f_1_local);
      Point z_02 = 0.5 * (f_0_local + f_2_local);
      Point z_12 = 0.5 * (f_1_local + f_2_local);

      for (int i = 0; i < 2; ++i)
        for (int j = 0; j < size(); ++j) {
          VectorField shape_value =
              S.LocalVector(z_01, j) + S.LocalVector(z_02, j) + S.LocalVector(z_12, j);
          VectorField Tshape_value;

          for (int k = 0; k < Tshape_value.Dim(); ++k) {
            Tshape_value[k] = 0.0;

            for (int l = 0; l < T.Dim(); ++l)
              Tshape_value[k] += T[k][l] * shape_value[l];
          }

          matrix[2 * (face + n_edges) + i][j] =
              (Tshape_value * tangent[2 * (face + n_edges) + i]) / 3.0;
        }
    }

    matrix.Invert();

    for (int q = 0; q < nQ(); ++q) {
      for (int i = 0; i < size(); ++i) {
        curlfield[q][i] = Point(0.0, 0.0, 0.0);
        vectorfield[q][i] = Point(0.0, 0.0, 0.0);

        for (int j = 0; j < size(); ++j) {
          VectorField tmp_curl = matrix[j][i] * S.LocalCurl(LocalQPoint(q), j);
          VectorField tmp_value = matrix[j][i] * S.LocalVector(LocalQPoint(q), j);

          for (int k = 0; k < F.Dim(); ++k)
            for (int l = 0; l < F.Dim(); ++l) {
              curlfield[q][i][k] += F[k][l] / detF * tmp_curl[l];
              vectorfield[q][i][k] += T[k][l] * tmp_value[l];
            }
        }
      }
    }
  } else {
    for (int i = 0; i < c.Edges(); ++i) {
      sign[i] = 2 * (c.EdgeCorner(i, 0) < c.EdgeCorner(i, 1)) - 1;
      tangent[i] = sign[i] * (c.EdgeCorner(i, 1) - c.EdgeCorner(i, 0));
    }
    if (base.DoFsName() == "Curl2DoF") {
      int n = c.Edges();
      for (int i = n; i < size(); ++i)
        sign[i] = 1;
      for (int i = 0; i < c.Faces(); ++i)
        sign[2 * n + 4 * i] = FaceOrientation(c, i);
    }
    for (int q = 0; q < nQ(); ++q) {
      double idet = 1 / GetTransformationDet(q);
      for (int i = 0; i < size(); ++i) {
        vectorfield[q][i] = sign[i] * (GetTransformation(q) * S.LocalVector(LocalQPoint(q), i));
        curlfield[q][i] =
            sign[i] * idet * GetTransformation(q).ApplyJ(S.LocalCurl(LocalQPoint(q), i));
      }
      /*
      for (int i=0; i<12; ++i) {
          vectorfield[q][i] = sign[i] *
              (GetTransformation(q) * S.LocalVector(LocalQPoint(q),i));
          curlfield[q][i] = sign[i] * idet *
              GetTransformation(q).ApplyJ(S.LocalCurl(LocalQPoint(q),i));
      for (int f=0; f<4; ++f) {
          vectorfield[q][12+2*f] =
              (GetTransformation(q) * (
              TT11[f] * S.LocalVector(LocalQPoint(q),12+2*f)
              + TT12[f] * S.LocalVector(LocalQPoint(q),12+2*f+1) );
          vectorfield[q][12+2*f+1] =
              (GetTransformation(q) * (
              TT21[f] * S.LocalVector(LocalQPoint(q),12+2*f)
              + TT22[f] * S.LocalVector(LocalQPoint(q),12+2*f+1) );

      }

      */
    }
  }
}

Point CurlElement::ReferenceTangent(int i) const {
  switch (i) {
  case (0):
    return Point(1, 0, 0);
  case (1):
    return Point(1, 0, 0);
  case (2):
    return Point(1, -1, 0);
  case (3):
    return Point(1, -1, 0);
  case (4):
    return Point(0, 1, 0);
  case (5):
    return Point(0, 1, 0);
  case (6):
    return Point(0, 0, 1);
  case (7):
    return Point(0, 0, 1);
  case (8):
    return Point(1, 0, -1);
  case (9):
    return Point(1, 0, -1);
  case (10):
    return Point(0, 1, -1);
  case (11):
    return Point(0, 1, -1);
  case (12):
    return Point(0, 1, 0);
  case (13):
    return Point(1, 0, 0);
  case (14):
    return Point(0, 1, -1);
  case (15):
    return Point(1, 0, -1);
  case (16):
    return Point(0, 0, 1);
  case (17):
    return Point(0, 1, 0);
  case (18):
    return Point(0, 0, 1);
  case (19):
    return Point(1, 0, 0);
  }
  THROW("Not implemented!")
  return Origin;
}

Point CurlElement::DofShape_int(double x, double y, double z, int i) const {
  // functions for DOF calculation of Curl2DoF
  switch (i) {
  // 2-nd order face based
  // face (1, 0, 3, 2)
  case (24):
    return Point(1.5, 0, 0);
  case (25):
    return Point(0, -1.5, 0);
  case (26):
    return Point(2.25 - 4.5 * x, -2.25 + 4.5 * y, 0);
  case (27):
    return Point(2.25 - 4.5 * x, 2.25 - 4.5 * y);
    // face (0, 1, 5, 4)
  case (28):
    return Point(-1.5, 0, 0);
  case (29):
    return Point(0, 0, -1.5);
  case (30):
    return Point(2.25 - 4.5 * x, 0, -2.25 + 4.5 * z);
  case (31):
    return Point(2.25 - 4.5 * x, 0, 2.25 - 4.5 * z);
    // face (1, 2, 6, 5)
  case (32):
    return Point(0, -1.5, 0);
  case (33):
    return Point(0, 0, -1.5);
  case (34):
    return Point(0, 2.25 - 4.5 * y, -2.25 + 4.5 * z);
  case (35):
    return Point(0, 2.25 - 4.5 * y, 2.25 - 4.5 * z);
    // face (2, 3, 7, 6)
  case (36):
    return Point(1.5, 0, 0);
  case (37):
    return Point(0, 0, -1.5);
  case (38):
    return Point(2.25 - 4.5 * x, 0, -2.25 + 4.5 * z);
  case (39):
    return Point(2.25 - 4.5 * x, 0, 2.25 - 4.5 * z);
    // face (3, 0, 4, 7)
  case (40):
    return Point(0, 1.5, 0);
  case (41):
    return Point(0, 0, -1.5);
  case (42):
    return Point(0, 2.25 - 4.5 * y, -2.25 + 4.5 * z);
  case (43):
    return Point(0, 2.25 - 4.5 * y, 2.25 - 4.5 * z);
    // face (4, 5, 6, 7)
  case (44):
    return Point(-1.5, 0, 0);
  case (45):
    return Point(0, -1.5, 0);
  case (46):
    return Point(2.25 - 4.5 * x, -2.25 + 4.5 * y, 0);
  case (47):
    return Point(2.25 - 4.5 * x, 2.25 - 4.5 * y, 0);
    // 2-nd order cell-based
  case (48):
    return Point(9, 0, 0);
  case (49):
    return Point(0, 9, 0);
  case (50):
    return Point(0, 0, 9);
  case (51):
    return Point(-6.75 + 13.5 * x, 6.75 - 13.5 * y, 0);
  case (52):
    return Point(-6.75 + 13.5 * x, 0, 6.75 - 13.5 * z);
  case (53):
    return Point(0, -6.75 + 13.5 * y, -6.75 + 13.5 * z);
  }
  THROW("Not implemented!")
  return Origin;
}

Point CurlElement::Edge_Corner(int _i, int _j) const {
  vector<Point> vertices(4);

  for (int i = 0; i < 4; ++i)
    vertices[i] = c.Corner(i);

  for (int i = 0; i < 4; ++i)
    for (int j = i + 1; j < 4; ++j)
      if (vertices[i] > vertices[j]) {
        Point tmp = vertices[i];

        vertices[i] = vertices[j];
        vertices[j] = tmp;
      }

  return vertices[EdgeCornerIndex(_i, _j)];
}

Point CurlElement::Face_Corner(int _i, int _j) const {
  vector<Point> vertices(4);

  for (int i = 0; i < 4; ++i)
    vertices[i] = c.Corner(i);

  for (int i = 0; i < 4; ++i)
    for (int j = i + 1; j < 4; ++j)
      if (vertices[i] > vertices[j]) {
        Point tmp = vertices[i];

        vertices[i] = vertices[j];
        vertices[j] = tmp;
      }

  return vertices[FaceCornerIndex(_i, _j)];
}

Point CurlElement::Q_Point(int i) const {
  const Point q = this->LocalQPoint(i);
  vector<Point> vertices(4);

  for (int ii = 0; ii < 4; ++ii)
    vertices[ii] = c.Corner(ii);

  for (int ii = 0; ii < 4; ++ii)
    for (int j = ii + 1; j < 4; ++j)
      if (vertices[ii] > vertices[j]) {
        Point tmp = vertices[ii];

        vertices[ii] = vertices[j];
        vertices[j] = tmp;
      }
  return (1 - q[0] - q[1] - q[2]) * vertices[0] + q[0] * vertices[1] + q[1] * vertices[2]
         + q[2] * vertices[3];
}

Point CurlElement::LocalQ_Point(int i) const {
  const Point q = this->LocalQPoint(i);
  vector<Point> vertices(4);

  for (int ii = 0; ii < 4; ++ii)
    vertices[ii] = c.Corner(ii);

  for (int ii = 0; ii < 4; ++ii)
    for (int j = ii + 1; j < 4; ++j)
      if (vertices[ii] > vertices[j]) {
        Point tmp = vertices[ii];

        vertices[ii] = vertices[j];
        vertices[j] = tmp;
      }
  return (1 - q[0] - q[1] - q[2]) * vertices[0] + q[0] * vertices[1] + q[1] * vertices[2]
         + q[2] * vertices[3];
}

Point CurlElement::Nodal_Point(int i) const {
  if (size() == 20) {
    vector<Point> vertices(4);

    for (int ii = 0; ii < 4; ++ii)
      vertices[ii] = c.Corner(ii);

    for (int ii = 0; ii < 4; ++ii)
      for (int j = ii + 1; j < 4; ++j)
        if (vertices[ii] > vertices[j]) {
          Point tmp = vertices[ii];

          vertices[ii] = vertices[j];
          vertices[j] = tmp;
        }
    if (i < 2 * c.Edges()) {
      int edge = i / 2;
      if (i == 2 * edge)
        return 0.5 * (vertices[EdgeCornerIndex(edge, 1)] + vertices[EdgeCornerIndex(edge, 0)]);
      else
        return 0.25 * (vertices[EdgeCornerIndex(edge, 1)] + vertices[EdgeCornerIndex(edge, 0)])
               + 0.5 * vertices[EdgeCornerIndex(edge, 1)];
    }
    int face = (i - 2 * c.Edges()) / 2;
    Point f_0 = vertices[FaceCornerIndex(face, 0)];
    Point f_1 = vertices[FaceCornerIndex(face, 1)];
    Point f_2 = vertices[FaceCornerIndex(face, 2)];
    if (2 * face == (i - 2 * c.Edges())) return (1 / 3.0) * (f_0 + f_1 + f_2);
    else return (0.5 / 3.0) * (f_0 + f_1 + f_2) + 0.5 * f_0;
  } else return (*this)[i]();
}

VectorField CurlElement::VectorValue(const Point &z, int i) const {
  if (size() != 20) return sign[i] * (c.GetTransformation(z) * S.LocalVector(z, i));

  vector<Point> vertices(4);

  for (int j = 0; j < 4; ++j)
    vertices[j] = c.Corner(j);

  for (int j = 0; j < 4; ++j)
    for (int k = j + 1; k < 4; ++k)
      if (vertices[j] > vertices[k]) {
        Point tmp = vertices[j];

        vertices[j] = vertices[k];
        vertices[k] = tmp;
      }

  Tensor F;

  for (int j = 0; j < F.Dim(); ++j)
    for (int k = 0; k < F.Dim(); ++k)
      F[j][k] = vertices[k + 1][j] - vertices[0][j];

  Scalar detF = F[0][0] * F[1][1] * F[2][2] + F[0][1] * F[1][2] * F[2][0]
                + F[0][2] * F[1][0] * F[2][1] - F[0][1] * F[1][0] * F[2][2]
                - F[0][2] * F[1][1] * F[2][0] - F[0][0] * F[1][2] * F[2][1];
  Tensor T;

  T[0][0] = (F[1][1] * F[2][2] - F[1][2] * F[2][1]) / detF;
  T[0][1] = (F[1][2] * F[2][0] - F[1][0] * F[2][2]) / detF;
  T[0][2] = (F[1][0] * F[2][1] - F[1][1] * F[2][0]) / detF;
  T[1][0] = (F[0][2] * F[2][1] - F[0][1] * F[2][2]) / detF;
  T[1][1] = (F[0][0] * F[2][2] - F[0][2] * F[2][0]) / detF;
  T[1][2] = (F[0][1] * F[2][0] - F[0][0] * F[2][1]) / detF;
  T[2][0] = (F[0][1] * F[1][2] - F[0][2] * F[1][1]) / detF;
  T[2][1] = (F[0][2] * F[1][0] - F[0][0] * F[1][2]) / detF;
  T[2][2] = (F[0][0] * F[1][1] - F[0][1] * F[1][0]) / detF;
  int n_edges = c.Edges();
  int sz = int(size());
  RMatrix matrix(sz);
  // 	SmallMatrix matrix (size ());

  for (int edge = 0; edge < n_edges; ++edge) {
    Point e_0 = c.LocalCorner(EdgeCornerIndex(edge, 0));
    Point e_1 = c.LocalCorner(EdgeCornerIndex(edge, 1));
    Point z_0 = (0.5 + 0.5 / sqrt(3.0)) * e_0 + (0.5 - 0.5 / sqrt(3.0)) * e_1;
    Point z_1 = (0.5 - 0.5 / sqrt(3.0)) * e_0 + (0.5 + 0.5 / sqrt(3.0)) * e_1;

    for (int j = 0; j < size(); ++j) {
      matrix[2 * edge][j] = S.LocalVector(z_0, j) * tangent[2 * edge];
      matrix[2 * edge + 1][j] = S.LocalVector(z_1, j) * tangent[2 * edge];
    }
  }

  for (int face = 0; face < c.Faces(); ++face) {
    Point f_0 = c.LocalCorner(FaceCornerIndex(face, 0));
    Point f_1 = c.LocalCorner(FaceCornerIndex(face, 1));
    Point f_2 = c.LocalCorner(FaceCornerIndex(face, 2));
    Point z_01 = 0.5 * (f_0 + f_1);
    Point z_02 = 0.5 * (f_0 + f_2);
    Point z_12 = 0.5 * (f_1 + f_2);

    for (int j = 0; j < 2; ++j)
      for (int k = 0; k < size(); ++k)
        matrix[2 * (face + n_edges) + j][k] =
            ((S.LocalVector(z_01, k) + S.LocalVector(z_02, k) + S.LocalVector(z_12, k))
             * tangent[2 * (face + n_edges) + j])
            / 3.0;
  }

  matrix.Invert();

  VectorField value = VectorField(0.0, 0.0, 0.0);

  for (int j = 0; j < size(); ++j) {
    VectorField tmp_value = matrix[j][i] * S.LocalVector(z, j);

    for (int k = 0; k < T.Dim(); ++k)
      for (int l = 0; l < T.Dim(); ++l)
        value[k] += T[k][l] * tmp_value[l];
  }

  return value;
}

VectorField CurlElement::VectorValue(const Point &z, const Vector &u, int k) const {
  VectorField V;
  for (int i = 0; i < size(); ++i)
    V += u(r(i), k) * sign[i] * (c.GetTransformation(z) * S.LocalVector(z, i));
  return V;
}

VectorField CurlElement::VectorValue(int q, const Vector &u, int k) const {
  VectorField V;
  for (int i = 0; i < size(); ++i)
    V += u(Nodal_Point(i), k) * vectorfield[q][i];
  //	V += u(r(i),k) * vectorfield[q][i];
  return V;
}

VectorField CurlElement::CurlVector(int q, const Vector &u, int k) const {
  VectorField V;
  for (int i = 0; i < size(); ++i)
    V += u(Nodal_Point(i), k) * curlfield[q][i];
  //    for (int i=0; i<size(); ++i) V += u(r(i),k) * curlfield[q][i];
  return V;
}

std::ostream &operator<<(std::ostream &s, const CurlElement &E) { return s << E.S; }
