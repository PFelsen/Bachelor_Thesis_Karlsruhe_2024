#include "NodalPointProvider.hpp"
#include "GaussLobattoNodalPoints.hpp"
#include "LagrangeNodalPoints.hpp"

std::vector<Point> EquidistantNodalPointProvider::GetNodalPoints(CELLTYPE ctype, int degree) const {
  vector<Point> points(EquidistantNodalPointCount(ctype, degree));
  const Cell *c = ReferenceCell(ctype);
  EquidistantNodalPoints(*c, degree, points);
  return points;
}

std::vector<Point> GaussLobattoNodalPointProvider::GetNodalPoints(CELLTYPE ctype,
                                                                  int degree) const {
  const Cell *c = ReferenceCell(ctype);
  switch (SpaceCellType(ctype)) {
  case INTERVAL:
    return GaussLobattoNodalPointsInterval<double, SpaceDimension, TimeDimension>(*c, degree + 1);
  case SPACETIME_INTERVAL:
  case QUADRILATERAL:
    return GaussLobattoNodalPointsQuadrilateral<double, SpaceDimension, TimeDimension>(*c,
                                                                                       degree + 1);
  case SPACETIME_QUADRILATERAL:
  case HEXAHEDRON:
    return GaussLobattoNodalPointsHexaedron<double, SpaceDimension, TimeDimension>(*c, degree + 1);
  default:
    THROW("Celltype " + std::to_string(ctype) + " not implemented.")
  }
}

int EquidistantNodalPointCount(CELLTYPE type, int degree) {
  if (degree == 0) return 1;
  switch (type) {
  case INTERVAL:
    return degree + 1;
  case TRIANGLE:
    return ((degree + 1) * (degree + 2)) / 2;
  case SPACETIME_INTERVAL:
  case QUADRILATERAL:
    return (degree + 1) * (degree + 1);
  case TETRAHEDRON:
    return ((degree + 1) * (degree + 2) * (degree + 3)) / 6;
  case SPACETIME_QUADRILATERAL:
  case HEXAHEDRON:
    return (degree + 1) * (degree + 1) * (degree + 1);
  case SPACETIME_TRIANGLE:
    return EquidistantNodalPointCount(TRIANGLE, degree)
           * EquidistantNodalPointCount(INTERVAL, degree);
  default:
    THROW("Cell type not implemented in "
          "EquidistantNodalPointCount(CELLTYPE type, int degree) for "
          "type="
          + std::to_string(type) + " degree=" + std::to_string(degree));
  }
}

int EquidistantNodalPointCountOnFace(CELLTYPE type, int degree, int face) {
  if (degree == 0) return 0;
  switch (type) {
  case INTERVAL:
    return 1;
  case TRIANGLE:
  case QUADRILATERAL:
    return degree + 1;
  case TETRAHEDRON:
    return ((degree + 1) * (degree + 2)) / 2;
  case SPACETIME_QUADRILATERAL:
  case HEXAHEDRON:
    return (degree + 1) * (degree + 1);
  }
  THROW("Cell type not implemented in LagrangeDoF::EquidistantNodalPointCountOnFace")
}

int EquidistantNodalPointCountOnEdge(CELLTYPE type, int degree, int edge) {
  if (degree == 0) return 0;
  switch (type) {
  case INTERVAL:
    return 1;
  case TRIANGLE:
  case QUADRILATERAL:
  case TETRAHEDRON:
  case SPACETIME_QUADRILATERAL:
  case HEXAHEDRON:
    return (degree + 1);
  }
  THROW("Cell type not implemented in LagrangeDoF::EquidistantNodalPointCountOnEdge")
}

int EquidistantNodalPointIdOnEdge(CELLTYPE type, int degree, int edgeId, int k) {
  if (degree == 0) return -1;
  switch (type) {
  case INTERVAL:
  case TRIANGLE:
  case QUADRILATERAL:
    return EquidistantNodalPointIdOnFace(type, degree, edgeId, k);
  case TETRAHEDRON:
    if (degree < 5) return tetrahdronEdgeNodalMap[degree - 1][edgeId][k];
    else THROW("NodalPointsOnFace not known for degree higher than 3")
    // THROW("EquidistantNodalPointIdOnEdge not known for TETRAHEDRON")
  case SPACETIME_QUADRILATERAL:
  case HEXAHEDRON:
    int N = degree + 1;
    int nodalpointsOnFace = EquidistantNodalPointCountOnFace(HEXAHEDRON, degree, 0);
    int nodalpointIdQuad = EquidistantNodalPointIdOnEdge(QUADRILATERAL, degree, edgeId % 4, k);
    switch (edgeId) {
    case 0:
    case 1:
    case 2:
    case 3:
      return nodalpointIdQuad;
    case 4:
      return k * nodalpointsOnFace;
    case 5:
      return k * nodalpointsOnFace + degree;
    case 6:
      return (k + 1) * nodalpointsOnFace - 1;
    case 7:
      return (k + 1) * nodalpointsOnFace - N;
    case 8:
    case 9:
    case 10:
    case 11:
      return degree * nodalpointsOnFace + nodalpointIdQuad;
    default:
      THROW("Edge index " + std::to_string(edgeId) + " unknown.")
    }
  }
  THROW("Cell type not implemented in LagrangeDoF")
}

int EquidistantNodalPointIdOnFace(CELLTYPE type, int degree, int face, int k) {
  if (degree == 0) return -1;
  int N = degree + 1;
  switch (type) {
  case INTERVAL:
    return face * degree;
  case TRIANGLE:
    if (k == degree) return (face + 1) % 3;
    return 3 * k + face;
  case QUADRILATERAL:
    switch (face) {
    case 0:
      return k;
    case 1:
      return degree + (degree + 1) * k;
    case 2:
      return degree * (degree + 1) + k;
    case 3:
      return (degree + 1) * k;
    default:
      THROW("Face index " + std::to_string(face) + " unknown.")
    }
  case TETRAHEDRON:
    if (degree < 5) return tetrahdronNodalMap[degree - 1][face][k];
    else THROW("NodalPointsOnFace not known for degree higher than 4")
  case SPACETIME_QUADRILATERAL:
    // mapping from face of hexahedron to spacetime_quad
    // int f2f[] = {4,0,1,2,3,5};
    switch (face) {
    case 4:
      return (k % N) + (k / N) * N;
    case 0:
      return (k % N) + (k / N) * N * N;
    case 1:
      return degree + (k % N) * N + (k / N) * N * N;
    case 2:
      return degree * N + (k % N) + (k / N) * N * N;
    case 3:
      return (k % N) * N + (k / N) * N * N;
    case 5:
      return degree * N * N + (k % N) + (k / N) * N;
    default:
      THROW("Face index " + std::to_string(face) + " unknown.")
    }
  case HEXAHEDRON:
    switch (face) {
    case 0:
      return (k % N) + (k / N) * N;
    case 1:
      return (k % N) + (k / N) * N * N;
    case 2:
      return degree + (k % N) * N + (k / N) * N * N;
    case 3:
      return degree * N + (k % N) + (k / N) * N * N;
    case 4:
      return (k % N) * N + (k / N) * N * N;
    case 5:
      return degree * N * N + (k % N) + (k / N) * N;
    default:
      THROW("Face index " + std::to_string(face) + " unknown.")
    }
  }
  THROW("Cell type not implemented in LagrangeDoF")
}

void EquidistantNodalPoints(const Cell &c, int degree, std::vector<Point> &z) {
  z.resize(EquidistantNodalPointCount(c.ReferenceType(), degree));
  if (degree == 0) {
    z[0] = c();
    return;
  }
  switch (c.ReferenceType()) {
  case INTERVAL:
    LagrangeNodalPointsInterval(c, z, degree);
    break;
  case TRIANGLE:
    LagrangeNodalPointsTriangle(c, z, degree);
    break;
  case QUADRILATERAL:
  case SPACETIME_INTERVAL:
    LagrangeNodalPointsQuadrilateral(c, z, degree);
    break;
  case TETRAHEDRON:
    LagrangeNodalPointsTetrahedron(c, z, degree);
    break;
  case HEXAHEDRON:
  case SPACETIME_QUADRILATERAL:
    LagrangeNodalPointsHexahedron(c, z, degree);
    break;
  case SPACETIME_TRIANGLE:
    LagrangeNodalPointsTrianglePrism(c, z, degree);
    break;
  default:
    THROW("Celltype for EquidistantNodalPoints not implemented: "
          + std::to_string(c.ReferenceType()))
  }
}

std::vector<Point> EquidistantNodalPoints(const Cell &c, int degree) {
  std::vector<Point> nodalPoints;
  EquidistantNodalPoints(c, degree, nodalPoints);
  return nodalPoints;
}
