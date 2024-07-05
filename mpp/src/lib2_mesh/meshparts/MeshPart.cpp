#include "MeshPart.hpp"

#include "Buffer.hpp"
#include "Mesh.hpp"

Buffer &operator>>(Buffer &b, BoundaryFaces &B) {
  Point z;
  int part;
  b >> z >> part;
  B.Insert(z, part);
  return b;
}

CellBoundaryFaces::CellBoundaryFaces(const Mesh &M, const cell &c) : N(0) {
  for (int i = 0; i < c.Faces(); ++i) {
    bf[N] = M.find_bnd_face(c.Face(i));
    if (bf[N] != M.bnd_faces_end()) ++N;
  }
}

Buffer &operator<<(Buffer &b, const CellBoundaryFaces &bf) {
  b << bf.size();
  for (int i = 0; i < bf.size(); ++i)
    b << bf[i];
  return b;
}