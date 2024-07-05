#include "Face.hpp"

#include "Buffer.hpp"

std::ostream &operator<<(std::ostream &s, const Face &face) {
  return s << "(" << face.centerLeftCell << "," << face.centerRightCell << ")";
}

std::ostream &operator<<(std::ostream &s, const Face *F) { return s << *F; }

Buffer &operator<<(Buffer &b, const Face &face) {
  return b << face.centerLeftCell << face.centerRightCell;
}

Buffer &operator>>(Buffer &b, Face &face) {
  //  Point x, y;
  //  b >> x >> y;
  //  F = Face(x, y);
  return b >> face.centerLeftCell >> face.centerRightCell;
}

Buffer &operator<<(Buffer &b, const bnd_face &B) { return b << B() << B.Part(); }