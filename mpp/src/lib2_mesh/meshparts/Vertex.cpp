#include "Vertex.hpp"

#include "Buffer.hpp"

Vertex::Vertex() {
  part = 0;
  cnt = 0;
}

Vertex::Vertex(short p) {
  part = p;
  cnt = 1;
}

Vertex::Vertex(const Vertex &v) {
  part = v.part;
  cnt = v.cnt;
}

std::ostream &operator<<(std::ostream &s, const Vertex &V) {
  return s << V.part << " [" << V.cnt << "]";
}

std::ostream &operator<<(std::ostream &s, const Vertex *V) { return s << *V; }

Buffer &operator<<(Buffer &b, const Vertex &V) { return b << V.part << V.cnt; }

Buffer &operator>>(Buffer &b, Vertex &V) { return b >> V.part >> V.cnt; }