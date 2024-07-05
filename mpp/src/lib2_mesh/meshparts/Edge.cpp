#include "Edge.hpp"

#include "Buffer.hpp"

Edge::Edge(const Point &l, const Point &r) : cnt(1) {
  if (l < r) {
    left = l;
    right = r;
  } else {
    left = r;
    right = l;
  }
}

std::ostream &operator<<(std::ostream &s, const Edge &E) {
  return s << "(" << E.left << "," << E.right << ") [" << E.cnt << "]";
}

std::ostream &operator<<(std::ostream &s, const Edge *E) { return s << *E; }

Buffer &operator<<(Buffer &b, const Edge &E) { return b << E.left << E.right << E.cnt; }

Buffer &operator>>(Buffer &b, Edge &E) { return b >> E.left >> E.right >> E.cnt; }
