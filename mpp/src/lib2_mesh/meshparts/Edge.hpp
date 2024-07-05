#ifndef _EDGE_H_
#define _EDGE_H_

#include <unordered_map>

#include "Point.hpp"

class Buffer;

class Edge {
  Point left;
  Point right;
  short cnt;
public:
  Edge() : left(Infty), right(Infty), cnt(0) {}

  Edge(const Point &l, const Point &r);

  int n() const { return cnt; }

  void inc() { ++cnt; }

  short dec() { return --cnt; }

  const Point &Left() const { return left; }

  const Point &Right() const { return right; }

  Point operator()() const { return 0.5 * (left + right); };

  friend std::ostream &operator<<(std::ostream &s, const Edge &E);

  friend std::ostream &operator<<(std::ostream &s, const Edge *E);

  friend Buffer &operator<<(Buffer &b, const Edge &F);

  friend Buffer &operator>>(Buffer &b, Edge &F);
};

class edge : public std::unordered_map<Point, Edge *>::const_iterator {
  typedef std::unordered_map<Point, Edge *>::const_iterator Iterator;
public:
  edge(Iterator e) : Iterator(e) {}

  const Point &operator()() const { return (*this)->first; }

  const Point &Center() const { return (*this)->first; }

  const Point &Left() const { return (*this)->second->Left(); }

  const Point &Right() const { return (*this)->second->Right(); }

  Point Tangent() const { return Right() - Left(); }

  double length() const { return dist(Left(), Right()); }

  friend std::ostream &operator<<(std::ostream &s, const edge &E) {
    return s << E->first << " : " << *(E->second) << endl;
  }

  bool isEqual(const edge &other) {
    return Center() == other.Center() && Left() == other.Left() && Right() == other.Right();
  }
};

#endif // of #ifndef _EDGE_H_
