#ifndef _VERTEX_H_
#define _VERTEX_H_

#include <unordered_map>

#include "Point.hpp"

#ifdef USE_DATAMESH

#include "DataSet.hpp"

#endif

class Buffer;

class Vertex {
  short part;
  short cnt;
public:
  Vertex();

  Vertex(short p);

  Vertex(const Vertex &v);

  int n() const { return cnt; }

  void inc() { ++cnt; }

  short dec() { return --cnt; }

  friend std::ostream &operator<<(std::ostream &s, const Vertex &V);

  friend std::ostream &operator<<(std::ostream &s, const Vertex *V);

  friend Buffer &operator<<(Buffer &b, const Vertex &V);

  friend Buffer &operator>>(Buffer &b, Vertex &V);

#ifdef USE_DATAMESH
private:
  DataContainer data{};
public:
  const DataContainer &GetData() { return data; }

  void SetData(const DataContainer &container) { data = container; }

#endif
};

class vertex : public std::unordered_map<Point, Vertex *>::const_iterator {
  typedef std::unordered_map<Point, Vertex *>::const_iterator Iterator;
public:
  vertex(Iterator v) : Iterator(v) {}

  const Point &operator()() const { return (*this)->first; }

  double operator[](int i) const { return (*this)->first[i]; }

#ifdef USE_DATAMESH

  const DataContainer &GetData() const { return (*this)->second->GetData(); }

  void SetData(const DataContainer &container) { (*this)->second->SetData(container); }

#endif
};

inline std::ostream &operator<<(std::ostream &s, const vertex &v) {
  return s << v->first << " : " << *(v->second) << endl;
}

#endif
