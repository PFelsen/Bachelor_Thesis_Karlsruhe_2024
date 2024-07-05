#ifndef MESHPART_H
#define MESHPART_H

#include "Cell.hpp"
#include "Edge.hpp"
#include "Face.hpp"
#include "Vertex.hpp"

#include <unordered_map>
#include "Parallel.hpp"

template<typename P, typename T>
class MeshPart : public std::unordered_map<P, T> {
public:
  typename std::unordered_map<P, T>::const_iterator Begin() const { return (*this).begin(); }

  typename std::unordered_map<P, T>::const_iterator End() const { return (*this).end(); }

  typename std::unordered_map<P, T>::const_iterator Find(const P &key) const {
    return (*this).find(key);
  }

  typename std::unordered_map<P, T>::const_iterator Insert(const P &key, T value) {
    return (*this).try_emplace(key, value).first;
  };

  typename std::unordered_map<P, T>::const_iterator Replace(const P &key, T value) {
    Remove(key);
    return (*this).try_emplace(key, value).first;
  };

  bool Remove(const P &key) {
    auto elem = (*this).find(key);
    if (elem == (*this).end()) { return false; }
    if constexpr (std::is_pointer<T>::value) { delete elem->second; }

    (*this).erase(elem);
    return true;
  }

  const auto &ref() const { return *this; }

  auto &ref() { return *this; }

  ~MeshPart() {
    if constexpr (std::is_pointer<T>::value) {
      for (auto item : *this) {
        if (item.second) delete item.second;
      }
    }
    (*this).clear();
  }
};

template<typename t, class MeshPart>
class meshpart {
  MeshPart _data;

  t begin() { return t(_data.Begin()); }

  t end() { return t(_data.End()); }

  template<typename P>
  t find(const P &p) {
    return t(_data.Find(p));
  }
};

using Vertices = MeshPart<Point, Vertex *>;

using Edges = MeshPart<Point, Edge *>;

using Faces = MeshPart<Point, Face>;

using BoundaryFaces = MeshPart<Point, int>;

Buffer &operator>>(Buffer &b, BoundaryFaces &B);

using Cells = MeshPart<Point, Cell *>;

class Mesh;

class CellBoundaryFaces {
  bnd_face bf[6];
  short N;
public:
  CellBoundaryFaces(const Mesh &M, const cell &c);

  short size() const { return N; }

  const bnd_face &operator[](int i) const { return bf[i]; }
};

Buffer &operator<<(Buffer &b, const CellBoundaryFaces &bf);

#endif // MESHPART_H
