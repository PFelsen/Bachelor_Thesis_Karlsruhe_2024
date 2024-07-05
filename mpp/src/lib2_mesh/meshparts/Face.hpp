#ifndef _FACE_H_
#define _FACE_H_

#include <unordered_map>

#include "Point.hpp"

class Buffer;

// TODO Should be const ref, not possible yet bc of buffer
class Face {
  Point centerLeftCell;

  Point centerRightCell;
public:
  explicit Face(const Point &leftCell) : centerLeftCell(leftCell), centerRightCell(Infty) {}

  Face(const Point &leftCell, const Point &rightCell) :
      centerLeftCell((leftCell < rightCell) ? leftCell : rightCell),
      centerRightCell((leftCell < rightCell) ? rightCell : leftCell) {}

  Face() : centerLeftCell(Infty), centerRightCell(Infty) {}

  const Point &Left() const { return centerLeftCell; }

  const Point &Right() const { return centerRightCell; }

  friend std::ostream &operator<<(std::ostream &s, const Face &face);

  friend Buffer &operator<<(Buffer &b, const Face &face);

  friend Buffer &operator>>(Buffer &b, Face &face);
};

class face : public std::unordered_map<Point, Face>::const_iterator {
  typedef std::unordered_map<Point, Face>::const_iterator Iterator;
public:
  face() {}

  face(Iterator f) : Iterator(f) {}

  const Point &operator()() const { return (*this)->first; }

  const Point &Center() const { return (*this)->first; }

  const Face &operator*() const { return (*this)->second; }

  const Point &Left() const { return (*this)->second.Left(); }

  const Point &Right() const { return (*this)->second.Right(); }

  bool isEqual(const face &other) {
    return Center() == other.Center() && Left() == other.Left() && Right() == other.Right();
  }
};

inline std::ostream &operator<<(std::ostream &s, const face &F) {
  return s << F->first << " : " << F->second << endl;
}

class bnd_face : public std::unordered_map<Point, int>::const_iterator {
  typedef std::unordered_map<Point, int>::const_iterator Iterator;
public:
  bnd_face() = default;

  bnd_face(Iterator b) : Iterator(b) {}

  const Point &operator()() const { return (*this)->first; }

  int Part() const { return (*this)->second; }
};

inline std::ostream &operator<<(std::ostream &s, const bnd_face &b) {
  return s << b() << " : " << b.Part() << endl;
}

Buffer &operator<<(Buffer &b, const bnd_face &B);

#endif // of #ifndef _FACE_H_
