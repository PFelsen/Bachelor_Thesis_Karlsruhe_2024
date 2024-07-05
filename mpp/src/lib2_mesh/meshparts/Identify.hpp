#ifndef _IDENTIFY_H_
#define _IDENTIFY_H_

#include "MeshPart.hpp"

class IdentifySet : public std::vector<Point> {
public:
  IdentifySet(){};

  void Add(const Point &x);

  friend std::ostream &operator<<(std::ostream &s, const IdentifySet &I);
};

class identifyset : public std::unordered_map<Point, IdentifySet>::const_iterator {
  typedef std::unordered_map<Point, IdentifySet>::const_iterator Iterator;
public:
  identifyset(const Iterator &I) : Iterator(I) {}

  const Point &operator()() const { return (*this)->first; }

  const IdentifySet &operator*() const { return (*this)->second; }

  long unsigned int size() const { return (*this)->second.size(); }

  Point operator[](int i) const { return (*this)->second[i]; }

  bool master() const { return ((*this)->first < (*this)->second[0]); }

  friend std::ostream &operator<<(std::ostream &s, const identifyset &I);
};

class IdentifySets : public MeshPart<Point, IdentifySet> {
public:
  identifyset identifysets() const { return identifyset(begin()); }

  identifyset identifysets_end() const { return identifyset(end()); }

  identifyset find_identifyset(const Point &z) const { return identifyset(find(z)); }

  void Insert(const Point &x, const Point &y);

  void Insert(const identifyset &is);

  void Append(const Point &x, const identifyset &is);

  void Identify(const Point &y, int mode);
};

#endif // of #ifndef _IDENTIFY_H_
