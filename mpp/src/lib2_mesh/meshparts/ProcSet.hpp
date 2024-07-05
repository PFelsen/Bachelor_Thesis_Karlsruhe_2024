#ifndef _PROCSET_H_
#define _PROCSET_H_

#include "MeshPart.hpp"

typedef short intproc;  // number of processors
typedef short intproci; // number of processors
typedef short intprocs;

class ProcSet : public std::vector<intproc> {
  int commSplit = 0;
public:
  ProcSet() {}

  ProcSet(int commSplit) : commSplit(commSplit) {}

  ProcSet(intproc q, int commSplit) : commSplit(commSplit), std::vector<intproc>(1) {
    (*this)[0] = q;
  }

  ProcSet(intproc p, intproc q, int commSplit) : commSplit(commSplit), std::vector<intproc>(2) {
    (*this)[0] = p;
    (*this)[1] = q;
  }

  ProcSet(const ProcSet &P) : commSplit(P.CommSplit()), std::vector<intproc>(0) {
    for (long unsigned int i = 0; i < P.size(); ++i)
      Add(P[i]);
  }

  int CommSplit() const { return commSplit; }

  void Add(intproc q);

  void Append(intproc q);

  void Append(const ProcSet &PS);

  void Add(const ProcSet &PS);

  bool Erase();

  void erase(intproc q);

  intproc master() const;

  void SetMaster(intproc q);

  bool const equalProcSet(const ProcSet &P) const;

  bool existselementof(const ProcSet &Q) const;

  bool existselementof(const intproc &q) const;

  bool subset(const ProcSet &Q) const;

  void Sort();

  void InsertFront(const ProcSet &PS);

  std::string ToString() const;

  friend Buffer &operator<<(Buffer &b, const ProcSet &P);

  friend Buffer &operator>>(Buffer &b, ProcSet &P);
};

inline std::ostream &operator<<(std::ostream &s, const ProcSet &P) {
  for (long unsigned int i = 0; i < P.size(); ++i)
    s << " " << P[i];
  return s;
}

class procset : public std::unordered_map<Point, ProcSet>::const_iterator {
  typedef std::unordered_map<Point, ProcSet>::const_iterator Iterator;
public:
  procset(Iterator p) : Iterator(p) {}

  const Point &operator()() const { return (*this)->first; }

  const ProcSet &operator*() const { return (*this)->second; }

  long unsigned int size() const { return (*this)->second.size(); }

  intproc operator[](int i) const { return (*this)->second[i]; }

  intproc master() const { return (*this)->second[0]; }

  bool in(intproc q) const {
    for (long unsigned int i = 0; i < size(); ++i)
      if ((*this)->second[i] == q) return true;
    return false;
  }
};

inline Buffer &operator<<(Buffer &buffer, const procset &P) { return buffer << P() << *P; }

inline Buffer &operator>>(Buffer &buffer, std::pair<Point, ProcSet> &ps) {
  return buffer >> ps.first >> ps.second;
}

inline std::ostream &operator<<(std::ostream &s, const procset &p) {
  return s << p() << " : " << *p;
}

class ProcSets : public MeshPart<Point, ProcSet> {
  int commSplit = 0;
public:
  ProcSets() {}

  procset procsets() const { return {begin()}; }

  procset procsets_end() const { return {end()}; }

  procset find_procset(const Point &z) const { return procset(find(z)); }

  template<class C>
  procset find_procset(const C &c) const {
    return procset(find(c()));
  }

  void SetCommSplit(int _commSplit) { commSplit = _commSplit; }

  void Copy(const procset &p, const Point &z);

  void Copy(const procset &p);

  void Add(const Point &z, intproc q);

  void Add(const Point &z, const procset &q);

  void AddInfty();

  void Append(const Point &z, intproc q);

  void Append(const Point &z, const procset &q);

  void Insert(const Point &z, const ProcSet &PS);

  void InsertFront(const Point &z, const ProcSet &PS);

  void Replace(const Point &z, const ProcSet &PS);

  bool on(std::unordered_map<Point, ProcSet>::iterator p, intproc q);

  void RemoveIfOn(const Point &z, intproc q);

  void RemoveSingle();

  void Clean();

  bool master(const Point &z) const;

  void CheckConsistency() const;
};

#endif // of #ifndef _PROCSET_H_
