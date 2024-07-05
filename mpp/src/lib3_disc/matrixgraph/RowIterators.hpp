#ifndef ROWITERATORS_HPP
#define ROWITERATORS_HPP

#include "Row.hpp"

class entry : public std::unordered_map<Point, RowEntry>::const_iterator {
public:
  entry() {}

  entry(std::unordered_map<Point, RowEntry>::const_iterator e) :
      std::unordered_map<Point, RowEntry>::const_iterator(e) {}

  const Point &operator()() const { return (*this)->first; }

  int Id() const { return (*this)->second.Id(); }

  int EntryIndex() const { return (*this)->second.EntryIndex(); }

  friend std::ostream &operator<<(std::ostream &s, const entry &e) {
    return s << e->first << " | " << e->second;
  }

  // ================================
  // TODO: Use methods above instead
  // ================================
  int GetEntry() const { return EntryIndex(); }
};

class row : public std::unordered_map<Point, Row>::const_iterator {
public:
  row() {}

  row(std::unordered_map<Point, Row>::const_iterator r) :
      std::unordered_map<Point, Row>::const_iterator(r) {}

  const Point &operator()() const { return (*this)->first; }

  int NumberOfEntries() const { return (*this)->second.NumberOfEntries(); }

  int Id() const { return (*this)->second.Id(); }

  int Index() const { return (*this)->second.Index(); }

  int NumberOfDofs() const { return (*this)->second.NumberOfDofs(); }

  int EntryIndex(const Point &z) const { return (*this)->second.EntryIndex(z); }

  int EntryIndexChecked(const Point &z) const { return (*this)->second.EntryIndexChecked(z); }

  template<typename S>
  friend LogTextStream<S> &operator<<(LogTextStream<S> &s, const row &r) {
    return s << r->first << " | " << r->second;
  }

  friend bool operator<(const row &r0, const row &r1) { return (r0() < r1()); }

  // ================================
  // TODO: Use methods above instead
  // ================================
  entry entries() const { return entry((*this)->second.EntriesBegin()); }

  entry entries_end() const { return entry((*this)->second.EntriesEnd()); }

  entry find_entry(const Point &z) const { return entry((*this)->second.FindEntry(z)); }

  int n() const { return (*this)->second.NumberOfDofs(); }

  int size() const { return (*this)->second.NumberOfEntries(); }

  int GetEntry() const { return Index(); }

  int GetEntry(const Point &z) const { return EntryIndex(z); }

  int GetEntryX(const Point &z) const { return EntryIndexChecked(z); }
};

#endif // ROWITERATORS_HPP
