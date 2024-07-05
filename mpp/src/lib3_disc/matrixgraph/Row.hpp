#ifndef ROW_HPP
#define ROW_HPP

#include "Point.hpp"

class RowEntry {
  /// number of the Entry with respect to all Entries (will be set in Rows)
  int id = -1;
  /// index later used as start index in the BasicVector for a Matrix
  int entryIndex = -1;
public:
  RowEntry() = default;

  void Id(int _id) { id = _id; }

  int Id() const { return id; }

  void EntryIndex(int index) { entryIndex = index; }

  int EntryIndex() const { return entryIndex; }

  friend std::ostream &operator<<(std::ostream &s, const RowEntry &entry);
};

class Row {
public:
  using RowEntryMap = std::unordered_map<Point, RowEntry>;
  using RowEntryIterator = std::unordered_map<Point, RowEntry>::iterator;
  using RowEntryIterator_const = std::unordered_map<Point, RowEntry>::const_iterator;
private:
  /// number of the entry with respect to all entries
  int id = -1;
  /// index later used as start index in the BasicVector for a Vector
  int index = -1;
  /// degrees of freedom in the Row, i.e., for each Entry
  int dofs = 0;
  /// map containing the entries
  RowEntryMap entries{};
public:
  Row() = default;

  explicit Row(int numberOfDofs) : dofs(numberOfDofs) {}

  void Id(int _id) { id = _id; }

  int Id() const { return id; }

  void Index(int _index) { index = _index; }

  int Index() const { return index; }

  void NumberOfDofs(int _dofs) { dofs = _dofs; }

  int NumberOfDofs() const { return dofs; }

  void InsertEntry(const Point &z);

  RowEntryMap &Entries() { return entries; }

  const RowEntryMap &Entries() const { return entries; }

  RowEntryIterator_const EntriesBegin() const { return entries.begin(); }

  RowEntryIterator_const EntriesEnd() const { return entries.end(); }

  RowEntryIterator_const FindEntry(const Point &z) const { return entries.find(z); }

  int NumberOfEntries() const { return entries.size(); }

  int EntryId(const Point &z) const;

  int EntryIdChecked(const Point &z) const;

  int EntryIndex(const Point &z) const;

  int EntryIndexChecked(const Point &z) const;

  template<typename S>
  friend LogTextStream<S> &operator<<(LogTextStream<S> &s, const Row &R);
};

class Rows : public std::unordered_map<Point, Row> {
public:
  using RowMap = std::unordered_map<Point, Row>;
  using RowIterator = std::unordered_map<Point, Row>::iterator;
  using RowIterator_const = std::unordered_map<Point, Row>::const_iterator;
private:
  /// total number of dofs (sum over all entries) -> size of Matrix
  size_t totalDofs = 0;
  /// total number of entries
  int totalEntries = 0;
  /// sum of number of dofs over all rows
  int sumDofs = 0;
public:
  Rows() = default;

  /// size of Matrix
  size_t NumberOfDofsTotal() const { return totalDofs; }

  int SumOfDoFs() const { return sumDofs; }

  int NumberOfEntriesTotal() const { return totalEntries; }

  void Insert(const Point &x, int n);

  void InsertInfty(int numberOfDofs);

  void AddEntryFixedOrder(const Point &_row, const Point &_entry);

  void AddEntry(const Point &x, const Point &y);

  /// initializes the ids and indices of the entries. Call after ALL entries are inserted
  void InitIndices(bool);

  int RowId(const Point &z) const;

  int RowIdChecked(const Point &z) const;

  int RowIndex(const Point &z) const;

  int RowNumberOfDofs(const Point &z) const;

  int RowNumberOfEntries(const Point &z) const;

  int Index(const RowIterator_const &r0, const RowIterator_const &r1) const;

  int IndexChecked(const RowIterator_const &r0, const RowIterator_const &r1) const;

  int DoubleIndexChecked(const RowIterator_const &r0, const RowIterator_const &r1) const;
};

#endif // ROW_HPP
