#include "Row.hpp"

#include <algorithm>
#include <vector>

using std::vector;

std::ostream &operator<<(std::ostream &s, const RowEntry &e) {
  return s << "Entry: id= " << e.id << "; index= " << e.entryIndex;
}

void Row::InsertEntry(const Point &z) {
  if (entries.find(z) == entries.end()) entries[z] = RowEntry();
}

int Row::EntryId(const Point &z) const { return entries.find(z)->second.Id(); }

int Row::EntryIdChecked(const Point &z) const {
  auto entry = entries.find(z);
  if (entry == entries.end()) return -1;
  return entry->second.Id();
}

int Row::EntryIndex(const Point &z) const { return entries.find(z)->second.EntryIndex(); }

int Row::EntryIndexChecked(const Point &z) const {
  auto entry = entries.find(z);
  if (entry == entries.end()) return -1;
  return entry->second.EntryIndex();
}

template<typename S>
LogTextStream<S> &operator<<(LogTextStream<S> &s, const Row &R) {
  return s << "Row: id= " << R.id << "; index= " << R.index << "; dofs= " << R.dofs << endl
           << R.entries;
}

int Rows::RowId(const Point &z) const { return find(z)->second.Id(); }

int Rows::RowIdChecked(const Point &z) const {
  auto r = find(z);
  if (r == end()) return -1;
  return r->second.Id();
}

int Rows::RowIndex(const Point &z) const { return find(z)->second.Index(); }

int Rows::RowNumberOfDofs(const Point &z) const { return find(z)->second.NumberOfDofs(); }

int Rows::RowNumberOfEntries(const Point &z) const { return find(z)->second.NumberOfEntries(); }

void Rows::Insert(const Point &x, int n) {
  if (find(x) == end()) (*this)[x] = Row(n);
}

void Rows::InsertInfty(int numberOfDofs) {
  if (numberOfDofs == 0) return;
  Insert(Infty, numberOfDofs);
  auto r_infty = find(Infty);
  for (auto r = begin(); r != end(); ++r) {
    if (r != r_infty) AddEntry(r->first, Infty);
  }
}

void Rows::AddEntryFixedOrder(const Point &_row, const Point &_entry) {
  (*this)[_row].InsertEntry(_entry);
}

void Rows::AddEntry(const Point &x, const Point &y) {
  if (y < x) AddEntryFixedOrder(x, y);
  else AddEntryFixedOrder(y, x);
}

void Rows::InitIndices(bool SingleEntries) {
  using Iterator = std::unordered_map<Point, Row>::iterator;

  std::vector<Iterator> rowsVec{};
  rowsVec.reserve(size());
  for (Iterator rowIt = begin(); rowIt != end(); ++rowIt)
    rowsVec.push_back(rowIt);
  std::sort(rowsVec.begin(), rowsVec.end(),
            [](const Iterator &r0, const Iterator &r1) { return r0->first < r1->first; });
  for (size_t id = 0; id < rowsVec.size(); ++id)
    rowsVec[id]->second.Id(id);

  totalEntries = 0;
  totalDofs = 0;
  sumDofs = 0;
  for (Iterator &it : rowsVec) {
    Row &r0 = it->second;
    r0.Index(totalDofs);
    ++totalEntries;
    totalDofs += r0.NumberOfDofs() * r0.NumberOfDofs();
    sumDofs += r0.NumberOfDofs();
    double t; // TODO: uninitialized time!!!
    if (it->first.isTimeDep()) t = it->first.t();
    for (auto &[p1, e1] : r0.Entries()) {
      const auto iterator = find(p1);
      // Row corresponding to RowEntry point e1
      if constexpr (DebugLevel > 0) {
        if (iterator == end()) {
          std::stringstream s;
          s << "Point " << p1 << " could not be found in iterator!";
          THROW(s.str());
        }
      }
      const Row &r1 = iterator->second;
      e1.Id(r1.Id());
      e1.EntryIndex(totalDofs);
      ++totalEntries;
      if (SingleEntries && (p1.t() != t)) {
        totalDofs += r0.NumberOfDofs() * r1.NumberOfDofs();
      } else {
        totalDofs += 2 * r0.NumberOfDofs() * r1.NumberOfDofs();
      }
    }
  }
}

int Rows::Index(const RowIterator_const &r0, const RowIterator_const &r1) const {
  if (r0 == r1) /// index of the diagonal block corresponding to r0 (and r1)
    return r0->second.Index();
  if (r1->first < r0->first) /// index of the non-diagonal block corresponding to r0 and r1
    return r0->second.EntryIndex(r1->first);

  return r1->second.EntryIndex(r0->first) + r0->second.NumberOfDofs() * r1->second.NumberOfDofs();
  // TODO: Check, if this is faster
  //  if (r0 == r1) /// index of the diagonal block corresponding to r0 (and r1)
  //    return r0->second.Index();
  //  if (r1->second.Id() < r0->second.Id()) /// index of the non-diagonal block
  //  corresponding to r0 and r1
  //    return r0->second.EntryIndex(r1->first);
  //  return r1->second.EntryIndex(r0->first) + r0->second.NumberOfDofs() *
  //  r1->second.NumberOfDofs();
}

int Rows::IndexChecked(const RowIterator_const &r0, const RowIterator_const &r1) const {
  if (r0 == r1) return r0->second.Index();
  if (r1->first < r0->first) return r0->second.EntryIndexChecked(r1->first);
  int k = r1->second.EntryIndexChecked(r0->first);
  if (k == -1) return -1;
  return k + r0->second.NumberOfDofs() * r1->second.NumberOfDofs();
  // TODO: Check, if this is faster
  //  if (r0 == r1)
  //    return r0->second.Index();
  //  if (r1->second.Id() < r0->second.Id())
  //    return r0->second.EntryIndexChecked(r1->first);
  //  int k = r1->second.EntryIndexChecked(r0->first);
  //  if (k == -1) return -1;
  //  return k + r0->second.NumberOfDofs() * r1->second.NumberOfDofs();
}

int Rows::DoubleIndexChecked(const RowIterator_const &r0, const RowIterator_const &r1) const {
  if (r0->second.Id() == r1->second.Id()) return r0->second.Index();
  if (r1->second.Id() < r0->second.Id()) return r0->second.EntryIndexChecked(r1->first);
  else return r1->second.EntryIndexChecked(r0->first);
}
