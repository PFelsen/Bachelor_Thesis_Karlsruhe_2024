#ifndef ROWLIST_HPP
#define ROWLIST_HPP

#include <cstddef>
#include <ostream>
#include <vector>


class Vector;

struct RowWeight {
  int Id;
  double Weight;
};

struct RowList {
  explicit RowList(const Vector &coarseVector, const Vector &fineVector);

  RowList(RowList &&other) noexcept;

  void AddRowWeight(int fineRow, int coarseRow, double coarseWeight);

  bool Empty() const;

  bool Empty(int fineRowId) const;

  size_t IdSize() const;

  size_t WeightSize() const;

  const std::vector<RowWeight> &CoarseRows(int fineRowId) const;

  double CoarseWeight(int coarseRowId) const;

  double CalculateFineEntry(int fineId, const Vector &coarse, int j = 0) const;
private:
  std::vector<double> coarseWeights{};
  std::vector<std::vector<RowWeight>> coarseRowIds{};
};

std::ostream &operator<<(std::ostream &s, const RowList &rowlist);

#endif // ROWLIST_HPP
