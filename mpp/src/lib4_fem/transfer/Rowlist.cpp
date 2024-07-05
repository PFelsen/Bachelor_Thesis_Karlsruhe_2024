#include "Rowlist.hpp"
#include <Vector.hpp>

RowList::RowList(const Vector &coarseVector, const Vector &fineVector) :
    coarseWeights(coarseVector.nR()), coarseRowIds(fineVector.nR()) {
  std::uninitialized_fill(coarseWeights.begin(), coarseWeights.end(), 0.0);
  std::uninitialized_fill(coarseRowIds.begin(), coarseRowIds.end(), std::vector<RowWeight>{});
}

RowList::RowList(RowList &&other) noexcept :
    coarseWeights(std::move(other.coarseWeights)), coarseRowIds(std::move(other.coarseRowIds)) {}

void RowList::AddRowWeight(int fineRow, int coarseRow, double coarseWeight) {
  coarseRowIds[fineRow].emplace_back(RowWeight{coarseRow, coarseWeight});
  coarseWeights[coarseRow] += coarseWeight;
}

// TODO: why only first entry is checked?
bool RowList::Empty() const { return coarseRowIds[0].empty(); }

bool RowList::Empty(int fineRowId) const { return coarseRowIds[fineRowId].empty(); }

size_t RowList::IdSize() const { return coarseRowIds.size(); }

size_t RowList::WeightSize() const { return coarseWeights.size(); }

const std::vector<RowWeight> &RowList::CoarseRows(int fineRowId) const {
  return coarseRowIds[fineRowId];
}

double RowList::CoarseWeight(int coarseRowId) const { return coarseWeights[coarseRowId]; }

double RowList::CalculateFineEntry(int fineId, const Vector &coarse, int j) const {
  const auto &coarseList = coarseRowIds[fineId];
  return std::accumulate(coarseList.begin(), coarseList.end(), 0.0,
                         [&coarse, &j](double sum, RowWeight r) {
                           return sum + r.Weight * coarse(r.Id, j);
                         });
}

std::ostream &operator<<(std::ostream &s, const RowList &rowlist) {
  for (int r = 0; r < rowlist.IdSize(); ++r) {
    auto item = rowlist.CoarseRows(r);
    s << "Id(" << r << ") = {";
    double w = 0.0;
    std::for_each(item.begin(), item.end(), [&w, &s](const RowWeight &i) {
      w += i.Weight;
      s << i.Id << ": " << i.Weight << " - ";
    });
    s << "Sum: " << w << "}" << endl;
  }
  return s;
}
