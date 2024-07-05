#ifndef LINEARTRANSFER_HPP
#define LINEARTRANSFER_HPP

#include "ITransfer.hpp"
#include "Rowlist.hpp"

namespace mpp {
std::vector<int> nodalIds(const Vector &vec, const vector<Point> &nodalPoints);

/// Calculates rowList entries for a single cell of the given type.
void linearTransfer(CELLTYPE celltype, const vector<Point> &coarseNodalPoints,
                    const vector<int> &coarseIds, double lvl, const Vector &fine, RowList &rowList);

RowList constructLinearRowList(const Vector &coarse, const Vector &fine, int coarseDegree,
                               int fineDegree);


} // namespace mpp

#endif // LINEARTRANSFER_HPP
