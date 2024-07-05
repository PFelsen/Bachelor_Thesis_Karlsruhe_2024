#include "MixedLagrangeTransfer.hpp"
#include "LinearTransfer.hpp"
#include "Vector.hpp"

MixedLagrangeTransfer::MixedLagrangeTransfer(const Vector &coarse, const Vector &fine) :
    discs((const MixedLagrangeDiscretization &)coarse.GetDisc(),
          (const MixedLagrangeDiscretization &)fine.GetDisc()),
    rowList1(
        mpp::constructLinearRowList(coarse, fine, discs.first.Degree1(), discs.second.Degree1())),
    rowList2(
        mpp::constructLinearRowList(coarse, fine, discs.first.Degree2(), discs.second.Degree2())) {
  if (discs.first.Dim() != discs.second.Dim()) { THROW("Dimension of Discretizations do no match") }
  if (discs.first.Size() != discs.second.Size()) { THROW("Sizes of Discretizations do no match") }
}

void MixedLagrangeTransfer::Prolongate(const Vector &coarse, Vector &fine) const {
  if (&coarse.GetDisc() != &discs.first || &fine.GetDisc() != &discs.second) {
    THROW("Discretizations of transfer and vectors do not match")
  }

  int nDoF1 = discs.first.Dim() * discs.first.Size();

  fine = 0.0;
  for (int i = 0; i < rowList1.IdSize(); ++i) {
    int j0 = 0;
    if (!rowList1.Empty(i)) {
      for (int j = 0; j < nDoF1; ++j) {
        fine(i, j) = rowList1.CalculateFineEntry(i, coarse, j);
      }
      j0 = nDoF1;
    }
    if (!rowList2.Empty(i)) {
      for (int j = j0; j < fine.Dof(i); ++j) {
        fine(i, j) = rowList1.CalculateFineEntry(i, coarse, j);
      }
    }
  }
  fine.ClearDirichletValues();
}

void MixedLagrangeTransfer::Restrict(Vector &coarse, const Vector &fine) const {
  if (&coarse.GetDisc() != &discs.first || &fine.GetDisc() != &discs.second) {
    THROW("Discretizations of transfer and vectors do not match")
  }

  coarse = 0.0;
  //  for (int i = 0; i < rowList.IdSize(); ++i) {
  //    for (const auto &cRow: rowList.CoarseRows(i)) {
  //      for (int j = 0; j < fine.Dof(i); ++j) {
  //        coarse(cRow.Id, j) += cRow.Weight / rowList.CoarseWeight(cRow.Id) * fine(i, j);
  //      }
  //    }
  //  }
  coarse.ClearDirichletValues();
}

void MixedLagrangeTransfer::Project(Vector &coarse, const Vector &fine) const {
  if (&coarse.GetDisc() != &discs.first || &fine.GetDisc() != &discs.second) {
    THROW("Discretizations of transfer and vectors do not match")
  }

  coarse = 0.0;
  //  for (row r = coarse.rows(); r != coarse.rows_end(); ++r) {
  //    for (int j = 0; j < r.n(); ++j) {
  //      coarse(r, j) = fine(r(), j);
  //    }
  //  }
  coarse.ClearDirichletValues();
}
