#include "LagrangeTransfer.hpp"
#include "LinearTransfer.hpp"
#include "Vector.hpp"

LagrangeTransfer::LagrangeTransfer(const Vector &coarse, const Vector &fine) :
    discs((const LagrangeDiscretization &)coarse.GetDisc(),
          (const LagrangeDiscretization &)fine.GetDisc()),
    rowList(
        mpp::constructLinearRowList(coarse, fine, discs.first.Degree(), discs.second.Degree())) {}

void LagrangeTransfer::Prolongate(const Vector &coarse, Vector &fine) const {
  if (&coarse.GetDisc() != &discs.first || &fine.GetDisc() != &discs.second) {
    THROW("Discretizations of transfer and vectors do not match")
  }

  fine = 0.0;
  for (int i = 0; i < rowList.IdSize(); ++i) {
    for (int j = 0; j < fine.Dof(i); ++j) {
      fine(i, j) = rowList.CalculateFineEntry(i, coarse, j);
    }
  }
  fine.ClearDirichletValues();
}

void LagrangeTransfer::ProlongateTransposed(Vector &coarse, const Vector &fine) const {
  if (&coarse.GetDisc() != &discs.first || &fine.GetDisc() != &discs.second) {
    THROW("Discretizations of transfer and vectors do not match")
  }

  coarse = 0.0;
  for (int i = 0; i < rowList.IdSize(); ++i) {
    for (const auto &cRow : rowList.CoarseRows(i)) {
      for (int j = 0; j < fine.Dof(i); ++j) {
        coarse(cRow.Id, j) += cRow.Weight * fine(i, j);
      }
    }
  }
  coarse.ClearDirichletValues();
  coarse.Collect();
}

void LagrangeTransfer::Restrict(Vector &coarse, const Vector &fine) const {
  if (&coarse.GetDisc() != &discs.first || &fine.GetDisc() != &discs.second) {
    THROW("Discretizations of transfer and vectors do not match")
  }

  coarse = 0.0;
  for (int i = 0; i < rowList.IdSize(); ++i) {
    for (const auto &cRow : rowList.CoarseRows(i)) {
      for (int j = 0; j < fine.Dof(i); ++j) {
        coarse(cRow.Id, j) += cRow.Weight / rowList.CoarseWeight(cRow.Id) * fine(i, j);
      }
    }
  }
  coarse.ClearDirichletValues();
}

void LagrangeTransfer::Project(Vector &coarse, const Vector &fine) const {
  if (&coarse.GetDisc() != &discs.first || &fine.GetDisc() != &discs.second) {
    THROW("Discretizations of transfer and vectors do not match")
  }

  coarse = 0.0;
  for (row r = coarse.rows(); r != coarse.rows_end(); ++r) {
    for (int j = 0; j < r.n(); ++j) {
      coarse(r, j) = fine(r(), j);
    }
  }
  coarse.ClearDirichletValues();
}

RMatrix LagrangeTransfer::AsMatrix() const {
  RMatrix transferMatrix((int)rowList.IdSize(), (int)rowList.WeightSize());
  for (int i = 0; i < rowList.IdSize(); ++i) {
    for (auto r : rowList.CoarseRows(i)) {
      transferMatrix(i, r.Id) = r.Weight;
    }
  }
  return transferMatrix;
}
