#ifndef LAGRANGETRANSFER_HPP
#define LAGRANGETRANSFER_HPP

#include "ITransfer.hpp"
#include "LagrangeDiscretization.hpp"
#include "Rowlist.hpp"

class LagrangeTransfer {
  std::pair<const LagrangeDiscretization &, const LagrangeDiscretization &> discs;
  RowList rowList;
public:
  LagrangeTransfer(const Vector &coarse, const Vector &fine);

  void Prolongate(const Vector &coarse, Vector &fine) const;

  void ProlongateTransposed(Vector &coarse, const Vector &fine) const;

  void Restrict(Vector &coarse, const Vector &fine) const;

  // TODO: Currently only works if fine vector contains all coarse rows.
  void Project(Vector &coarse, const Vector &fine) const;

  RMatrix AsMatrix() const;
};
#endif // LAGRANGETRANSFER_HPP
