#ifndef MIXEDLAGRANGETRANSFER_HPP
#define MIXEDLAGRANGETRANSFER_HPP

#include "ITransfer.hpp"
#include "MixedLagrangeDiscretization.hpp"
#include "Rowlist.hpp"

class MixedLagrangeTransfer {
  std::pair<const MixedLagrangeDiscretization &, const MixedLagrangeDiscretization &> discs;
  RowList rowList1;
  RowList rowList2;
public:
  MixedLagrangeTransfer(const Vector &coarse, const Vector &fine);

  void Prolongate(const Vector &coarse, Vector &fine) const;

  void Restrict(Vector &coarse, const Vector &fine) const;

  void Project(Vector &coarse, const Vector &fine) const;
};

#endif // MIXEDLAGRANGETRANSFER_HPP
