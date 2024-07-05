#ifndef IAEIGENVALUEBOUNDS_H
#define IAEIGENVALUEBOUNDS_H

#include "basicalgebra/RVector.hpp"


typedef RVector EigenvalueBounds;

class IALowerEigenvalueBounds : public EigenvalueBounds {
public:
  IALowerEigenvalueBounds() = default;

  IALowerEigenvalueBounds(int length);

  IALowerEigenvalueBounds &operator+=(const double &a);

  IALowerEigenvalueBounds &operator-=(const double &a);
};

class IAUpperEigenvalueBounds : public EigenvalueBounds {
public:
  IAUpperEigenvalueBounds() = default;

  IAUpperEigenvalueBounds(int length);

  IAUpperEigenvalueBounds &operator+=(const double &a);

  IAUpperEigenvalueBounds &operator-=(const double &a);
};

class IAEigenvalueEnclosures : public IARVector {
public:
  IAEigenvalueEnclosures() = default;

  IAEigenvalueEnclosures(int length);

  IAEigenvalueEnclosures(const IALowerEigenvalueBounds &, const IAUpperEigenvalueBounds &);

  void setLowerBounds(const IALowerEigenvalueBounds &);

  void setUpperBounds(const IAUpperEigenvalueBounds &);

  IALowerEigenvalueBounds getLowerBounds();

  IAUpperEigenvalueBounds getUpperBounds();

  void getLowerBounds(IALowerEigenvalueBounds &);

  void getUpperBounds(IAUpperEigenvalueBounds &);

  double getLowerBound(int n) { return inf((*this)[n]); }

  double getUpperBound(int n) { return sup((*this)[n]); }
};

#endif // IAEIGENVALUEBOUNDS_H
