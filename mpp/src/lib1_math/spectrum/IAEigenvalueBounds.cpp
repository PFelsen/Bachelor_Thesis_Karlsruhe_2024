#include "IAEigenvalueBounds.hpp"

IALowerEigenvalueBounds::IALowerEigenvalueBounds(int length) : EigenvalueBounds(length) {}

IALowerEigenvalueBounds &IALowerEigenvalueBounds::operator+=(const double &b) {
  for (int i = 0; i < size(); ++i)
    (*this)[i] = inf(IAInterval((*this)[i]) + b);
  return *this;
}

IALowerEigenvalueBounds &IALowerEigenvalueBounds::operator-=(const double &b) {
  for (int i = 0; i < size(); ++i)
    (*this)[i] = inf(IAInterval((*this)[i]) - b);
  return *this;
}

IAUpperEigenvalueBounds::IAUpperEigenvalueBounds(int length) : EigenvalueBounds(length) {}

IAUpperEigenvalueBounds &IAUpperEigenvalueBounds::operator+=(const double &b) {
  for (int i = 0; i < size(); ++i)
    (*this)[i] = sup(IAInterval((*this)[i]) + b);
  return *this;
}

IAUpperEigenvalueBounds &IAUpperEigenvalueBounds::operator-=(const double &b) {
  for (int i = 0; i < size(); ++i)
    (*this)[i] = sup(IAInterval((*this)[i]) - b);
  return *this;
}

IAEigenvalueEnclosures::IAEigenvalueEnclosures(int length) : IARVector(length) {}

IAEigenvalueEnclosures::IAEigenvalueEnclosures(const IALowerEigenvalueBounds &lambda,
                                               const IAUpperEigenvalueBounds &Lambda) :
    IARVector(std::min((int)lambda.size(), (int)Lambda.size())) {
  for (int i = 0; i < size(); ++i)
    (*this)[i] = IAInterval(lambda[i], Lambda[i]);
}

void IAEigenvalueEnclosures::setLowerBounds(const IALowerEigenvalueBounds &lambda) {
  if (lambda.size() != size()) THROW("Size does not fit!")
  for (int i = 0; i < size(); ++i)
    ((*this)[i]).setInf(lambda[i]);
}

void IAEigenvalueEnclosures::setUpperBounds(const IAUpperEigenvalueBounds &Lambda) {
  if (Lambda.size() != size()) THROW("Size does not fit!")
  for (int i = 0; i < size(); ++i)
    ((*this)[i]).setSup(Lambda[i]);
}

IALowerEigenvalueBounds IAEigenvalueEnclosures::getLowerBounds() {
  IALowerEigenvalueBounds lb(size());
  for (int i = 0; i < size(); ++i)
    lb[i] = inf((*this)[i]);
  return lb;
}

IAUpperEigenvalueBounds IAEigenvalueEnclosures::getUpperBounds() {
  IAUpperEigenvalueBounds ub(size());
  for (int i = 0; i < size(); ++i)
    ub[i] = sup((*this)[i]);
  return ub;
}

void IAEigenvalueEnclosures::getLowerBounds(IALowerEigenvalueBounds &lb) {
  for (int i = 0; i < size(); ++i)
    lb[i] = inf((*this)[i]);
}

void IAEigenvalueEnclosures::getUpperBounds(IAUpperEigenvalueBounds &ub) {
  for (int i = 0; i < size(); ++i)
    ub[i] = sup((*this)[i]);
}
