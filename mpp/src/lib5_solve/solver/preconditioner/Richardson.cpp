#include "Richardson.hpp"
#include "Vector.hpp"

void Richardson::Construct(const Matrix &A) {}

void Richardson::Destruct() {}

void Richardson::multiply(Vector &u, const Vector &b) const {
  u = damp * b;
  u.SetAccumulateFlag(false);
  u.Accumulate();
}