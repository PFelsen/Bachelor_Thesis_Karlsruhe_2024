#include "Identify.hpp"
#include "Parallel.hpp"

void IdentifySet::Add(const Point &x) {
  long unsigned int m = size();
  if (m == 0) {
    resize(1);
    (*this)[0] = x;
    return;
  }
  for (short i = 0; i < m; ++i)
    if ((*this)[i] == x) return;
  resize(m + 1);
  if (x < (*this)[0]) {
    (*this)[m] = (*this)[0];
    (*this)[0] = x;
  } else (*this)[m] = x;
}

std::ostream &operator<<(std::ostream &s, const IdentifySet &I) {
  for (short i = 0; i < I.size(); ++i)
    s << " " << I[i];
  return s;
}

std::ostream &operator<<(std::ostream &s, const identifyset &I) {
  return s << I() << " : " << *I << endl;
}

/*
 * IdentifySets
 */

void IdentifySets::Insert(const Point &x, const Point &y) {
  (*this)[x];
  auto i = find(x);
  if (i != end()) i->second.Add(y);
}

void IdentifySets::Insert(const identifyset &is) { Append(is(), is); }

void IdentifySets::Append(const Point &x, const identifyset &is) {
  (*this)[x];
  auto i = find(x);
  if (is() != x) i->second.Add(is());
  for (int j = 0; j < is.size(); ++j)
    if (is[j] != x) i->second.Add(is[j]);
}

void IdentifySets::Identify(const Point &y, int mode) {
  if (find(y) != end()) return;
  if (mode == 999) {
    Point X1(1, 0);
    Point X2(0, 1);
    if (y[0] < GeometricTolerance) Insert(y, y + X1);
    if (abs(y[0] - 1) < GeometricTolerance) Insert(y, y - X1);
    if (y[1] < GeometricTolerance) Insert(y, y + X2);
    if (abs(y[1] - 1) < GeometricTolerance) Insert(y, y - X2);
    if (y[0] < GeometricTolerance && y[1] < GeometricTolerance) Insert(y, y + X1 + X2);
    if (abs(y[0] - 1) < GeometricTolerance && y[1] < GeometricTolerance) Insert(y, y - X1 + X2);
    if (y[0] < GeometricTolerance && abs(y[1] - 1) < GeometricTolerance) Insert(y, y + X1 - X2);
    if (abs(y[0] - 1) < GeometricTolerance && abs(y[1] - 1) < GeometricTolerance)
      Insert(y, y - X1 - X2);
  } else if (mode == 9999) {
    Point s[3];
    int f[3] = {0, 0, 0};

    s[0] = Point(1, 0, 0);
    s[1] = Point(0, 1, 0);
    s[2] = Point(0, 0, 1);
    for (int c = 0; c < y.SpaceDim(); ++c)
      for (int t = 0; t < 2; ++t)
        if (abs(y[c] - t) < GeometricTolerance) {
          f[c] = (t == 0) ? 1 : -1;
          Insert(y, y + f[c] * s[c]);
        };
    for (int c1 = 1; c1 < 3; ++c1)
      for (int c2 = 0; c2 < c1; ++c2)
        if ((f[c1] != 0) && (f[c2] != 0)) Insert(y, y + f[c1] * s[c1] + f[c2] * s[c2]);
    if ((f[0] != 0) && (f[1] != 0) && (f[2] != 0))
      Insert(y, y + f[0] * s[0] + f[1] * s[1] + f[2] * s[2]);
  }
}
