#include "RVector.hpp"
#include "Parallel.hpp"
#include "SaveLoad.hpp"

namespace mpp_ba {
double tolerance = 1e-12;

void SetTolerance(double tol) { tolerance = tol; }

double GetTolerance() { return tolerance; }

bool isNear(double a, double b) { return abs(a - b) < tolerance; }

bool isNearZero(double a) { return abs(a) < tolerance; }

#ifdef BUILD_IA

bool isNear(const IAInterval &a, const IAInterval &b) {
  return (abs(inf(a) - inf(b)) < tolerance) && (abs(sup(a) - sup(b)) < tolerance);
}

bool isNearZero(const IAInterval &a) {
  return (abs(inf(a)) < tolerance) && (abs(sup(a)) < tolerance);
}

#endif
} // namespace mpp_ba

template<typename REAL>
REAL RVectorT<REAL>::normSqr() const {
  REAL sum{};
  for (int i = 0; i < size(); ++i)
    sum += z[i] * z[i];
  return sum;
}

template<typename REAL>
REAL RVectorT<REAL>::norm() const {
  return sqrt(normSqr());
}

template<typename REAL>
REAL RVectorT<REAL>::normAccumulated(int commSplit) const {
  return sqrt(PPM->SumOnCommSplit(normSqr(), commSplit));
}

template<typename REAL>
REAL RVectorT<REAL>::Min() const {
  if constexpr (std::is_same_v<REAL, double>) {
    double s = infty;
    for (int i = 0; i < z.size(); ++i)
      s = std::min(s, z[i]);
    return s;
  } else THROW("Min not implemented")
}

template<typename REAL>
REAL RVectorT<REAL>::MinAccumulated(int commSplit) const {
  if constexpr (std::is_same_v<REAL, double>) {
    return PPM->Min(Min(), commSplit);
  } else THROW("Min not implemented")
}

template<typename REAL>
REAL RVectorT<REAL>::Max() const {
  if constexpr (std::is_same_v<REAL, double>) {
    double s = -infty;
    for (int i = 0; i < z.size(); ++i)
      s = std::max(s, z[i]);
    return s;
  } else THROW("Max not implemented")
}

template<typename REAL>
REAL RVectorT<REAL>::MaxAccumulated(int commSplit) const {
  if constexpr (std::is_same_v<REAL, double>) {
    return PPM->Max(Max(), commSplit);
  } else THROW("Max not implemented")
}

template<typename REAL>
void RVectorT<REAL>::Accumulate(int commSplit) {
  PPM->SumOnCommSplit(z, commSplit);
}

template<typename REAL>
void RVectorT<REAL>::CleanZero() {
  for (int i = 0; i < z.size(); ++i)
    if (mpp_ba::isNear(z[i], REAL{})) z[i] = REAL{};
}

template<typename REAL>
Saver &RVectorT<REAL>::save(Saver &saver) const {
  saver << size();
  for (int i = 0; i < z.size(); ++i)
    saver << z[i];
  return saver;
}

template<typename REAL>
Loader &RVectorT<REAL>::load(Loader &loader) {
  int N;
  loader >> N;
  resize(N);
  for (int i = 0; i < z.size(); ++i)
    loader >> z[i];
  return loader;
}

template<typename REAL>
REAL RVectorT<REAL>::Variance() const {
  REAL m2{};
  REAL diff{};
  REAL correction{};
  REAL mean = this->Mean();
  for (auto const &value : z) {
    diff = (value - mean);
    m2 += diff * diff;
    correction += diff;
  }
  return (m2 - correction * correction / z.size()) / (z.size() - 1);
}

template<typename REAL>
REAL RVectorT<REAL>::Mean() const {
  REAL mean{};
  for (auto const &value : z)
    mean += value;
  return mean / REAL(z.size());
}

template<typename REAL>
REAL RVectorT<REAL>::normOnes() const {
  REAL sum{};
  for (auto const &value : z)
    sum += value;
  return sum;
}

template class RVectorT<double>;

#ifdef BUILD_IA

template class RVectorT<IAInterval>;

RVector mid(const IARVector &IA_v) {
  RVector v(IA_v.size());
  for (int i = 0; i < v.size(); ++i)
    v[i] = mid(IA_v[i]);
  return v;
}

RVector sup(const IARVector &IA_v) {
  RVector v(IA_v.size());
  for (int i = 0; i < v.size(); ++i)
    v[i] = sup(IA_v[i]);
  return v;
}

RVector inf(const IARVector &IA_v) {
  RVector v(IA_v.size());
  for (int i = 0; i < v.size(); ++i)
    v[i] = inf(IA_v[i]);
  return v;
}

#endif // BUILD_IA