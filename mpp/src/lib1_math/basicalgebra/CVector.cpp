#include "CVector.hpp"
#include "CFunctions.hpp"
#include "Parallel.hpp"
#include "SaveLoad.hpp"

template<typename REAL>
CVectorT<REAL>::CVectorT(int length) : z(length, 0.0) {}

template<typename REAL>
CVectorT<REAL> &CVectorT<REAL>::conj() {
  for (int i = 0; i < z.size(); ++i)
    z[i] = baConj(z[i]);
  return *this;
}

template<typename REAL>
RVectorT<REAL> CVectorT<REAL>::real() const {
  RVectorT<REAL> re(z.size());
  for (int i = 0; i < z.size(); ++i)
    re[i] = baReal(z[i]);
  return re;
}

template<typename REAL>
RVectorT<REAL> CVectorT<REAL>::imag() const {
  RVectorT<REAL> im(z.size());
  for (int i = 0; i < z.size(); ++i)
    im[i] = baImag(z[i]);
  return im;
}

template<typename REAL>
REAL CVectorT<REAL>::normSqr() const {
  REAL sum = REAL{};
  for (int i = 0; i < z.size(); ++i)
    sum += baAbsSqr(z[i]);
  return sum;
}

template<typename REAL>
REAL CVectorT<REAL>::norm() const {
  return sqrt(normSqr());
}

template<typename REAL>
void CVectorT<REAL>::Accumulate(int commSplit) {
  PPM->SumOnCommSplit(z, commSplit);
}

template<typename REAL>
Saver &CVectorT<REAL>::save(Saver &saver) const {
  saver << size();
  for (int i = 0; i < z.size(); ++i)
    saver << z[i];
  return saver;
}

template<typename REAL>
Loader &CVectorT<REAL>::load(Loader &loader) {
  int N;
  loader >> N;
  resize(N);
  for (int i = 0; i < z.size(); ++i)
    loader >> z[i];
  return loader;
}

template<typename REAL>
COMPLEX_TYPE<REAL> CVectorT<REAL>::Mean() const {
  COMPLEX_TYPE<REAL> mean{};
  for (auto const &value_vec : z)
    mean += value_vec;
  return mean / z.size();
}

template<typename REAL>
REAL CVectorT<REAL>::Variance() const {
  REAL var_real = real().Variance();
  REAL var_imag = imag().Variance();
  return var_real + var_imag;
}

template class CVectorT<double>;

#ifdef BUILD_IA

template class CVectorT<IAInterval>;

CVector mid(const IACVector &IA_v) {
  CVector v(IA_v.size());
  for (int i = 0; i < v.size(); ++i)
    v[i] = mid(IA_v[i]);
  return v;
}

#endif // BUILD_IA