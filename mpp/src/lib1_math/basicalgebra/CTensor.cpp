#include "CTensor.hpp"
#include "CFunctions.hpp"
#include "Parallel.hpp"
#include "SaveLoad.hpp"

#include <cmath>

template<typename REAL>
RTensorT<REAL> CTensorT<REAL>::real() const {
  RTensorT<REAL> re(NFirstDimension, NSecondDimension, NThirdDimension);
  for (int i = 0; i < z.size(); ++i)
    re.Data()[i] = baReal(z[i]);
  return re;
}

template<typename REAL>
RTensorT<REAL> CTensorT<REAL>::imag() const {
  RTensorT<REAL> im(NFirstDimension, NSecondDimension, NThirdDimension);
  for (int i = 0; i < z.size(); ++i)
    im.Data()[i] = baImag(z[i]);
  return im;
}

template<typename REAL>
void CTensorT<REAL>::resize(int FirstComponents, int SecondComponents, int ThirdComponents) {
  z = std::vector<COMPLEX>(FirstComponents * SecondComponents * ThirdComponents, COMPLEX{});
  NFirstDimension = FirstComponents;
  NSecondDimension = SecondComponents;
  NThirdDimension = ThirdComponents;
}

template<typename REAL>
COMPLEX_TYPE<REAL> CTensorT<REAL>::Mean() const {
  COMPLEX mean{};
  for (auto const &value_vec : z)
    mean += value_vec;
  return mean / z.size();
}

template<typename REAL>
REAL CTensorT<REAL>::Variance() const {
  REAL var_real = real().Variance();
  REAL var_imag = imag().Variance();
  return var_real + var_imag;
}

template class CTensorT<double>;
