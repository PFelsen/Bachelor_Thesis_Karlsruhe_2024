#include "RTensor.hpp"
#include "GlobalDefinitions.hpp"
#include "Lapdef.hpp"
#include "Parallel.hpp"
#include "SaveLoad.hpp"

#include <cmath>

template<typename REAL>
void RTensorT<REAL>::resize(int FirstComponents, int SecondComponents, int ThirdComponents) {
  z = std::vector<REAL>(FirstComponents * SecondComponents * ThirdComponents, REAL{});
  NFirstDimension = FirstComponents;
  NSecondDimension = SecondComponents;
  NThirdDimension = ThirdComponents;
}

template<typename REAL>
REAL RTensorT<REAL>::Mean() const {
  REAL mean{};
  for (auto const &value : z)
    mean += value;
  return mean / REAL(z.size());
}

template<typename REAL>
REAL RTensorT<REAL>::Variance() const {
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

template class RTensorT<double>;
