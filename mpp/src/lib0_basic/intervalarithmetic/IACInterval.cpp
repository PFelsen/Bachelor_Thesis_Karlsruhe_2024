#include "IACInterval.hpp"

#include <cidot.hpp>

IAInterval IACInterval::real() const { return IAInterval(cxsc::Re(CINTERVAL)); }

void IACInterval::setReal(const IAInterval &a) { cxsc::SetRe(CINTERVAL, a.INTERVAL); }

IAInterval IACInterval::imag() const { return IAInterval(cxsc::Im(CINTERVAL)); }

void IACInterval::setImag(const IAInterval &a) { cxsc::SetIm(CINTERVAL, a.INTERVAL); }

IACInterval IACInterval::conj() const { return {cxsc::conj(CINTERVAL)}; }

IAInterval IACInterval::absSqr() const {
  return IAInterval(cxsc::sqr(cxsc::Re(CINTERVAL)) + cxsc::sqr(cxsc::Im(CINTERVAL)));
}

IAInterval IACInterval::abs() const { return IAInterval(cxsc::abs(CINTERVAL)); }

std::complex<double> IACInterval::mid() const {
  return std::complex<double>(cxsc::_double(cxsc::mid(cxsc::Re(CINTERVAL))),
                              cxsc::_double(cxsc::mid(cxsc::Im(CINTERVAL))));
}

IACInterval operator-(const IACInterval &a) { return {-a.CINTERVAL}; }

IACInterval operator+(const IACInterval &a, const IACInterval &b) {
  return {a.CINTERVAL + b.CINTERVAL};
}

IACInterval operator+(const IACInterval &a, const std::complex<double> &b) {
  return {a.CINTERVAL + cxsc::complex(std::real(b), std::imag(b))};
}

IACInterval operator+(const std::complex<double> &a, const IACInterval &b) {
  return {cxsc::complex(std::real(a), std::imag(a)) + b.CINTERVAL};
}

IACInterval operator+(const IACInterval &a, const IAInterval &b) {
  return {a.CINTERVAL + IACInterval::extractInterval(b)};
}

IACInterval operator+(const IAInterval &a, const IACInterval &b) {
  return {IACInterval::extractInterval(a) + b.CINTERVAL};
}

IACInterval operator+(const IACInterval &a, const double &b) {
  return {a.CINTERVAL + cxsc::real(b)};
}

IACInterval operator+(const double &a, const IACInterval &b) {
  return {cxsc::real(a) + b.CINTERVAL};
}

IACInterval operator-(const IACInterval &a, const IACInterval &b) {
  return {a.CINTERVAL - b.CINTERVAL};
}

IACInterval operator-(const IACInterval &a, const std::complex<double> &b) {
  return {a.CINTERVAL - cxsc::complex(std::real(b), std::imag(b))};
}

IACInterval operator-(const std::complex<double> &a, const IACInterval &b) {
  return {cxsc::complex(std::real(a), std::imag(a)) - b.CINTERVAL};
}

IACInterval operator-(const IACInterval &a, const IAInterval &b) {
  return {a.CINTERVAL - IACInterval::extractInterval(b)};
}

IACInterval operator-(const IAInterval &a, const IACInterval &b) {
  return {IACInterval::extractInterval(a) - b.CINTERVAL};
}

IACInterval operator-(const IACInterval &a, const double &b) {
  return {a.CINTERVAL - cxsc::real(b)};
}

IACInterval operator-(const double &a, const IACInterval &b) {
  return {cxsc::real(a) - b.CINTERVAL};
}

IACInterval operator*(const IACInterval &a, const IACInterval &b) {
  return {a.CINTERVAL * b.CINTERVAL};
}

IACInterval operator*(const IACInterval &a, const std::complex<double> &b) {
  return {a.CINTERVAL * cxsc::complex(std::real(b), std::imag(b))};
}

IACInterval operator*(const std::complex<double> &a, const IACInterval &b) {
  return {cxsc::complex(std::real(a), std::imag(a)) * b.CINTERVAL};
}

IACInterval operator*(const IACInterval &a, const IAInterval &b) {
  return {a.CINTERVAL * IACInterval::extractInterval(b)};
}

IACInterval operator*(const IAInterval &a, const IACInterval &b) {
  return {IACInterval::extractInterval(a) * b.CINTERVAL};
}

IACInterval operator*(const IACInterval &a, const double &b) {
  return {a.CINTERVAL * cxsc::real(b)};
}

IACInterval operator*(const double &a, const IACInterval &b) {
  return {cxsc::real(a) * b.CINTERVAL};
}

IACInterval operator/(const IACInterval &a, const IACInterval &b) {
  return {a.CINTERVAL / b.CINTERVAL};
}

IACInterval operator/(const IACInterval &a, const std::complex<double> &b) {
  return {a.CINTERVAL / cxsc::complex(std::real(b), std::imag(b))};
}

IACInterval operator/(const std::complex<double> &a, const IACInterval &b) {
  return {cxsc::complex(std::real(a), std::imag(a)) / b.CINTERVAL};
}

IACInterval operator/(const IACInterval &a, const IAInterval &b) {
  return {a.CINTERVAL / IACInterval::extractInterval(b)};
}

IACInterval operator/(const IAInterval &a, const IACInterval &b) {
  return {IACInterval::extractInterval(a) / b.CINTERVAL};
}

IACInterval operator/(const IACInterval &a, const double &b) {
  return {a.CINTERVAL / cxsc::real(b)};
}

IACInterval operator/(const double &a, const IACInterval &b) {
  return {cxsc::real(a) / b.CINTERVAL};
}

bool operator==(const IACInterval &a, const IACInterval &b) { return a.CINTERVAL == b.CINTERVAL; }

bool operator!=(const IACInterval &a, const IACInterval &b) { return !(a == b); }

IACInterval accumulate(const std::vector<IACInterval> &a,
                       const std::vector<std::complex<double>> &b) {
  cxsc::cidotprecision sum(0.0);
  for (int i = 0; i < std::min(a.size(), b.size()); ++i)
    cxsc::accumulate(sum, a[i].CINTERVAL, cxsc::complex(std::real(b[i]), std::imag(b[i])));
  return cxsc::cinterval(sum);
}

IACInterval accumulate(const std::vector<std::complex<double>> &a,
                       const std::vector<IACInterval> &b) {
  cxsc::cidotprecision sum(0.0);
  for (int i = 0; i < std::min(a.size(), b.size()); ++i)
    cxsc::accumulate(sum, cxsc::complex(std::real(a[i]), std::imag(a[i])), b[i].CINTERVAL);
  return cxsc::cinterval(sum);
}

std::ostream &operator<<(std::ostream &os, const IACInterval &a) { return os << a.CINTERVAL; }