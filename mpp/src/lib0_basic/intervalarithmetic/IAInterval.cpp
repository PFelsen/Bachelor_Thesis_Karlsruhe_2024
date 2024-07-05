#include "IAInterval.hpp"

#include <fi_lib.hpp>
#include <idot.hpp>
#include <imath.hpp>

void operator>>(std::string &input, IAInterval &a) { input >> a.INTERVAL; }

void operator>>(const char *input, IAInterval &a) {
  std::string input_s(input);
  input_s >> a;
}

IAInterval &IAInterval::operator+=(const IAInterval &a) {
  *this = *this + a;
  return *this;
}

IAInterval &IAInterval::operator+=(const double &a) {
  *this = *this + a;
  return *this;
}

IAInterval &IAInterval::operator-=(const IAInterval &a) {
  *this = *this - a;
  return *this;
}

IAInterval &IAInterval::operator-=(const double &a) {
  *this = *this - a;
  return *this;
}

IAInterval &IAInterval::operator*=(const IAInterval &a) {
  *this = *this * a;
  return *this;
}

IAInterval &IAInterval::operator*=(const double &a) {
  *this = *this * a;
  return *this;
}

IAInterval &IAInterval::operator/=(const IAInterval &a) {
  *this = *this / a;
  return *this;
}

IAInterval &IAInterval::operator/=(const double &a) {
  *this = *this / a;
  return *this;
}

double IAInterval::inf() const { return cxsc::_double(cxsc::Inf(INTERVAL)); }

double IAInterval::sup() const { return cxsc::_double(cxsc::Sup(INTERVAL)); }

void IAInterval::setInf(double inf) { cxsc::SetInf(INTERVAL, inf); }

void IAInterval::setSup(double sup) { cxsc::SetSup(INTERVAL, sup); }

IAInterval IAInterval::abs() const { return IAInterval(cxsc::abs(INTERVAL)); }

double IAInterval::mid() const { return cxsc::_double(cxsc::mid(INTERVAL)); }

double IAInterval::diam() const { return cxsc::_double(cxsc::diam(INTERVAL)); }

IAInterval operator-(const IAInterval &a) { return IAInterval(-a.INTERVAL); }

IAInterval operator+(const IAInterval &a, const IAInterval &b) {
  return IAInterval(a.INTERVAL + b.INTERVAL);
}

IAInterval operator+(const double &a, const IAInterval &b) { return IAInterval(a + b.INTERVAL); }

IAInterval operator+(const IAInterval &a, const double &b) { return IAInterval(a.INTERVAL + b); }

IAInterval operator-(const IAInterval &a, const IAInterval &b) {
  return IAInterval(a.INTERVAL - b.INTERVAL);
}

IAInterval operator-(const double &a, const IAInterval &b) { return IAInterval(a - b.INTERVAL); }

IAInterval operator-(const IAInterval &a, const double &b) { return IAInterval(a.INTERVAL - b); }

IAInterval operator*(const IAInterval &a, const IAInterval &b) {
  if (&a == &b) return sqr(a);
  return IAInterval(a.INTERVAL * b.INTERVAL);
}

IAInterval operator*(const double &a, const IAInterval &b) { return IAInterval(a * b.INTERVAL); }

IAInterval operator*(const IAInterval &a, const double &b) { return IAInterval(a.INTERVAL * b); }

IAInterval operator/(const IAInterval &a, const IAInterval &b) {
  return IAInterval(a.INTERVAL / b.INTERVAL);
}

IAInterval operator/(const double &a, const IAInterval &b) { return IAInterval(a / b.INTERVAL); }

IAInterval operator/(const IAInterval &a, const double &b) { return IAInterval(a.INTERVAL / b); }

bool operator==(const IAInterval &a, const IAInterval &b) {
  return (cxsc::Inf(a.INTERVAL) == cxsc::Inf(b.INTERVAL))
         && (cxsc::Sup(a.INTERVAL) == cxsc::Sup(b.INTERVAL));
}

bool operator!=(const IAInterval &a, const IAInterval &b) { return !(a == b); }

bool operator==(const IAInterval &a, const double &b) {
  return (cxsc::Inf(a.INTERVAL) == b) && (cxsc::Sup(a.INTERVAL) == b);
}

bool operator!=(const IAInterval &a, const double &b) { return !(a == b); }

bool operator==(const double &a, const IAInterval &b) { return b == a; }

bool operator!=(const double &a, const IAInterval &b) { return b != a; }

IAInterval operator|(const IAInterval &a, const IAInterval &b) {
  return IAInterval(a.INTERVAL | b.INTERVAL);
}

IAInterval hull(const IAInterval &a, const IAInterval &b) {
  return IAInterval(a.INTERVAL | b.INTERVAL);
}

IAInterval operator&(const IAInterval &a, const IAInterval &b) {
  return IAInterval(a.INTERVAL & b.INTERVAL);
}

IAInterval intsec(const IAInterval &a, const IAInterval &b) {
  return IAInterval(a.INTERVAL & b.INTERVAL);
}

bool disjoint(const IAInterval &a, const IAInterval &b) {
  return cxsc::dis_ii(a.INTERVAL, b.INTERVAL);
}

bool operator<=(const double &a, const IAInterval &b) { return a <= b.INTERVAL; }

bool in(const double &a, const IAInterval &b) { return cxsc::in(a, b.INTERVAL); }

bool operator<=(const IAInterval &a, const IAInterval &b) { return a.INTERVAL <= b.INTERVAL; }

bool in(const IAInterval &a, const IAInterval &b) { return cxsc::in(a.INTERVAL, b.INTERVAL); }

bool operator<(const IAInterval &a, const IAInterval &b) { return a.INTERVAL < b.INTERVAL; }

bool operator<(const double &a, const IAInterval &b) { return a < b.INTERVAL; }

bool operator>=(const IAInterval &a, const double &b) { return a.INTERVAL >= b; }

bool operator>=(const IAInterval &a, const IAInterval &b) { return a.INTERVAL >= b.INTERVAL; }

bool operator>(const IAInterval &a, const double &b) { return a.INTERVAL > b; }

bool operator>(const IAInterval &a, const IAInterval &b) { return a.INTERVAL > b.INTERVAL; }

IAInterval pow(const IAInterval &a, const IAInterval &b) {
  if (inf(a) <= 0.0) THROW("Non positive exponent error")
  return IAInterval(cxsc::pow(a.INTERVAL, b.INTERVAL));
}

IAInterval pow(const IAInterval &a, int n) { return IAInterval(cxsc::power(a.INTERVAL, n)); }

std::ostream &operator<<(std::ostream &os, const IAInterval &a) { return os << a.INTERVAL; }

IAInterval accumulate(const std::vector<IAInterval> &a, const std::vector<IAInterval> &b) {
  cxsc::idotprecision sum_inf(0.0);
  cxsc::idotprecision sum_sup(0.0);
  for (int i = 0; i < std::min(a.size(), b.size()); ++i) {
    cxsc::accumulate(sum_inf, cxsc::Inf(a[i].INTERVAL), b[i].INTERVAL);
    cxsc::accumulate(sum_sup, cxsc::Sup(a[i].INTERVAL), b[i].INTERVAL);
  }
  return IAInterval(cxsc::hull(cxsc::interval(sum_inf), cxsc::interval(sum_sup)));
}

IAInterval accumulate(const std::vector<IAInterval> &a, const std::vector<double> &b) {
  cxsc::idotprecision sum(0.0);
  for (int i = 0; i < std::min(a.size(), b.size()); ++i)
    cxsc::accumulate(sum, a[i].INTERVAL, cxsc::real(b[i]));
  return IAInterval(cxsc::interval(sum));
}

IAInterval accumulate(const std::vector<double> &a, const std::vector<IAInterval> &b) {
  cxsc::idotprecision sum(0.0);
  for (int i = 0; i < std::min(a.size(), b.size()); ++i)
    cxsc::accumulate(sum, cxsc::real(a[i]), b[i].INTERVAL);
  return IAInterval(cxsc::interval(sum));
}

IAInterval accumulate(const std::vector<double> &a, const std::vector<double> &b) {
  cxsc::idotprecision sum(0.0);
  for (int i = 0; i < std::min(a.size(), b.size()); ++i)
    cxsc::accumulate(sum, cxsc::real(a[i]), cxsc::real(b[i]));
  return IAInterval(cxsc::interval(sum));
}

IAInterval exp(const IAInterval &a) { return IAInterval(fi_lib::j_exp(a.INTERVAL)); }

IAInterval sinh(const IAInterval &a) { return IAInterval(fi_lib::j_sinh(a.INTERVAL)); }

IAInterval cosh(const IAInterval &a) { return IAInterval(fi_lib::j_cosh(a.INTERVAL)); }

IAInterval coth(const IAInterval &a) { return IAInterval(fi_lib::j_coth(a.INTERVAL)); }

IAInterval tanh(const IAInterval &a) { return IAInterval(fi_lib::j_tanh(a.INTERVAL)); }

IAInterval log(const IAInterval &a) { return IAInterval(fi_lib::j_log(a.INTERVAL)); }

IAInterval ln(const IAInterval &a) { return IAInterval(fi_lib::j_log(a.INTERVAL)); }

IAInterval sqrt(const IAInterval &a) { return IAInterval(fi_lib::j_sqrt(a.INTERVAL)); }

IAInterval sqr(const IAInterval &a) { return IAInterval(fi_lib::j_sqr(a.INTERVAL)); }

IAInterval asnh(const IAInterval &a) { return IAInterval(fi_lib::j_asnh(a.INTERVAL)); }

IAInterval asinh(const IAInterval &a) { return IAInterval(fi_lib::j_asnh(a.INTERVAL)); }

IAInterval acsh(const IAInterval &a) { return IAInterval(fi_lib::j_acsh(a.INTERVAL)); }

IAInterval acosh(const IAInterval &a) { return IAInterval(fi_lib::j_acsh(a.INTERVAL)); }

IAInterval acth(const IAInterval &a) { return IAInterval(fi_lib::j_acth(a.INTERVAL)); }

IAInterval acoth(const IAInterval &a) { return IAInterval(fi_lib::j_acth(a.INTERVAL)); }

IAInterval atnh(const IAInterval &a) { return IAInterval(fi_lib::j_atnh(a.INTERVAL)); }

IAInterval atanh(const IAInterval &a) { return IAInterval(fi_lib::j_atnh(a.INTERVAL)); }

IAInterval asin(const IAInterval &a) { return IAInterval(fi_lib::j_asin(a.INTERVAL)); }

IAInterval acos(const IAInterval &a) { return IAInterval(fi_lib::j_acos(a.INTERVAL)); }

IAInterval acot(const IAInterval &a) { return IAInterval(fi_lib::j_acot(a.INTERVAL)); }

IAInterval atan(const IAInterval &a) { return IAInterval(fi_lib::j_atan(a.INTERVAL)); }

IAInterval sin(const IAInterval &a) { return IAInterval(fi_lib::j_sin(a.INTERVAL)); }

IAInterval cos(const IAInterval &a) { return IAInterval(fi_lib::j_cos(a.INTERVAL)); }

IAInterval cot(const IAInterval &a) { return IAInterval(fi_lib::j_cot(a.INTERVAL)); }

IAInterval tan(const IAInterval &a) { return IAInterval(fi_lib::j_tan(a.INTERVAL)); }

IAInterval exp2(const IAInterval &a) { return IAInterval(fi_lib::j_exp2(a.INTERVAL)); }

IAInterval ex10(const IAInterval &a) { return IAInterval(fi_lib::j_ex10(a.INTERVAL)); }

IAInterval log2(const IAInterval &a) { return IAInterval(fi_lib::j_log2(a.INTERVAL)); }

IAInterval lg10(const IAInterval &a) { return IAInterval(fi_lib::j_lg10(a.INTERVAL)); }

IAInterval erf(const IAInterval &a) { return IAInterval(fi_lib::j_erf(a.INTERVAL)); }

IAInterval erfc(const IAInterval &a) { return IAInterval(fi_lib::j_erfc(a.INTERVAL)); }

const IAInterval &IAInterval::c2r() {
  static const IAInterval c_2r = IAInterval(1.0) / 2.0;
  return c_2r;
}

const IAInterval &IAInterval::c3r() {
  static const IAInterval c_3r = IAInterval(1.0) / 3.0;
  return c_3r;
}

const IAInterval &IAInterval::c5r() {
  static const IAInterval c_5r = IAInterval(1.0) / 5.0;
  return c_5r;
}

const IAInterval &IAInterval::c6r() {
  static const IAInterval c_6r = IAInterval(1.0) / 6.0;
  return c_6r;
}

IAInterval IAInterval::Pi() { return IAInterval(cxsc::Pi_interval); }

IAInterval IAInterval::Pi2() { return IAInterval(cxsc::Pi2_interval); }

IAInterval IAInterval::Pi3() { return IAInterval(cxsc::Pi3_interval); }

IAInterval IAInterval::Pid2() { return IAInterval(cxsc::Pid2_interval); }

IAInterval IAInterval::Pid3() { return IAInterval(cxsc::Pid3_interval); }

IAInterval IAInterval::Pid4() { return IAInterval(cxsc::Pid4_interval); }

IAInterval IAInterval::Pir() { return IAInterval(cxsc::Pir_interval); }

IAInterval IAInterval::Pi2r() { return IAInterval(cxsc::Pi2r_interval); }

IAInterval IAInterval::Pip2() { return IAInterval(cxsc::Pip2_interval); }

IAInterval IAInterval::PiSqr() { return IAInterval::Pip2(); }

IAInterval IAInterval::SqrtPi() { return IAInterval(cxsc::SqrtPi_interval); }

IAInterval IAInterval::Sqrt2Pi() { return IAInterval(cxsc::Sqrt2Pi_interval); }

IAInterval IAInterval::SqrtPir() { return IAInterval(cxsc::SqrtPir_interval); }

IAInterval IAInterval::Sqrt2Pir() { return IAInterval(cxsc::Sqrt2Pir_interval); }

IAInterval IAInterval::Sqrt2() { return IAInterval(cxsc::Sqrt2_interval); }

IAInterval IAInterval::Sqrt3() { return IAInterval(cxsc::Sqrt3_interval); }

IAInterval IAInterval::Sqrt5() { return IAInterval(cxsc::Sqrt5_interval); }

IAInterval IAInterval::Sqrt7() { return IAInterval(cxsc::Sqrt7_interval); }

IAInterval IAInterval::Sqrt2r() { return IAInterval(cxsc::Sqrt2r_interval); }

IAInterval IAInterval::Sqrt3r() { return IAInterval(cxsc::Sqrt3r_interval); }

IAInterval IAInterval::Sqrt3d2() { return IAInterval(cxsc::Sqrt3d2_interval); }

IAInterval IAInterval::Ln2() { return IAInterval(cxsc::Ln2_interval); }

IAInterval IAInterval::Ln10() { return IAInterval(cxsc::Ln10_interval); }

IAInterval IAInterval::LnPi() { return IAInterval(cxsc::LnPi_interval); }

IAInterval IAInterval::Ln2Pi() { return IAInterval(cxsc::Ln2Pi_interval); }

IAInterval IAInterval::Ln2r() { return IAInterval(cxsc::Ln2r_interval); }

IAInterval IAInterval::Ln10r() { return IAInterval(cxsc::Ln10r_interval); }

IAInterval IAInterval::E() { return IAInterval(cxsc::E_interval); }

IAInterval IAInterval::Er() { return IAInterval(cxsc::Er_interval); }

IAInterval IAInterval::Ep2() { return IAInterval(cxsc::Ep2_interval); }

IAInterval IAInterval::ESqr() { return IAInterval::Ep2(); }

IAInterval IAInterval::Ep2r() { return IAInterval(cxsc::Ep2r_interval); }

IAInterval IAInterval::EpPi() { return IAInterval(cxsc::EpPi_interval); }

IAInterval IAInterval::Ep2Pi() { return IAInterval(cxsc::Ep2Pi_interval); }

IAInterval IAInterval::EpPid2() { return IAInterval(cxsc::EpPid2_interval); }

IAInterval IAInterval::EpPid4() { return IAInterval(cxsc::EpPid4_interval); }
