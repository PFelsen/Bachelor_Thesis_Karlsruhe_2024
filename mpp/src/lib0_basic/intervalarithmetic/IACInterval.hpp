#ifndef IACINTERVAL_HPP
#define IACINTERVAL_HPP

#include "IAInterval.hpp"

#include <complex>
#include <cinterval.hpp>
#include <json.hpp>

class IACInterval {

  cxsc::cinterval CINTERVAL;

  IACInterval(const cxsc::cinterval &z) : CINTERVAL(z) {}

  static const cxsc::interval &extractInterval(const IAInterval &a) { return a.INTERVAL; }
public:
  constexpr IACInterval() : CINTERVAL(0.0) {}

  constexpr IACInterval(const IACInterval &z) = default;

  constexpr IACInterval(IACInterval &&z) : CINTERVAL(std::move(z.CINTERVAL)) {}

  IACInterval(const std::complex<double> &z) :
      CINTERVAL(cxsc::interval(std::real(z)), cxsc::interval(std::imag(z))) {}

  constexpr IACInterval(const IAInterval &re, const IAInterval &im = IAInterval()) :
      CINTERVAL(re.INTERVAL, im.INTERVAL) {}

  IACInterval(double re, double im = 0.0) : CINTERVAL(cxsc::interval(re), cxsc::interval(im)) {}

  constexpr IACInterval &operator=(const IACInterval &z) {
    CINTERVAL = z.CINTERVAL;
    return *this;
  }

  constexpr IACInterval &operator=(const std::complex<double> &z) {
    CINTERVAL = cxsc::cinterval(cxsc::interval(std::real(z)), cxsc::interval(std::imag(z)));
    return *this;
  }

  IACInterval &operator=(const IAInterval &re) {
    CINTERVAL = re.INTERVAL;
    return *this;
  }

  IACInterval &operator=(const double &re) {
    CINTERVAL = re;
    return *this;
  }

  template<typename T>
  IACInterval &operator+=(const T &a) {
    *this = *this + a;
    return *this;
  }

  template<typename T>
  IACInterval &operator-=(const T &a) {
    *this = *this - a;
    return *this;
  }

  template<typename T>
  IACInterval &operator*=(const T &a) {
    *this = *this * a;
    return *this;
  }

  template<typename T>
  IACInterval &operator/=(const T &a) {
    *this = *this / a;
    return *this;
  }

  IAInterval real() const;

  void setReal(const IAInterval &a);

  IAInterval imag() const;

  void setImag(const IAInterval &a);

  IACInterval conj() const;

  IAInterval absSqr() const;

  IAInterval abs() const;

  std::complex<double> mid() const;

  friend IACInterval operator-(const IACInterval &a);

  friend IACInterval operator+(const IACInterval &a, const IACInterval &b);

  friend IACInterval operator+(const IACInterval &a, const std::complex<double> &b);

  friend IACInterval operator+(const std::complex<double> &a, const IACInterval &b);

  friend IACInterval operator+(const IACInterval &a, const IAInterval &b);

  friend IACInterval operator+(const IAInterval &a, const IACInterval &b);

  friend IACInterval operator+(const IACInterval &a, const double &b);

  friend IACInterval operator+(const double &a, const IACInterval &b);

  friend IACInterval operator-(const IACInterval &a, const IACInterval &b);

  friend IACInterval operator-(const IACInterval &a, const std::complex<double> &b);

  friend IACInterval operator-(const std::complex<double> &a, const IACInterval &b);

  friend IACInterval operator-(const IACInterval &a, const IAInterval &b);

  friend IACInterval operator-(const IAInterval &a, const IACInterval &b);

  friend IACInterval operator-(const IACInterval &a, const double &b);

  friend IACInterval operator-(const double &a, const IACInterval &b);

  friend IACInterval operator*(const IACInterval &a, const IACInterval &b);

  friend IACInterval operator*(const IACInterval &a, const std::complex<double> &b);

  friend IACInterval operator*(const std::complex<double> &a, const IACInterval &b);

  friend IACInterval operator*(const IACInterval &a, const IAInterval &b);

  friend IACInterval operator*(const IAInterval &a, const IACInterval &b);

  friend IACInterval operator*(const IACInterval &a, const double &b);

  friend IACInterval operator*(const double &a, const IACInterval &b);

  friend IACInterval operator/(const IACInterval &a, const IACInterval &b);

  friend IACInterval operator/(const IACInterval &a, const std::complex<double> &b);

  friend IACInterval operator/(const std::complex<double> &a, const IACInterval &b);

  friend IACInterval operator/(const IACInterval &a, const IAInterval &b);

  friend IACInterval operator/(const IAInterval &a, const IACInterval &b);

  friend IACInterval operator/(const IACInterval &a, const double &b);

  friend IACInterval operator/(const double &a, const IACInterval &b);

  friend bool operator==(const IACInterval &a, const IACInterval &b);

  friend bool operator!=(const IACInterval &a, const IACInterval &b);

  friend IACInterval accumulate(const std::vector<IACInterval> &a,
                                const std::vector<std::complex<double>> &b);

  friend IACInterval accumulate(const std::vector<std::complex<double>> &a,
                                const std::vector<IACInterval> &b);

  friend std::ostream &operator<<(std::ostream &os, const IACInterval &a);

  friend inline void to_json(nlohmann::json &j, const IACInterval &i);
};

inline IAInterval real(const IACInterval &z) { return z.real(); }

inline IAInterval imag(const IACInterval &z) { return z.imag(); }

inline IACInterval conj(const IACInterval &z) { return z.conj(); }

inline IAInterval absSqr(const IACInterval &z) { return z.absSqr(); }

inline IAInterval abs(const IACInterval &z) { return z.abs(); }

inline std::complex<double> mid(const IACInterval &z) { return z.mid(); }

inline void to_json(nlohmann::json &j, const IACInterval &i) {
  nlohmann::json imagJson(imag(i));
  nlohmann::json realJson(real(i));

  j = nlohmann::json{{"im", imagJson}, {"rel", realJson}};
}

#endif // IACINTERVAL_HPP
