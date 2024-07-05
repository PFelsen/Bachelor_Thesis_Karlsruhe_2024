#ifndef IAINTERVAL_HPP
#define IAINTERVAL_HPP

#include "Assertion.hpp"

#include <cmath>
#include <iosfwd>
#include <memory>
#include <vector>
#include <interval.hpp>
#include <json.hpp>

// To avoid implicit type cast problems
inline double pow(const double &a, int n) { return std::pow(a, double(n)); }

class IAInterval {
  friend class IACInterval;

  cxsc::interval INTERVAL;

  // only for type casts
  explicit IAInterval(const cxsc::interval &a) : INTERVAL(a) {}
public:
  constexpr IAInterval() : INTERVAL(0.0) {}

  constexpr IAInterval(const IAInterval &a) = default;

  constexpr IAInterval(IAInterval &&a) : INTERVAL(std::move(a.INTERVAL)) {}

  /// Note that this only yields a correct IAInterval type if point is an exactly representable
  /// float. To ensure a correct enclosure use operator>>(std::string &input, IAInterval &a) instead
  constexpr IAInterval(double point) : INTERVAL(point) {}

  /// Note that this only yields a correct IAInterval type if inf and sup are exactly representable
  /// floats. To ensure a correct enclosure use operator>>(std::string &input, IAInterval &a)
  /// instead
  constexpr IAInterval(double inf, const double sup) : INTERVAL(inf, sup) {}

  constexpr IAInterval &operator=(const IAInterval &a) {
    INTERVAL = a.INTERVAL;
    return *this;
  }

  constexpr IAInterval &operator=(IAInterval &&a) {
    INTERVAL = std::move(a.INTERVAL);
    return *this;
  }

  /// Note that this only yields a correct IAInterval type if a are exactly representable floats.
  constexpr IAInterval &operator=(double a) {
    INTERVAL = a;
    return *this;
  }

  /// converts strings like [0.1, 0.3] to IAInterval type with correct roundings
  friend void operator>>(std::string &input, IAInterval &a);

  /// converts chars like [0.1, 0.3] to IAInterval type with correct roundings
  friend void operator>>(const char *input, IAInterval &a);

  IAInterval &operator+=(const IAInterval &a);

  IAInterval &operator+=(const double &a);

  IAInterval &operator-=(const IAInterval &a);

  IAInterval &operator-=(const double &a);

  IAInterval &operator*=(const IAInterval &a);

  IAInterval &operator*=(const double &a);

  IAInterval &operator/=(const IAInterval &a);

  IAInterval &operator/=(const double &a);

  double inf() const;

  double sup() const;

  void setInf(double inf);

  void setSup(double sup);

  IAInterval abs() const;

  double mid() const;

  double diam() const;

  friend IAInterval operator-(const IAInterval &a);

  friend IAInterval operator+(const IAInterval &a, const IAInterval &b);

  friend IAInterval operator+(const double &a, const IAInterval &b);

  friend IAInterval operator+(const IAInterval &a, const double &b);

  friend IAInterval operator-(const IAInterval &a, const IAInterval &b);

  friend IAInterval operator-(const double &a, const IAInterval &b);

  friend IAInterval operator-(const IAInterval &a, const double &b);

  friend IAInterval operator*(const IAInterval &a, const IAInterval &b);

  friend IAInterval operator*(const double &a, const IAInterval &b);

  friend IAInterval operator*(const IAInterval &a, const double &b);

  friend IAInterval operator/(const IAInterval &a, const IAInterval &b);

  friend IAInterval operator/(const double &a, const IAInterval &b);

  friend IAInterval operator/(const IAInterval &a, const double &b);

  friend bool operator==(const IAInterval &a, const IAInterval &b);

  friend bool operator!=(const IAInterval &a, const IAInterval &b);

  friend bool operator==(const IAInterval &a, const double &b);

  friend bool operator!=(const IAInterval &a, const double &b);

  friend bool operator==(const double &a, const IAInterval &b);

  friend bool operator!=(const double &a, const IAInterval &b);

  friend IAInterval operator|(const IAInterval &a, const IAInterval &b);

  friend IAInterval hull(const IAInterval &a, const IAInterval &b);

  friend IAInterval operator&(const IAInterval &a, const IAInterval &b);

  friend IAInterval intsec(const IAInterval &a, const IAInterval &b);

  friend bool disjoint(const IAInterval &a, const IAInterval &b);

  friend bool operator<=(const double &a, const IAInterval &b);

  friend void to_json(nlohmann::json &j, const IAInterval &i);

  friend bool in(const double &a, const IAInterval &b);

  friend bool operator<=(const IAInterval &a, const IAInterval &b);

  friend bool in(const IAInterval &a, const IAInterval &b);

  friend bool operator<(const IAInterval &a, const IAInterval &b);

  friend bool operator<(const double &a, const IAInterval &b);

  friend bool operator>=(const IAInterval &a, const double &b);

  friend bool operator>=(const IAInterval &a, const IAInterval &b);

  friend bool operator>(const IAInterval &a, const double &b);

  friend bool operator>(const IAInterval &a, const IAInterval &b);

  friend IAInterval pow(const IAInterval &a, const IAInterval &b);

  friend IAInterval pow(const IAInterval &a, int b);

  friend std::ostream &operator<<(std::ostream &os, const IAInterval &a);

  friend IAInterval accumulate(const std::vector<IAInterval> &a, const std::vector<IAInterval> &b);

  friend IAInterval accumulate(const std::vector<IAInterval> &a, const std::vector<double> &b);

  friend IAInterval accumulate(const std::vector<double> &a, const std::vector<IAInterval> &b);

  friend IAInterval accumulate(const std::vector<double> &a, const std::vector<double> &b);

  friend IAInterval exp(const IAInterval &a);

  friend IAInterval sinh(const IAInterval &a);

  friend IAInterval cosh(const IAInterval &a);

  friend IAInterval coth(const IAInterval &a);

  friend IAInterval tanh(const IAInterval &a);

  friend IAInterval log(const IAInterval &a);

  friend IAInterval ln(const IAInterval &a);

  friend IAInterval sqrt(const IAInterval &a);

  friend IAInterval sqr(const IAInterval &a);

  friend IAInterval asnh(const IAInterval &a);

  friend IAInterval asinh(const IAInterval &a);

  friend IAInterval acsh(const IAInterval &a);

  friend IAInterval acosh(const IAInterval &a);

  friend IAInterval acth(const IAInterval &a);

  friend IAInterval acoth(const IAInterval &a);

  friend IAInterval atnh(const IAInterval &a);

  friend IAInterval atanh(const IAInterval &a);

  friend IAInterval asin(const IAInterval &a);

  friend IAInterval acos(const IAInterval &a);

  friend IAInterval acot(const IAInterval &a);

  friend IAInterval atan(const IAInterval &a);

  friend IAInterval sin(const IAInterval &a);

  friend IAInterval cos(const IAInterval &a);

  friend IAInterval cot(const IAInterval &a);

  friend IAInterval tan(const IAInterval &a);

  friend IAInterval exp2(const IAInterval &a);

  friend IAInterval ex10(const IAInterval &a);

  friend IAInterval log2(const IAInterval &a);

  friend IAInterval lg10(const IAInterval &a);

  friend IAInterval erf(const IAInterval &a);

  friend IAInterval erfc(const IAInterval &a);

  /// constant 1/2
  static const IAInterval &c2r();

  /// constant 1/3
  static const IAInterval &c3r();

  /// constant 1/5
  static const IAInterval &c5r();

  /// constant 1/6
  static const IAInterval &c6r();

  /// constant Pi
  static IAInterval Pi();

  /// constant 2*Pi
  static IAInterval Pi2();

  /// constant 3*Pi
  static IAInterval Pi3();

  /// constant Pi/2
  static IAInterval Pid2();

  /// constant Pi/3
  static IAInterval Pid3();

  /// constant Pi/4
  static IAInterval Pid4();

  /// constant 1/Pi
  static IAInterval Pir();

  /// constant 1/(2*Pi)
  static IAInterval Pi2r();

  /// constant Pi^2
  static IAInterval Pip2();
  static IAInterval PiSqr();

  /// constant sqrt(Pi)
  static IAInterval SqrtPi();

  /// constant sqrt(2*Pi)
  static IAInterval Sqrt2Pi();

  /// constant 1/sqrt(Pi)
  static IAInterval SqrtPir();

  /// constant 1/sqrt(2*Pi)
  static IAInterval Sqrt2Pir();

  /// constant sqrt(2)
  static IAInterval Sqrt2();

  /// constant sqrt(3)
  static IAInterval Sqrt3();

  /// constant sqrt(5)
  static IAInterval Sqrt5();

  /// constant sqrt(7)
  static IAInterval Sqrt7();

  /// constant 1/sqrt(2)
  static IAInterval Sqrt2r();

  /// constant 1/sqrt(3)
  static IAInterval Sqrt3r();

  /// constant sqrt(3)/2
  static IAInterval Sqrt3d2();

  /// constant ln(2)
  static IAInterval Ln2();

  /// constant ln(10)
  static IAInterval Ln10();

  /// constant ln(Pi)
  static IAInterval LnPi();

  /// constant ln(2*Pi)
  static IAInterval Ln2Pi();

  /// constant 1/ln(2)
  static IAInterval Ln2r();

  /// constant 1/ln(10)
  static IAInterval Ln10r();

  /// constant e
  static IAInterval E();

  /// constant 1/e
  static IAInterval Er();

  /// constant e^2
  static IAInterval Ep2();
  static IAInterval ESqr();

  /// constant 1/e^2
  static IAInterval Ep2r();

  /// constant e^Pi
  static IAInterval EpPi();

  /// constant e^(2*Pi)
  static IAInterval Ep2Pi();

  /// constant e^(Pi/2)
  static IAInterval EpPid2();

  /// constant e^(Pi/4)
  static IAInterval EpPid4();
};

inline double inf(const IAInterval &a) { return a.inf(); }

inline double sup(const IAInterval &a) { return a.sup(); }

inline IAInterval abs(const IAInterval &a) { return a.abs(); }

inline double mid(const IAInterval &a) { return a.mid(); }

inline double diam(const IAInterval &a) { return a.diam(); }

template<>
struct std::is_scalar<IAInterval> : std::true_type {};

inline void to_json(nlohmann::json &j, const IAInterval &i) {
  j["inf"] = inf(i);
  j["sup"] = sup(i);
}
#endif // IAINTERVAL_HPP
