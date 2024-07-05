#ifndef COMPLEXFUNCTIONS_HPP
#define COMPLEXFUNCTIONS_HPP

#include <complex>

inline double baAbsSqr(const std::complex<double> &z) { return std::norm(z); }

inline double baReal(const std::complex<double> &z) { return std::real(z); }

inline double baImag(const std::complex<double> &z) { return std::imag(z); }

inline std::complex<double> baConj(const std::complex<double> &z) { return std::conj(z); }

#ifdef BUILD_IA

inline IAInterval baAbsSqr(const IACInterval &z) { return absSqr(z); }

inline IAInterval baReal(const IACInterval &z) { return real(z); }

inline IAInterval baImag(const IACInterval &z) { return imag(z); }

inline IACInterval baConj(const IACInterval &z) { return conj(z); }

#endif

#endif // COMPLEXFUNCTIONS_HPP
