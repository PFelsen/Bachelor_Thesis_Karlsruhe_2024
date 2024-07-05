#ifndef __NEW_POLYNOMIALS
#define __NEW_POLYNOMIALS

#include <cinttypes>
#include <Assertion.hpp>
#include "PolyCalculationType.hpp"

namespace mpp::shape {

#define FUNCTION_FOR_DEGREE(targetDegree)                                                          \
  static_assert(targetDegree != degree || targetDegree >= caseid,                                  \
                "Degree must be equal or less then the function!");                                \
  if constexpr (targetDegree == degree)

#define CASE_FOR(index, add)                                                                       \
  if constexpr (index == caseid) return Type(add)

template<uint32_t degree, uint32_t caseid, PolyCalculationType Type = double>
constexpr Type polynomial(const Type x) {
  FUNCTION_FOR_DEGREE(0) { CASE_FOR(0, 1.0); }
  FUNCTION_FOR_DEGREE(1) {
    CASE_FOR(0, 1.0 - x);
    CASE_FOR(1, x);
  }
  FUNCTION_FOR_DEGREE(2) {
    CASE_FOR(0, (2 * x - 1) * (x - 1));
    CASE_FOR(1, 4 * x * (1 - x));
    CASE_FOR(2, x * (2 * x - 1));
  }
  FUNCTION_FOR_DEGREE(3) {
    CASE_FOR(0, -1 / Type(2.0) * (3 * x - 1) * (3 * x - 2) * (x - 1));
    CASE_FOR(1, 9 / Type(2.0) * x * (3 * x - 2) * (x - 1));
    CASE_FOR(2, -9 / Type(2.0) * x * (3 * x - 1) * (x - 1));
    CASE_FOR(3, 1 / Type(2.0) * x * (3 * x - 1) * (3 * x - 2));
  }
  FUNCTION_FOR_DEGREE(4) {
    CASE_FOR(0, 1 / Type(6.0) * (4 * x - 1) * (4 * x - 2) * (4 * x - 3) * (x - 1));
    CASE_FOR(1, -8 / Type(3.0) * x * (4 * x - 2) * (4 * x - 3) * (x - 1));
    CASE_FOR(2, 4 * x * (4 * x - 1) * (4 * x - 3) * (x - 1));
    CASE_FOR(3, -8 / Type(3.0) * x * (4 * x - 1) * (4 * x - 2) * (x - 1));
    CASE_FOR(4, 1 / Type(6.0) * x * (4 * x - 1) * (4 * x - 2) * (4 * x - 3));
  }
  FUNCTION_FOR_DEGREE(5) {
    CASE_FOR(0,
             -1.0 / Type(24.0) * (x - 1) * (5 * x - 4) * (5 * x - 3) * (5 * x - 2) * (5 * x - 1));
    CASE_FOR(1, 25.0 / Type(24.0) * x * (x - 1) * (5 * x - 4) * (5 * x - 3) * (5 * x - 2));
    CASE_FOR(2, -25.0 / Type(12.0) * x * (x - 1) * (5 * x - 4) * (5 * x - 3) * (5 * x - 1));
    CASE_FOR(3, 25.0 / Type(12.0) * x * (x - 1) * (5 * x - 4) * (5 * x - 2) * (5 * x - 1));
    CASE_FOR(4, -25.0 / Type(24.0) * x * (x - 1) * (5 * x - 3) * (5 * x - 2) * (5 * x - 1));
    CASE_FOR(5, 1.0 / Type(24.0) * x * (5 * x - 4) * (5 * x - 3) * (5 * x - 2) * (5 * x - 1));
  }
  FUNCTION_FOR_DEGREE(6) {
    CASE_FOR(0, 1.0 / Type(10.0) * (x - 1) * (2 * x - 1) * (3 * x - 2) * (3 * x - 1) * (6 * x - 5)
                    * (6 * x - 1));
    CASE_FOR(1, -18.0 / Type(5.0) * x * (x - 1) * (2 * x - 1) * (3 * x - 2) * (3 * x - 1)
                    * (6 * x - 5));
    CASE_FOR(2,
             9.0 / Type(2.0) * x * (x - 1) * (2 * x - 1) * (3 * x - 2) * (6 * x - 5) * (6 * x - 1));
    CASE_FOR(3, -4 * x * (x - 1) * (3 * x - 2) * (3 * x - 1) * (6 * x - 5) * (6 * x - 1));
    CASE_FOR(4,
             9.0 / Type(2.0) * x * (x - 1) * (2 * x - 1) * (3 * x - 1) * (6 * x - 5) * (6 * x - 1));
    CASE_FOR(5, -18.0 / Type(5.0) * x * (x - 1) * (2 * x - 1) * (3 * x - 2) * (3 * x - 1)
                    * (6 * x - 1));
    CASE_FOR(6, 1.0 / Type(10.0) * x * (2 * x - 1) * (3 * x - 2) * (3 * x - 1) * (6 * x - 5)
                    * (6 * x - 1));
  }
  FUNCTION_FOR_DEGREE(7) {
    CASE_FOR(0, -1.0 / Type(720.0) * (x - 1) * (7 * x - 6) * (7 * x - 5) * (7 * x - 4) * (7 * x - 3)
                    * (7 * x - 2) * (7 * x - 1));
    CASE_FOR(1, 49.0 / Type(720.0) * x * (x - 1) * (7 * x - 6) * (7 * x - 5) * (7 * x - 4)
                    * (7 * x - 3) * (7 * x - 2));
    CASE_FOR(2, -49.0 / Type(240.0) * x * (x - 1) * (7 * x - 6) * (7 * x - 5) * (7 * x - 4)
                    * (7 * x - 3) * (7 * x - 1));
    CASE_FOR(3, 49.0 / Type(144.0) * x * (x - 1) * (7 * x - 6) * (7 * x - 5) * (7 * x - 4)
                    * (7 * x - 2) * (7 * x - 1));
    CASE_FOR(4, -49.0 / Type(144.0) * x * (x - 1) * (7 * x - 6) * (7 * x - 5) * (7 * x - 3)
                    * (7 * x - 2) * (7 * x - 1));
    CASE_FOR(5, 49.0 / Type(240.0) * x * (x - 1) * (7 * x - 6) * (7 * x - 4) * (7 * x - 3)
                    * (7 * x - 2) * (7 * x - 1));
    CASE_FOR(6, -49.0 / Type(720.0) * x * (x - 1) * (7 * x - 5) * (7 * x - 4) * (7 * x - 3)
                    * (7 * x - 2) * (7 * x - 1));
    CASE_FOR(7, 1.0 / Type(720.0) * x * (7 * x - 6) * (7 * x - 5) * (7 * x - 4) * (7 * x - 3)
                    * (7 * x - 2) * (7 * x - 1));
  }
  FUNCTION_FOR_DEGREE(8) {
    CASE_FOR(0, 1.0 / Type(315.0) * (x - 1) * (2 * x - 1) * (4 * x - 3) * (4 * x - 1) * (8 * x - 7)
                    * (8 * x - 5) * (8 * x - 3) * (8 * x - 1));
    CASE_FOR(1, -64.0 / Type(315.0) * x * (x - 1) * (2 * x - 1) * (4 * x - 3) * (4 * x - 1)
                    * (8 * x - 7) * (8 * x - 5) * (8 * x - 3));
    CASE_FOR(2, 16.0 / Type(45.0) * x * (x - 1) * (2 * x - 1) * (4 * x - 3) * (8 * x - 7)
                    * (8 * x - 5) * (8 * x - 3) * (8 * x - 1));
    CASE_FOR(3, -64.0 / Type(45.0) * x * (x - 1) * (2 * x - 1) * (4 * x - 3) * (4 * x - 1)
                    * (8 * x - 7) * (8 * x - 5) * (8 * x - 1));
    CASE_FOR(4, 4.0 / Type(9.0) * x * (x - 1) * (4 * x - 3) * (4 * x - 1) * (8 * x - 7)
                    * (8 * x - 5) * (8 * x - 3) * (8 * x - 1));
    CASE_FOR(5, -64.0 / Type(45.0) * x * (x - 1) * (2 * x - 1) * (4 * x - 3) * (4 * x - 1)
                    * (8 * x - 7) * (8 * x - 3) * (8 * x - 1));
    CASE_FOR(6, 16.0 / Type(45.0) * x * (x - 1) * (2 * x - 1) * (4 * x - 1) * (8 * x - 7)
                    * (8 * x - 5) * (8 * x - 3) * (8 * x - 1));
    CASE_FOR(7, -64.0 / Type(315.0) * x * (x - 1) * (2 * x - 1) * (4 * x - 3) * (4 * x - 1)
                    * (8 * x - 5) * (8 * x - 3) * (8 * x - 1));
    CASE_FOR(8, 1.0 / Type(315.0) * x * (2 * x - 1) * (4 * x - 3) * (4 * x - 1) * (8 * x - 7)
                    * (8 * x - 5) * (8 * x - 3) * (8 * x - 1));
  }
  THROW("Not implemented!")
}

template<uint32_t degree, uint32_t caseid, PolyCalculationType Type = double>
constexpr Type polynomialDerivative(const Type x) noexcept {
  FUNCTION_FOR_DEGREE(0) { CASE_FOR(0, 0.0); }
  FUNCTION_FOR_DEGREE(1) {
    CASE_FOR(0, -1.0);
    CASE_FOR(1, 1.0);
  }
  FUNCTION_FOR_DEGREE(2) {
    CASE_FOR(0, 4 * x - 3);
    CASE_FOR(1, -8 * x + 4);
    CASE_FOR(2, 4 * x - 1);
  }
  FUNCTION_FOR_DEGREE(3) {
    CASE_FOR(0, -27 / Type(2.0) * x * x + 18 * x - 11 / Type(2.0));
    CASE_FOR(1, 81 / Type(2.0) * x * x - 45 * x + 9);
    CASE_FOR(2, -81 / Type(2.0) * x * x + 36 * x - 9 / Type(2.0));
    CASE_FOR(3, 27 / Type(2.0) * x * x - 9 * x + 1);
  }
  FUNCTION_FOR_DEGREE(4) {
    Type x2 = x * x;
    Type x3 = x2 * x;
    CASE_FOR(0, 1 / Type(3.0) * (128 * x3 - 240 * x2 + 140 * x - 25));
    CASE_FOR(1, -1 / Type(3.0) * (512 * x3 - 864 * x2 + 416 * x - 48));
    CASE_FOR(2, 256 * x3 - 384 * x2 + 152 * x - 12);
    CASE_FOR(3, -1 / Type(3.0) * (512 * x3 - 672 * x2 + 224 * x - 16));
    CASE_FOR(4, 1 / Type(3.0) * (128 * x3 - 144 * x2 + 44 * x - 3));
  }
  FUNCTION_FOR_DEGREE(5) {
    Type x2 = x * x;
    Type x3 = x2 * x;
    Type x4 = x2 * x2;
    CASE_FOR(0, -1.0 / Type(24.0) * (3125 * x4 - 7500 * x3 + 6375 * x2 - 2250 * x + 274));
    CASE_FOR(1, 25.0 / Type(24.0) * (625 * x4 - 1400 * x3 + 1065 * x2 - 308 * x + 24));
    CASE_FOR(2, -25.0 / Type(12.0) * (625 * x4 - 1300 * x3 + 885 * x2 - 214 * x + 12));
    CASE_FOR(3, 25.0 / Type(12.0) * (625 * x4 - 1200 * x3 + 735 * x2 - 156 * x + 8));
    CASE_FOR(4, -25.0 / Type(24.0) * (625 * x4 - 1100 * x3 + 615 * x2 - 122 * x + 6));
    CASE_FOR(5, 1.0 / Type(24.0) * (3125 * x4 - 5000 * x3 + 2625 * x2 - 500 * x + 24));
  }
  FUNCTION_FOR_DEGREE(6) {
    Type x2 = x * x;
    Type x3 = x2 * x;
    Type x4 = x2 * x2;
    Type x5 = x2 * x3;
    CASE_FOR(0, 1.0 / Type(10.0) * (12 * x - 7) * (324 * x4 - 756 * x3 + 609 * x2 - 196 * x + 21));
    CASE_FOR(1, -36.0 / Type(5.0) * (324 * x5 - 900 * x4 + 930 * x3 - 435 * x2 + 87 * x - 5));
    CASE_FOR(2, 9.0 / Type(2.0) * (1296 * x5 - 3420 * x4 + 3288 * x3 - 1383 * x2 + 234 * x - 10));
    CASE_FOR(3, -8 * (2 * x - 1) * (18 * x2 - 18 * x + 1) * (27 * x2 - 27 * x + 5));
    CASE_FOR(4, 9.0 / Type(2.0) * (1296 * x5 - 3060 * x4 + 2568 * x3 - 921 * x2 + 132 * x - 5));
    CASE_FOR(5, -36.0 / Type(5.0) * (324 * x5 - 720 * x4 + 570 * x3 - 195 * x2 + 27 * x - 1));
    CASE_FOR(6, 1.0 / Type(10.0) * (12 * x - 5) * (324 * x4 - 540 * x3 + 285 * x2 - 50 * x + 2));
  }
  FUNCTION_FOR_DEGREE(7) {

    Type x2 = x * x;
    Type x3 = x2 * x;
    Type x4 = x2 * x2;
    Type x5 = x2 * x3;
    Type x6 = x3 * x3;
    CASE_FOR(0, -1.0 / Type(720.0)
                    * (823543 * x6 - 2823576 * x5 + 3865610 * x4 - 2689120 * x3 + 995043 * x2
                       - 183848 * x + 13068));
    CASE_FOR(1, 49.0 / Type(720.0)
                    * (117649 * x6 - 388962 * x5 + 505925 * x4 - 326340 * x3 + 107184 * x2
                       - 16056 * x + 720));
    CASE_FOR(2, -49.0 / Type(240.0)
                    * (117649 * x6 - 374556 * x5 + 463050 * x4 - 278320 * x3 + 82509 * x2
                       - 10548 * x + 360));
    CASE_FOR(3, 49.0 / Type(144.0)
                    * (117649 * x6 - 360150 * x5 + 423605 * x4 - 238924 * x3 + 65352 * x2 - 7592 * x
                       + 240));
    CASE_FOR(4, -49.0 / Type(144.0)
                    * (117649 * x6 - 345744 * x5 + 387590 * x4 - 206976 * x3 + 53445 * x2 - 5904 * x
                       + 180));
    CASE_FOR(5, 49.0 / Type(240.0)
                    * (117649 * x6 - 331338 * x5 + 355005 * x4 - 181300 * x3 + 45024 * x2 - 4824 * x
                       + 144));
    CASE_FOR(6, -49.0 / Type(720.0)
                    * (117649 * x6 - 316932 * x5 + 325850 * x4 - 160720 * x3 + 38829 * x2 - 4076 * x
                       + 120));
    CASE_FOR(7, 1.0 / Type(720.0)
                    * (823543 * x6 - 2117682 * x5 + 2100875 * x4 - 1008420 * x3 + 238728 * x2
                       - 24696 * x + 720));
  }
  FUNCTION_FOR_DEGREE(8) {

    Type x2 = x * x;
    Type x3 = x2 * x;
    Type x4 = x2 * x2;
    Type x5 = x2 * x3;
    Type x6 = x3 * x3;
    Type x7 = x4 * x3;
    CASE_FOR(0, 1.0 / Type(315.0) * (16 * x - 9)
                    * (65536 * x6 - 221184 * x5 + 294912 * x4 - 196992 * x3 + 68784 * x2 - 11772 * x
                       + 761));
    CASE_FOR(1, -64.0 / Type(315.0)
                    * (131072 * x7 - 501760 * x6 + 784896 * x5 - 644000 * x4 + 294784 * x3
                       - 73290 * x2 + 8658 * x - 315));
    CASE_FOR(2, 16.0 / Type(45.0)
                    * (262144 * x7 - 974848 * x6 + 1468416 * x5 - 1145600 * x4 + 489248 * x3
                       - 110118 * x2 + 11178 * x - 315));
    CASE_FOR(3, -64.0 / Type(45.0)
                    * (131072 * x7 - 473088 * x6 + 686592 * x5 - 511200 * x4 + 205824 * x3
                       - 43038 * x2 + 4006 * x - 105));
    CASE_FOR(4, 4.0 / Type(9.0) * (2 * x - 1)
                    * (262144 * x6 - 786432 * x5 + 890880 * x4 - 471040 * x3 + 116256 * x2
                       - 11808 * x + 315));
    CASE_FOR(5, -64.0 / Type(45.0)
                    * (131072 * x7 - 444416 * x6 + 600576 * x5 - 412960 * x4 + 152704 * x3
                       - 29346 * x2 + 2538 * x - 63));
    CASE_FOR(6, 16.0 / Type(45.0)
                    * (262144 * x7 - 860160 * x6 + 1124352 * x5 - 748800 * x4 + 269088 * x3
                       - 50490 * x2 + 4286 * x - 105));
    CASE_FOR(7, -64.0 / Type(315.0)
                    * (131072 * x7 - 415744 * x6 + 526848 * x5 - 341600 * x4 + 120064 * x3
                       - 22134 * x2 + 1854 * x - 45));
    CASE_FOR(8, 1.0 / Type(315.0) * (16 * x - 7)
                    * (65536 * x6 - 172032 * x5 + 172032 * x4 - 81536 * x3 + 18480 * x2 - 1764 * x
                       + 45));
  }
}

template<uint32_t degree, uint32_t caseid, PolyCalculationType Type = double>
constexpr Type shapeFunction(const Type x) {
  if constexpr (degree < 7) {
    CASE_FOR(0, 1.0);
    CASE_FOR(1, degree * x);
    CASE_FOR(2, (degree * x) * (degree * x - 1) / Type(2.0));
    CASE_FOR(3, (degree * x) * (degree * x - 1) * (degree * x - 2) / Type(6.0));
    CASE_FOR(4, (degree * x) * (degree * x - 1) * (degree * x - 2) * (degree * x - 3) / Type(24.0));
    CASE_FOR(5, (degree * x) * (degree * x - 1) * (degree * x - 2) * (degree * x - 3)
                    * (degree * x - 4) / Type(120.0));
    CASE_FOR(6, (degree * x) * (degree * x - 1) * (degree * x - 2) * (degree * x - 3)
                    * (degree * x - 4) * (degree * x - 5) / Type(720.0));
  } else {
    THROW("Not implemented!");
  }
}

template<uint32_t degree, uint32_t caseid, PolyCalculationType Type = double>
constexpr Type shapeDerivative(const Type x) {
  if constexpr (degree < 7) {
    double di[degree]{};
    di[0] = degree;
    Type xi[degree]{};
    xi[0] = x;

    for (int power = 1; power < caseid; ++power) {
      xi[power] = xi[power - 1] * x;
      di[power] = di[power - 1] * degree;
    }
    CASE_FOR(0, 0.0);
    CASE_FOR(1, degree);
    CASE_FOR(2, di[1] / Type(1.0) * xi[0] - di[0] / Type(2.0));
    CASE_FOR(3, di[2] / Type(2.0) * xi[1] + (-6.0 * di[1] * xi[0] + 2.0 * di[0]) / Type(6.0));
    CASE_FOR(4, di[3] / Type(6.0) * xi[2]
                    + (-18.0 * di[2] * xi[1] + 22.0 * di[1] * xi[0] - 6.0 * di[0]) / Type(24.0));
    CASE_FOR(5, di[4] / Type(24.0) * xi[3]
                    + (-40.0 * di[3] * xi[2] + 105.0 * di[2] * xi[1] - 100.0 * di[1] * xi[0]
                       + 24.0 * di[0])
                          / Type(120.0));
    CASE_FOR(6, di[5] / Type(120.0) * xi[4]
                    + (-75.0 * di[4] * xi[3] + 340.0 * di[3] * xi[2] - 675.0 * di[2] * xi[1]
                       + 548.0 * di[1] * xi[0] - 120.0 * di[0])
                          / Type(720.0));
  } else {
    THROW("Not implemented!");
  }
}

#define MPP_IMPL_CASE_FUNC(name, index)                                                            \
  case index:                                                                                      \
    return name<index, Type>(degree, input)

#define MPP_IMPL_CASE_DEG(name, index)                                                             \
  case index: {                                                                                    \
    if constexpr (index <= function) {                                                             \
                                                                                                   \
      return name<function, index, Type>(input);                                                   \
    } else {                                                                                       \
      THROW("Not supported!");                                                                     \
    }                                                                                              \
  }

#define MPP_IMPL_LAYERS(name)                                                                      \
  template<uint32_t function, PolyCalculationType Type = double>                                   \
  constexpr Type name(const Type input, const uint32_t degree) {                                   \
    switch (degree) {                                                                              \
      MPP_IMPL_CASE_DEG(name, 0);                                                                  \
      MPP_IMPL_CASE_DEG(name, 1);                                                                  \
      MPP_IMPL_CASE_DEG(name, 2);                                                                  \
      MPP_IMPL_CASE_DEG(name, 3);                                                                  \
      MPP_IMPL_CASE_DEG(name, 4);                                                                  \
      MPP_IMPL_CASE_DEG(name, 5);                                                                  \
      MPP_IMPL_CASE_DEG(name, 6);                                                                  \
      MPP_IMPL_CASE_DEG(name, 7);                                                                  \
      MPP_IMPL_CASE_DEG(name, 8);                                                                  \
    }                                                                                              \
    THROW("Not implemented!");                                                                     \
  }                                                                                                \
                                                                                                   \
  template<PolyCalculationType Type = double>                                                      \
  constexpr Type name(const Type input, const uint32_t degree, const uint32_t function) {          \
    switch (function) {                                                                            \
      MPP_IMPL_CASE_FUNC(name, 0);                                                                 \
      MPP_IMPL_CASE_FUNC(name, 1);                                                                 \
      MPP_IMPL_CASE_FUNC(name, 2);                                                                 \
      MPP_IMPL_CASE_FUNC(name, 3);                                                                 \
      MPP_IMPL_CASE_FUNC(name, 4);                                                                 \
      MPP_IMPL_CASE_FUNC(name, 5);                                                                 \
      MPP_IMPL_CASE_FUNC(name, 6);                                                                 \
      MPP_IMPL_CASE_FUNC(name, 7);                                                                 \
      MPP_IMPL_CASE_FUNC(name, 8);                                                                 \
    }                                                                                              \
    THROW("Not implemented!");                                                                     \
  }

MPP_IMPL_LAYERS(polynomial);
MPP_IMPL_LAYERS(polynomialDerivative);
MPP_IMPL_LAYERS(shapeFunction);
MPP_IMPL_LAYERS(shapeDerivative);

#undef MPP_IMPL_LAYERS
#undef FUNCTION_FOR_DEGREE
#undef MPP_IMPL_CASE_DEG
#undef CASE_FOR
#undef MPP_IMPL_CASE_FUNC

} // namespace mpp::shape

#endif