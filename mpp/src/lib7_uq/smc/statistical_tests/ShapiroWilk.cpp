#include "ShapiroWilk.hpp"

#include <math.h>

/*
  Code taken from http://www.matrixscience.com/msparser/help/group___shapiro_wilk_source_code.html
  Changes regard the interface (and possibly some error handling) only.
  We assume not having any censored samples.
*/

bool ShapiroWilkTest::Run(std::vector<double> sample, double significance_level,
                          std::vector<double> weights) {
  /* Algorithm AS R94, Journal of the Royal Statistical Society Series C
   * (Applied Statistics) vol. 44, no. 4, pp. 547-551 (1995).
   */
  if (sample.size() < 3 || sample.size() > 5000)
    ; // error

  double w = test_statistic(sample);
  double p = p_value(w, sample.size());

  return p > significance_level;
}

std::vector<double> ShapiroWilkTest::get_coefficients_a(int n) {
  std::vector<double> a(n);
  int i1;
  double an25, summ2, ssumm2, rsn, a1, a2, fac;


  static const double c1[6] = {0., .221157, -.147981, -2.07119, 4.434685, -2.706056},
                      c2[6] = {0., .042981, -.293762, -1.752461, 5.682633, -3.582633};

  if (n == 3) {
    a[0] = M_SQRT1_2;
  } else {
    an25 = n + .25;
    summ2 = 0;
    for (int i = 1; i <= n / 2; ++i) {
      a[i - 1] = (double)ppnd7((i - .375f) / an25);
      summ2 += (a[i - 1] * a[i - 1]);
    }
    summ2 *= 2;
    ssumm2 = std::sqrt(summ2);
    rsn = 1 / std::sqrt(n);
    a1 = poly(c1, 6, rsn) - a[0] / ssumm2;

    /* Normalize a[] */
    if (n > 5) {
      i1 = 3;
      a2 = -a[1] / ssumm2 + poly(c2, 6, rsn);
      fac =
          std::sqrt((summ2 - 2 * a[0] * a[0] - 2 * a[1] * a[1]) / (1 - 2 * a1 * a1 - 2 * a2 * a2));
      a[1] = a2;
    } else {
      i1 = 2;
      fac = std::sqrt((summ2 - 2 * a[0] * a[0]) / (1 - 2 * a1 * a1));
    }
    a[0] = a1;
    for (int i = i1; i <= n / 2; ++i)
      a[i - 1] /= -fac;
  }
  return a;
}

double ShapiroWilkTest::p_value(double w, long n) {
  double pw, y, xx, s, m, gamma;

  static const double small_value = 1e-19;
  static const double pi6 = 1.909859;
  static const double stqr = 1.047198;

  static const double g[2] = {-2.273, .459};
  static const double c3[4] = {.544, -.39978, .025054, -6.714e-4};
  static const double c4[4] = {1.3822, -.77857, .062767, -.0020322};
  static const double c5[4] = {-1.5861, -.31082, -.083751, .0038915};
  static const double c6[3] = {-.4803, -.082676, .0030302};

  w = 1. - w;
  if (n == 3) { return pi6 * (std::asin(std::sqrt(w)) - stqr); }
  y = std::log(1 - w);
  xx = std::log(n);
  m = 0;
  s = 1;
  if (n <= 11) {
    gamma = poly(g, 2, n);
    if (y >= gamma) {
      pw = small_value; /* FIXME: rather use an even small_valueer value, or NA ? */
                        // error!!!
    }
    y = -std::log(gamma - y);
    m = poly(c3, 4, n);
    s = std::exp(poly(c4, 4, n));
  } else { /* n >= 12 */
    m = poly(c5, 4, xx);
    s = std::exp(poly(c6, 3, xx));
  }
  /*DBG printf("c(w1=%g, w=%g, y=%g, m=%g, s=%g)\n",w1,*w,y,m,s); */

  return alnorm((y - m) / s, true);
}

double ShapiroWilkTest::test_statistic(std::vector<double> x, std::vector<double> weights) {
  if (!weights.empty()) {
    throw std::invalid_argument("Shapiro Wilk test has no implementation for weighted samples.");
  }
  long n = x.size();
  std::vector<double> a = get_coefficients_a(n);
  int j;
  double sa, sx, xx, xi, w1, range;
  double asa, ssa, sax, ssx, xsx, ssassx;

  range = x[n - 1] - x[0];
  xx = x[0] / range;
  sx = xx;
  sa = -a[0];
  j = n - 1;
  for (int i = 2; i <= n; ++i) {
    xi = x[i - 1] / range;
    sx += xi;
    if (i != j) sa += sign(1, i - j) * a[std::min(i, j) - 1];
    xx = xi;
    --j;
  }

  sa /= n;
  sx /= n;
  ssa = ssx = sax = 0;
  j = n;
  for (int i = 1; i <= n; ++i, --j) {
    if (i != j) asa = sign(1, i - j) * a[std::min(i, j) - 1] - sa;
    else asa = -sa;
    xsx = x[i - 1] / range - sx;
    ssa += asa * asa;
    ssx += xsx * xsx;
    sax += asa * xsx;
  }

  /*  W1 equals (1-W) claculated to avoid excessive rounding error
      for W very near 1 (a potential problem in very large samples) */

  ssassx = std::sqrt(ssa * ssx);
  w1 = (ssassx - sax) * (ssassx + sax) / (ssa * ssx);
  return w1;
}

double ShapiroWilkTest::poly(const double *cc, int nord, double x) {
  /* Auxiliary procedure in algorithm AS 181, Journal of the Royal
   * Statistical Society Series C (Applied Statistics) vol. 31, no. 2,
   * pp. 176-180 (1982).
   */

  /* Local variables */
  int n2, j;
  double ret_val = cc[0];
  if (nord == 1) return ret_val;
  double p = x * cc[nord - 1];
  if (nord != 2) {
    n2 = nord - 2;
    j = n2 + 1;
    for (int i = 1; i <= n2; i++) {
      p = (p + cc[j - 1]) * x;
      j--;
    }
  }
  ret_val = ret_val + p;
  return ret_val;
} /* poly */

double ShapiroWilkTest::ppnd7(double p) {
  /* Algorithm AS 241, Journal of the Royal Statistical Society Series C
   * (Applied Statistics) vol. 26, no. 3, pp. 118-121 (1977).
   */

  static const double zero = 0.0;
  static const double one = 1.0;
  static const double half = 0.5;
  static const double split1 = 0.425;
  static const double split2 = 5.0;
  static const double const1 = 0.180625;
  static const double const2 = 1.6;
  static const double a0 = 3.3871327179E+00;
  static const double a1 = 5.0434271938E+01;
  static const double a2 = 1.5929113202E+02;
  static const double a3 = 5.9109374720E+01;
  static const double b1 = 1.7895169469E+01;
  static const double b2 = 7.8757757664E+01;
  static const double b3 = 6.7187563600E+01;
  static const double c0 = 1.4234372777E+00;
  static const double c1 = 2.7568153900E+00;
  static const double c2 = 1.3067284816E+00;
  static const double c3 = 1.7023821103E-01;
  static const double d1 = 7.3700164250E-01;
  static const double d2 = 1.2021132975E-01;
  static const double e0 = 6.6579051150E+00;
  static const double e1 = 3.0812263860E+00;
  static const double e2 = 4.2868294337E-01;
  static const double e3 = 1.7337203997E-02;
  static const double f1 = 2.4197894225E-01;
  static const double f2 = 1.2258202635E-02;

  double normal_dev;
  double q;
  double r;

  q = p - half;
  if (std::abs(q) <= split1) {
    r = const1 - q * q;
    normal_dev = q * (((a3 * r + a2) * r + a1) * r + a0) / (((b3 * r + b2) * r + b1) * r + one);
    return normal_dev;
  } else {
    if (q < zero) {
      r = p;
    } else {
      r = one - p;
    }
    if (r <= zero) {
      normal_dev = zero;
      return normal_dev;
    }
    r = std::sqrt(-std::log(r));
    if (r <= split2) {
      r = r - const2;
      normal_dev = (((c3 * r + c2) * r + c1) * r + c0) / ((d2 * r + d1) * r + one);
    } else {
      r = r - split2;
      normal_dev = (((e3 * r + e2) * r + e1) * r + e0) / ((f2 * r + f1) * r + one);
    }
    if (q < zero) { normal_dev = -normal_dev; }
    return normal_dev;
  }
}

double ShapiroWilkTest::alnorm(double x, bool upper) {
  /* Algorithm AS 66, Journal of the Royal Statistical Society Series C
   * (Applied Statistics) vol. 22, pp. 424-427 (1973).
   */

  static const double zero = 0;
  static const double one = 1;
  static const double half = 0.5;
  static const double con = 1.28;
  static const double ltone = 7.0;
  static const double utzero = 18.66;
  static const double p = 0.398942280444;
  static const double q = 0.39990348504;
  static const double r = 0.398942280385;
  static const double a1 = 5.75885480458;
  static const double a2 = 2.62433121679;
  static const double a3 = 5.92885724438;
  static const double b1 = -29.8213557807;
  static const double b2 = 48.6959930692;
  static const double c1 = -3.8052E-8;
  static const double c2 = 3.98064794E-4;
  static const double c3 = -0.151679116635;
  static const double c4 = 4.8385912808;
  static const double c5 = 0.742380924027;
  static const double c6 = 3.99019417011;
  static const double d1 = 1.00000615302;
  static const double d2 = 1.98615381364;
  static const double d3 = 5.29330324926;
  static const double d4 = -15.1508972451;
  static const double d5 = 30.789933034;

  double alnorm;
  double z;
  double y;
  bool up = upper;
  z = x;
  if (z < zero) {
    up = !up;
    z = -z;
  }
  if (z <= ltone || (up && z <= utzero)) {
    y = half * z * z;
    if (z > con) {
      alnorm = r * std::exp(-y)
               / (z + c1
                  + d1 / (z + c2 + d2 / (z + c3 + d3 / (z + c4 + d4 / (z + c5 + d5 / (z + c6))))));
    } else {
      alnorm = half - z * (p - q * y / (y + a1 + b1 / (y + a2 + b2 / (y + a3))));
    }
  } else {
    alnorm = zero;
  }

  if (!up) { alnorm = one - alnorm; }
  return alnorm;
}

long ShapiroWilkTest::sign(long x, long y) {
  if (y < 0) return -std::labs(x);
  else return std::labs(x);
}