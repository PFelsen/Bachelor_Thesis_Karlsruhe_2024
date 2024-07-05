#ifndef SHAPIRO_WILK_H
#define SHAPIRO_WILK_H

#include <typeinfo>
#include "StatisticalTest.hpp"

class ShapiroWilkTest : public StatisticalTest1D {
private:
  std::vector<double> get_coefficients_a(int n);
  double poly(const double *cc, int nord, double x);
  double ppnd7(double p);
  double alnorm(double x, bool upper);
  long sign(long x, long y);
public:
  ShapiroWilkTest(TargetDistributionConfig conf) : StatisticalTest1D(conf) {
    if (conf.name != "Normal")
      ; // error
  }

  double p_value(double w, long n);
  double test_statistic(std::vector<double> x, std::vector<double> weights = {});
  bool Run(std::vector<double> sample, double significance_level, std::vector<double> weights = {});
};

#endif // SHAPIRO_WILK_H