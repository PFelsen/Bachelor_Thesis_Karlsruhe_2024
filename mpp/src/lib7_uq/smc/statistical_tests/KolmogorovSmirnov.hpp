#ifndef KOLMOGOROV_SMIRNOV_H
#define KOLMOGOROV_SMIRNOV_H

#include <array>
#include <map>
#include "StatisticalTest.hpp"

class KolmogorovSmirnovTest : public StatisticalTest1D {
private:
  double weighted_test_statistic(std::vector<double> sample, std::vector<double> weights);
public:
  KolmogorovSmirnovTest(TargetDistributionConfig conf) : StatisticalTest1D(conf) {}

  double test_statistic(std::vector<double> sample, std::vector<double> weights = {});
  double p_value(double D, int sample_size);
  bool Run(std::vector<double> sample, const double significance_level,
           std::vector<double> weights = {});
};

void SortSampleAndWeights(std::vector<double> &sample, std::vector<double> &weights);

#endif // KOLMOGOROV_SMIRNOV_H
