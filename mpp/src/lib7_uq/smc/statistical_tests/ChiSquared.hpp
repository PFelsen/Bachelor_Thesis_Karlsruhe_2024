#ifndef CHI_SQUARE_H
#define CHI_SQUARE_H

#include "StatisticalTest.hpp"

class ChiSquaredTest : public StatisticalTest1D {
private:
  long sample_size;
  const double min_expected_samples_in_bin = 5;
  const double bin_threshold = 130;
  const int min_degrees_of_freedom = 6;
public:
  ChiSquaredTest(TargetDistributionConfig conf) : StatisticalTest1D(conf) {}

  std::vector<int> count_bin_samples(std::vector<double> sample,
                                     std::vector<double> binning_delimiters,
                                     std::vector<double> weights = {});
  double test_statistic(std::vector<double> sample, std::vector<double> weights = {});
  double chi2_statistic(std::vector<double> probabilities, std::vector<int> sample_counts);
  double p_value(double chi_squared, int degree_of_freedom);
  bool Run(std::vector<double> sample, double significance_level, std::vector<double> weights = {});
};


#endif // CHI_SQUARE_H