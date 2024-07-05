#include "ChiSquared.hpp"

#include <algorithm>
#include <numeric>
#include <stdexcept>
#include <utility>
#include <math.h>

/*
  TODO
    Warning, if not enough samples are drawn for binning.
    Does this even work this way for discrete distributions?
*/

bool ChiSquaredTest::Run(std::vector<double> sample, double significance_level,
                         std::vector<double> weights) {
  double chi_squared = test_statistic(sample, weights);
  double degree_of_freedom = binning->GetDegreeOfFreedom();
  double p = p_value(chi_squared, binning->GetDegreeOfFreedom());
  return p > significance_level;
}

double ChiSquaredTest::test_statistic(std::vector<double> sample, std::vector<double> weights) {
  std::vector<double> delimiters = binning->binning_delimiters(sample.size());
  std::vector<double> probabilities = binning->probabilities(delimiters);
  std::vector<int> sample_counts = count_bin_samples(sample, delimiters);
  return chi2_statistic(probabilities, sample_counts);
}

double ChiSquaredTest::chi2_statistic(std::vector<double> probabilities,
                                      std::vector<int> sample_counts) {
  double sample_size = accumulate(sample_counts.begin(), sample_counts.end(), 0);
  double chi_squared = 0;
  for (int i = 0; i < probabilities.size(); i++) {
    chi_squared += pow((sample_counts[i] - probabilities[i] * sample_size), 2)
                   / (probabilities[i] * sample_size);
  }
  return chi_squared;
}

double ChiSquaredTest::p_value(double chi_squared, int degree_of_freedom) {
  // TODO reference for approximate series of lower incomplete gamma devided by gamma
  if (chi_squared > 190) return 0;
  double current_p_value = 0;
  int i = 0;
  while (tgamma(i + degree_of_freedom / 2.0 + 1) < 1e250 && current_p_value < 1e100) {
    current_p_value += pow(chi_squared / 2, i) / tgamma(i + degree_of_freedom / 2.0 + 1);
    i++;
  }
  current_p_value *= pow(chi_squared / 2, degree_of_freedom / 2) * exp(-chi_squared / 2);
  return 1 - current_p_value;
}

std::vector<int> ChiSquaredTest::count_bin_samples(std::vector<double> sample,
                                                   std::vector<double> binning_delimiters,
                                                   std::vector<double> weights) {
  if (!weights.empty()) {
    throw std::invalid_argument("Chi-Squared test has no implementation for weighted samples.");
  }
  std::vector<int> sample_counts;
  std::sort(sample.begin(), sample.end());

  int sample_counter = 0;
  int previous = sample_counter;
  for (double bin_delimiter : binning_delimiters) {
    while (sample[sample_counter] <= bin_delimiter && sample_counter < sample.size())
      sample_counter++;
    sample_counts.push_back(sample_counter - previous);
    previous = sample_counter;
  }
  sample_counts.push_back(sample.size() - sample_counter);

  return sample_counts;
}
