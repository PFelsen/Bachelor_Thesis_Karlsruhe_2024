#include "KolmogorovSmirnov.hpp"

#include <algorithm>
#include <math.h>

bool KolmogorovSmirnovTest::Run(std::vector<double> sample, const double significance_level,
                                std::vector<double> weights) {
  double D = test_statistic(sample, weights);
  double p = p_value(D, sample.size());
  return p > significance_level;
}

double KolmogorovSmirnovTest::test_statistic(std::vector<double> sample,
                                             std::vector<double> weights) {
  if (!weights.empty()) { return weighted_test_statistic(sample, weights); }
  double local_statistic, max_statistic = 0;
  std::sort(sample.begin(), sample.end());
  for (int i = 0; i < sample.size(); i++) {
    while (i < sample.size() - 1 && sample[i] == sample[i + 1])
      i++;
    local_statistic =
        std::max<double>(std::fabs(double(i) / sample.size() - target_distribution->CDF(sample[i])),
                         std::fabs((i + 1.0) / sample.size()
                                   - target_distribution->CDF(sample[i])));
    if (local_statistic > max_statistic) { max_statistic = local_statistic; }
  }
  return max_statistic;
}

double KolmogorovSmirnovTest::weighted_test_statistic(std::vector<double> sample,
                                                      std::vector<double> weights) {
  if (weights.size() != sample.size())
    throw std::invalid_argument("Sample must be of same size as weights.");
  double local_statistic, max_statistic = 0;
  SortSampleAndWeights(sample, weights);

  for (int i = 1; i < weights.size(); i++) {
    weights[i] += weights[i - 1];
  }
  weights.insert(weights.begin(), 0.0);
  for (int i = 0; i < sample.size(); i++) {
    while (i < sample.size() - 1 && sample[i] == sample[i + 1])
      i++;
    local_statistic =
        std::max<double>(std::fabs(weights[i] - target_distribution->CDF(sample[i])),
                         std::fabs(weights[i + 1] - target_distribution->CDF(sample[i])));
    if (local_statistic > max_statistic) { max_statistic = local_statistic; }
  }
  return max_statistic;
}

double KolmogorovSmirnovTest::p_value(double D, int sample_size) {
  int n = sample_size;
  double half_p_value = 0;
  double gamma, x;
  gamma = pow(1 - D, n);
  half_p_value += gamma;
  for (int j = 1; j < n * (1 - D); j++) {
    x = (n - j + 1) * (D * n + j) / (j * (n - D * n - j))
        * pow(1 - 1 / (n - D * n - j + 1), n - j + 1) * pow(1 + 1 / (D * n + j - 1), j - 2);
    gamma *= x;
    half_p_value += gamma;
  }
  return half_p_value * 2;
}

void SortSampleAndWeights(std::vector<double> &sample, std::vector<double> &weights) {
  // Create an index vector for sorting
  std::vector<size_t> indices(sample.size());
  for (size_t i = 0; i < indices.size(); i++) {
    indices[i] = i;
  }

  // Custom comparator to sort "indices" based on "sample" values
  std::sort(indices.begin(), indices.end(),
            [&](size_t a, size_t b) { return sample[a] < sample[b]; });

  // Rearrange the "weights" vector based on the sorted indices
  std::vector<double> sortedWeights(sample.size());
  for (size_t i = 0; i < sample.size(); i++) {
    sortedWeights[i] = weights[indices[i]];
  }

  // Sort the "sample" vector
  std::sort(sample.begin(), sample.end());

  // Update the "weights" vector with the rearranged weights
  weights = sortedWeights;
}