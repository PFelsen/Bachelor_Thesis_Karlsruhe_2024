#include "OtherTools.hpp"

#include <numeric>
#include <stdexcept>
#include <math.h>

void normalizeWeights(std::vector<double> &weights) {
  double sum = std::accumulate(weights.begin(), weights.end(), 0.0);

  if (sum == 1) return;

  if (sum != 0.0) {
    for (double &weight : weights) {
      weight /= sum;
    }
  }
}

std::vector<double> mean_in_column(std::vector<std::vector<double>> samples) {
  std::vector<double> mean_vector;
  double sample_sum;
  for (int m = 0; m < samples[0].size(); m++) {
    sample_sum = 0;
    for (int i = 0; i < samples.size(); i++) {
      sample_sum += samples[i][m];
    }
    mean_vector.push_back(sample_sum / double(samples.size()));
  }
  return mean_vector;
}

std::vector<double> variance_in_column(std::vector<std::vector<double>> samples,
                                       std::vector<double> sample_mean) {
  std::vector<double> var_vector;
  double sample_var_sum;
  for (int m = 0; m < samples[0].size(); m++) {
    sample_var_sum = 0;
    for (int i = 0; i < samples.size(); i++) {
      sample_var_sum += std::pow(samples[i][m] - sample_mean[m], 2);
    }
    var_vector.push_back(sample_var_sum / double(samples.size() - 1));
  }
  return var_vector;
}
