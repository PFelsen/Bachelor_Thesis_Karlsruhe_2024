
#include "Binning.hpp"

void Binning::set_n_bins(int n_bins) {
  if (n_bins < min_bins) {
    throw std::invalid_argument(
        "Chi Squared: Cannot generate enough bins for this target distribution.");
  }
  this->n_bins = n_bins;
}

int Binning::GetDegreeOfFreedom() {
  if (n_bins < min_bins) {
    throw std::invalid_argument(
        "Chi Squared: Number of bins appears to be zero. Call set_n_bins. ");
  }
  return n_bins - 1;
}

void InverseCDFBinning::check_sufficient_sample_size(int sample_size) {
  int n_possible_bins = sample_size / min_expected_samples_in_bin;
  set_n_bins(n_possible_bins);
}

std::vector<double> InverseCDFBinning::binning_delimiters(int sample_size) {
  check_sufficient_sample_size(sample_size);

  double n_possible_bins = std::min(sample_size / min_expected_samples_in_bin, bin_threshold);

  std::vector<double> binning_delimiters;
  double new_delimiter;
  binning_delimiters.push_back(target_distribution->CDFInverse(1.0 / n_possible_bins));
  for (double quantile = 2 * 1.0 / n_possible_bins; quantile < 1.0;
       quantile += 1.0 / n_possible_bins) {
    new_delimiter = target_distribution->CDFInverse(quantile);
    if (binning_delimiters.back() != new_delimiter) { binning_delimiters.push_back(new_delimiter); }
  }
  binning_delimiters.pop_back();

  set_n_bins(binning_delimiters.size() + 1);
  return binning_delimiters;
}

std::vector<double> InverseCDFBinning::probabilities(std::vector<double> delimiters) {
  double quantile_size =
      target_distribution->CDF(delimiters[1]) - target_distribution->CDF(delimiters[0]);
  std::vector<double> bin_probabilities(delimiters.size() + 1, quantile_size);
  return bin_probabilities;
}

void EmpiricalBinning::check_sufficient_sample_size(int sample_size) {
  int n_possible_bins = sample_size / min_expected_samples_in_bin / 2;
  set_n_bins(n_possible_bins);
}

std::vector<double> EmpiricalBinning::binning_delimiters(int sample_size) {
  check_sufficient_sample_size(sample_size);

  double n_possible_bins = std::min(sample_size / min_expected_samples_in_bin, bin_threshold);

  std::vector<double> binning_delimiters;
  double cdf_left, cdf_right, cdf_new;
  double new_delimiter;
  std::vector<std::pair<double, double>> bins;
  double left = -1;
  double right = 1;
  while (target_distribution->CDF(left) > 1 / (10 * n_possible_bins))
    left *= 2;
  while (target_distribution->CDF(right) < 1 - 1 / (10 * n_possible_bins))
    right *= 2;
  bins.push_back({left, right});

  while (!bins.empty()) {
    left = bins.back().first;
    right = bins.back().second;
    cdf_left = target_distribution->CDF(left);
    cdf_right = target_distribution->CDF(right);
    if (cdf_right - cdf_left > 2.01 / n_possible_bins) {
      while (true) {
        new_delimiter = (left + right) / 2;
        cdf_new = target_distribution->CDF(new_delimiter);
        if (cdf_new - cdf_left < 1.0 / n_possible_bins) {
          left = new_delimiter;
        } else if (cdf_right - cdf_new < 1.0 / n_possible_bins) {
          right = new_delimiter;
        } else break;
      }
      left = bins.back().first;
      right = bins.back().second;
      bins.pop_back();
      bins.push_back({left, new_delimiter});
      bins.push_back({new_delimiter, right});
      binning_delimiters.push_back(new_delimiter);
    } else bins.pop_back();
  }
  std::sort(binning_delimiters.begin(), binning_delimiters.end());

  set_n_bins(binning_delimiters.size() + 1);
  return binning_delimiters;
}

std::vector<double> EmpiricalBinning::probabilities(std::vector<double> delimiters) {
  std::vector<double> bin_probabilities;
  bin_probabilities.push_back(target_distribution->CDF(delimiters[0]));
  for (int i = 1; i < delimiters.size(); i++) {
    bin_probabilities.push_back(target_distribution->CDF(delimiters[i])
                                - target_distribution->CDF(delimiters[i - 1]));
  }
  bin_probabilities.push_back(1 - target_distribution->CDF(delimiters.back()));
  return bin_probabilities;
}

Binning *CreateBinning(TargetDistributionConfig conf) {
  if (conf.name == "Normal" || conf.name == "Gauss" || conf.name == "Gaussian"
      || conf.name == "Multinomial" || conf.name == "Uniform") {
    return new InverseCDFBinning(conf);
  }
  return new EmpiricalBinning(conf);
}

std::shared_ptr<Binning> CreateBinningShared(TargetDistributionConfig conf) {
  return std::shared_ptr<Binning>(CreateBinning(conf));
}
