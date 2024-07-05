#include "QuantityOfInterest.hpp"
#include <tuple>
#include "KolmogorovSmirnov.hpp"
#include "Logging.hpp"
#include "OtherTools.hpp"

void WeightedSampleDistributions::AddSample(std::vector<RVector> sample, std::vector<double> weight,
                                            double distribution_weight) {
  if (distribution_weight != 0.0) {
    samples.push_back(sample);
    if (weight.empty()) { weight.resize(sample.size(), 1.0 / sample.size()); }
    weights.push_back(weight);
    inclusion_weights.push_back(distribution_weight);
  }
}

std::tuple<std::vector<RVector>, std::vector<double>>
WeightedSampleDistributions::Evaluate(std::string name, bool verbose) {
  std::vector<RVector> all_samples;
  std::vector<double> all_weights;

  normalizeWeights(inclusion_weights);

  for (int i = 0; i < inclusion_weights.size(); i++) {
    for (int m = 0; m < weights[i].size(); m++) {
      all_samples.push_back(samples[i][m]);
      all_weights.push_back(weights[i][m] * inclusion_weights[i]);
    }
  }

  mout.PrintInfo(name, verbose, PrintInfoEntry("Samples", all_samples),
                 PrintInfoEntry("Weights", all_weights));
  return std::make_tuple(all_samples, all_weights);
}

std::tuple<std::vector<RVector>, std::vector<double>>
WeightedSampleDistributions::EvaluateIteration(int i) {
  if (i == -1) i = samples.size() - 1;
  else if (i < 0 || i >= samples.size()) { throw std::out_of_range("Index out of range"); }
  return std::make_tuple(samples[i], weights[i]);
}

std::vector<RVector> WeightedSampleDistributions::EvaluateParticle(int m, int burnin) {
  std::vector<RVector> particleSamples;
  for (int i = burnin; i < samples.size(); ++i) {
    if (m < samples[i].size()) { particleSamples.push_back(samples[i][m]); }
  }
  return particleSamples;
}

std::tuple<std::vector<RVector>, std::vector<double>>
WeightedSampleDistributions::EvaluateAll(int burn_in, const std::vector<double> iteration_weights) {
  std::vector<RVector> combinedSamples;
  std::vector<double> combinedWeights;
  for (int i = burn_in; i < samples.size(); ++i) {
    if (iteration_weights[i - burn_in] != 0.0) {
      combinedSamples.insert(combinedSamples.end(), samples[i].begin(), samples[i].end());
      for (double weight : weights[i]) {
        combinedWeights.push_back(weight * iteration_weights[i - burn_in]);
      }
    }
  }
  return std::make_tuple(combinedSamples, combinedWeights);
}

RVector Mean(const std::vector<RVector> &samples, const std::vector<double> &weights) {
  if (samples.empty()) { throw std::invalid_argument("Samples vector is empty."); }

  RVector mean(0.0, samples[0].size());
  if (weights.empty()) {
    // Unweighted mean
    for (const RVector &s : samples) {
      mean = mean + s;
    }
    for (double &m : mean) {
      m /= samples.size();
    }
  } else {
    // Weighted mean
    if (samples.size() != weights.size()) {
      throw std::invalid_argument("Samples and weights vectors must be of the same size.");
    }
    for (size_t i = 0; i < samples.size(); ++i) {
      mean = mean + samples[i] * weights[i];
    }
  }
  return mean;
}

std::vector<double> ExtractCoordinateData(const std::vector<RVector> &sample, int k) {
  std::vector<double> coordinateData;
  for (const auto &vec : sample) {
    if (k >= 0 && k < vec.size()) {
      coordinateData.push_back(vec[k]);
    } else {
      throw std::out_of_range("Index k is out of range for the RVector");
    }
  }
  return coordinateData;
}

double TotalError1(const RVector &estimate, const RVector &true_parameter,
                   const RVector &sqrtEigenvalues) {
  double total_error = 0;
  for (int k = 0; k < sqrtEigenvalues.size(); k++) {
    total_error += std::fabs(sqrtEigenvalues[k] * (estimate[k] - true_parameter[k]));
  }
  return total_error;
}

double RelativeError1(const RVector &estimate, const RVector &true_parameter,
                      const RVector &sqrtEigenvalues) {
  double total_error = TotalError1(estimate, true_parameter, sqrtEigenvalues);
  double relative_constant = 0;
  for (int k = 0; k < sqrtEigenvalues.size(); k++) {
    relative_constant += std::fabs(sqrtEigenvalues[k] * true_parameter[k]);
  }
  return total_error / relative_constant;
}

double TwoSampleKLDistance(const std::vector<double> reference, const std::vector<double> solution,
                           const std::vector<double> weights) {
  std::vector<double> probabilities(reference.size(), 1.0 / double(reference.size()));
  std::map<std::string, std::vector<double>> params{{"probabilities", probabilities},
                                                    {"classes", reference}};
  KolmogorovSmirnovTest kstest{TargetDistributionConfig().WithParams(params).Type("Multinomial")};
  return kstest.test_statistic(solution, weights);
}
