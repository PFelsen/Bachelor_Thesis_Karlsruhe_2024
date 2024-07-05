#include "TargetDistributions.hpp"

#include <random>
#include <stdexcept>
#include <math.h>
#include "Config.hpp"

NormalDistribution::NormalDistribution(std::map<std::string, std::vector<double>> params) {
  mean = params.at("mean")[0];
  if (params.count("variance")) {
    standard_deviation = sqrt(params.at("variance")[0]);
  } else {
    standard_deviation = params.at("standard_deviation")[0];
  }
}

double NormalDistribution::CDF(double x) {
  return 0.5
         * erfc(-(x - mean) / standard_deviation
                * M_SQRT1_2); // TODO transform: variance and mean included
}

double NormalDistribution::CDFInverse(double u) {
  // Beasley-Springer-Moro algorithm from Glasserman [2004]
  static double a[4] = {2.50662823884, -18.61500062529, 41.39119773534, -25.44106049637};
  static double b[4] = {-8.47351093090, 23.08336743743, -21.06224101826, 3.13082909833};
  static double c[9] = {0.3374754822726147, 0.9761690190917186, 0.1607979714918209,
                        0.0276438810333863, 0.0038405729373609, 0.0003951896511919,
                        0.0000321767881768, 0.0000002888167364, 0.0000003960315187};
  double y = u - 0.5;
  double x;
  if (fabs(y) < 0.42) {
    double r = y * y;
    x = y * (((a[3] * r + a[2]) * r + a[1]) * r + a[0])
        / ((((b[3] * r + b[2]) * r + b[1]) * r + b[0]) * r + 1.0);
  } else {
    double r = u;
    if (y > 0) { r = 1 - u; }
    r = log(-log(r));
    x = c[0]
        + r
              * (c[1]
                 + r
                       * (c[2]
                          + r
                                * (c[3]
                                   + r
                                         * (c[4]
                                            + r * (c[5] + r * (c[6] + r * (c[7] + r * c[8])))))));
    if (y < 0) { x = -x; }
  }
  return x;
}

OrderedMultinomialDistribution::OrderedMultinomialDistribution(
    std::map<std::string, std::vector<double>> params) {
  acc_probabilities = params.at("probabilities");
  for (int i = 1; i < acc_probabilities.size(); i++) {
    acc_probabilities[i] += acc_probabilities[i - 1];
  }
  if (acc_probabilities.back() - 1 > 1e-3)
    ; // TODO return error probabilities need to add to one

  if (params.count("classes")) {
    outcomes = params.at("classes");
    for (size_t i = 1; i < outcomes.size(); i++) {
      if (outcomes[i - 1] > outcomes[i])
        ; // TODO return error: outcomes needs to be sorted
    }
  } else {
    for (int i = 0; i < acc_probabilities.size(); i++)
      outcomes.push_back(i);
  }
}

double OrderedMultinomialDistribution::CDF(double x) {
  for (int i = outcomes.size() - 1; i >= 0; i--) {
    if (x >= outcomes[i]) return acc_probabilities[i];
  }
  return 0;
}

double OrderedMultinomialDistribution::CDFInverse(double quantile) {
  for (int i = 0; i < outcomes.size(); i++) {
    if (acc_probabilities[i] >= quantile) return outcomes[i];
  }
  throw std::invalid_argument("Unexpected in Ordered Multinomial Distribution CDF inverse.");
}

UniformDistribution::UniformDistribution(std::map<std::string, std::vector<double>> params) {
  if (params.count("bounds")) {
    a = params.at("bounds")[0];
    b = params.at("bounds")[1];
  } else {
    a = -0.5;
    b = 0.5;
  }
}

double UniformDistribution::CDF(double x) {
  if (x < a) {
    return 0.0;
  } else if (x <= b) {
    return (x - a) / (b - a);
  } else {
    return 1.0;
  }
}

double UniformDistribution::CDFInverse(double x) {
  if (x < 0) {
    throw std::invalid_argument(
        "Inverse cumulative distribution function only takes inputs between 0 and 1.");
    return 0.0;
  } else if (x <= 1) {
    return a + x * (b - a);
  } else {
    throw std::invalid_argument(
        "Inverse cumulative distribution function only takes inputs between 0 and 1.");
  }
}

TargetDistributionConfig::TargetDistributionConfig() { Config::Get("TargetDistribution", name); }

TargetDistributionConfig
TargetDistributionConfig::WithParams(std::map<std::string, std::vector<double>> params) {
  this->params = params;
  return *this;
}

TargetDistributionConfig TargetDistributionConfig::Type(std::string name) {
  this->name = name;
  return *this;
}

TargetDistribution *CreateTargetDistribution(TargetDistributionConfig conf) {
  if (conf.name == "Normal" || conf.name == "Gauss" || conf.name == "Gaussian") {
    return new NormalDistribution(conf.params);
  }
  if (conf.name == "Multinomial") { return new OrderedMultinomialDistribution(conf.params); }
  if (conf.name == "Uniform") { return new UniformDistribution(conf.params); }
  Exit("Target distribution " + conf.name + " not found.")
}

std::shared_ptr<TargetDistribution> CreateTargetDistributionShared(TargetDistributionConfig conf) {
  return std::shared_ptr<TargetDistribution>(CreateTargetDistribution(conf));
}
