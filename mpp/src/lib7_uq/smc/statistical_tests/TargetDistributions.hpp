#ifndef TARGET_DISTRIBUTIONS_H
#define TARGET_DISTRIBUTIONS_H

#include <map>
#include <memory>
#include <vector>

struct TargetDistributionConfig {
  std::map<std::string, std::vector<double>> params;
  std::string name = "Uniform";

  TargetDistributionConfig();
  TargetDistributionConfig WithParams(std::map<std::string, std::vector<double>> params);
  TargetDistributionConfig Type(std::string name);
};

class TargetDistribution {
public:
  virtual double CDF(double x) = 0;
  virtual double CDFInverse(double quantile) = 0;
};

class NormalDistribution : public TargetDistribution {
public:
  double mean, standard_deviation;
  NormalDistribution(std::map<std::string, std::vector<double>> params);
  double CDF(double x);
  double CDFInverse(double quantile);
};

class OrderedMultinomialDistribution : public TargetDistribution {
public:
  std::vector<double> outcomes, acc_probabilities;
  OrderedMultinomialDistribution(std::map<std::string, std::vector<double>> params);
  double CDF(double x);
  double CDFInverse(double quantile);
};

class UniformDistribution : public TargetDistribution {
public:
  double a, b;
  UniformDistribution(std::map<std::string, std::vector<double>> params);
  double CDF(double x);
  double CDFInverse(double quantile);
};

TargetDistribution *CreateTargetDistribution(TargetDistributionConfig conf);
std::shared_ptr<TargetDistribution> CreateTargetDistributionShared(TargetDistributionConfig conf);

#endif // TARGET_DISTRIBUTIONS_H