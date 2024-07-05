#ifndef STATISTICAL_TEST_H
#define STATISTICAL_TEST_H

#include <memory>
#include <string>
#include <vector>
#include "Binning.hpp"
#include "TargetDistributions.hpp"

class StatisticalTest1D {
protected:
  std::shared_ptr<TargetDistribution> target_distribution; // TIDY
  std::shared_ptr<Binning> binning;
public:
  StatisticalTest1D(TargetDistributionConfig conf) {
    target_distribution = CreateTargetDistributionShared(conf);
    binning = CreateBinningShared(conf);
  }

  virtual bool Run(std::vector<double> sample, double significance_level,
                   std::vector<double> weights) = 0;
};

#endif // STATISTICAL_TEST_H