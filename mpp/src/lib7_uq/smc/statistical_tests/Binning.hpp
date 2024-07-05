#ifndef MLUQ_BINNING_H
#define MLUQ_BINNING_H

#include "TargetDistributions.hpp"

class Binning {
protected:
  const double min_expected_samples_in_bin = 5;
  const double bin_threshold = 130;
  const int min_bins = 7; // chi squared statistic succeeded from 6 degrees of freedom
  std::shared_ptr<TargetDistribution> target_distribution;

  int n_bins;
  void set_n_bins(int n_bins);
public:
  Binning(TargetDistributionConfig conf) {
    target_distribution = CreateTargetDistributionShared(conf);
  }

  int GetDegreeOfFreedom();
  virtual std::vector<double> binning_delimiters(int n_bins) = 0;
  virtual std::vector<double> probabilities(std::vector<double> delimiters) = 0;
};

class InverseCDFBinning : public Binning {
private:
  void check_sufficient_sample_size(int sample_size);
public:
  InverseCDFBinning(TargetDistributionConfig conf) : Binning(conf) {}

  std::vector<double> binning_delimiters(int n_bins);
  std::vector<double> probabilities(std::vector<double> delimiters);
};

class EmpiricalBinning : public Binning {
private:
  void check_sufficient_sample_size(int sample_size);
public:
  EmpiricalBinning(TargetDistributionConfig conf) : Binning(conf) {}

  std::vector<double> binning_delimiters(int n_bins);
  std::vector<double> probabilities(std::vector<double> delimiters);
};

Binning *CreateBinning(TargetDistributionConfig conf);
std::shared_ptr<Binning> CreateBinningShared(TargetDistributionConfig conf);

#endif // MLUQ_BINNING_H