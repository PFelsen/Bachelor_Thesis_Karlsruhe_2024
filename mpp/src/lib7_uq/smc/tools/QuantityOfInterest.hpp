#ifndef QUANTITY_OF_INTEREST_H
#define QUANTITY_OF_INTEREST_H

#include "InverseProblem.hpp"

class WeightedSampleDistributions {
private:
  std::vector<std::vector<RVector>> samples;
  std::vector<std::vector<double>> weights;
  std::vector<double> inclusion_weights;
public:
  void AddSample(std::vector<RVector> sample, std::vector<double> weight = {},
                 double distribution_weight = 1);
  std::tuple<std::vector<RVector>, std::vector<double>> Evaluate(std::string name,
                                                                 bool verbose = true);
  std::tuple<std::vector<RVector>, std::vector<double>> EvaluateIteration(int i);
  std::vector<RVector> EvaluateParticle(int m, int burnin);
  std::tuple<std::vector<RVector>, std::vector<double>>
  EvaluateAll(int burn_in, const std::vector<double> iteration_weights);
};

RVector Mean(const std::vector<RVector> &samples, const std::vector<double> &weights = {});

std::vector<double> ExtractCoordinateData(const std::vector<RVector> &sample, int k);

double TotalError1(const RVector &estimate, const RVector &true_parameter,
                   const RVector &sqrtEigenvalues);

double RelativeError1(const RVector &estimate, const RVector &true_parameter,
                      const RVector &sqrtEigenvalues);

double TwoSampleKLDistance(const std::vector<double> reference, const std::vector<double> solution,
                           const std::vector<double> weights = {});

#endif // QUANTITY_OF_INTEREST_H