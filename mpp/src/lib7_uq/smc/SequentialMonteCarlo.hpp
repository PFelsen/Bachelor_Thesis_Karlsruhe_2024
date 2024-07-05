#ifndef MLUQ_SMC_H
#define MLUQ_SMC_H

#include "ParticleSet.hpp"
#include "QuantityOfInterest.hpp"
#include "Resampling.hpp"

struct SequentialMonteCarloConfig {
  bool record_intermediate = false;
  int n_markov_steps_at_once = 1;
  double relative_ESS_min = 0.5;
  ParticleSetConfig particle_set_conf;
  std::string version = "SMC";
  int n_intermediate = 20;

  SequentialMonteCarloConfig();
  SequentialMonteCarloConfig WithRecordIntermediate();
  SequentialMonteCarloConfig WithNStepsAtOnce(int n_at_once);
  SequentialMonteCarloConfig WithRelativeESSMin(double r_ess_min);
  SequentialMonteCarloConfig WithParticleSetConfig(ParticleSetConfig conf);
  SequentialMonteCarloConfig WithVersion(std::string version);
  SequentialMonteCarloConfig WithNIntermediate(int n_intermediate);
};

class SequentialMonteCarlo {
private:
  int verbose = 1;
  bool adaptive = false;
  bool record_intermediate = false;
  int n_markov_steps_at_once;
  double relative_ESS_min;
  ParticleSet particle_set;
  double temperature_increase;
  int n_intermediate_measures;
  double temperature;
public:
  WeightedSampleDistributions target_dist;
  WeightedSampleDistributions intermediate_dist;
  WeightedSampleDistributions observation_dist;

  SequentialMonteCarlo(SequentialMonteCarloConfig conf);

  void Method();

  void EstimatorResults() {}
};

#endif // MLUQ_SWQ_H
