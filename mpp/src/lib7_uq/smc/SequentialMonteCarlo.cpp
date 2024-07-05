#include "SequentialMonteCarlo.hpp"

#include "OtherTools.hpp"

SequentialMonteCarlo::SequentialMonteCarlo(SequentialMonteCarloConfig conf) :
    particle_set(conf.particle_set_conf) {
  Config::Get("SMCVerbose", verbose);
  relative_ESS_min = conf.relative_ESS_min;
  n_markov_steps_at_once = conf.n_markov_steps_at_once;
  n_intermediate_measures = conf.n_intermediate;
  record_intermediate = conf.record_intermediate;

  if (conf.version == "SMC") {
    temperature_increase = 1.0 / (conf.n_intermediate + 1);
    temperature = 0.0;
  }
  if (conf.version == "MCMC") {
    temperature_increase = 0.0;
    temperature = 1.0;
    relative_ESS_min = 0.0;
  }
  if (conf.version == "IS") {
    temperature_increase = 1.0;
    temperature = 0.0;
    n_intermediate_measures = 0;
    relative_ESS_min = 0.0;
    n_markov_steps_at_once = 0;
  }
}

void SequentialMonteCarlo::Method() {
  mout.StartBlock("SequentialMonteCarlo");
  mout << "Start SMC with M=" << particle_set.ParticleNumber() << ", I=" << n_intermediate_measures
       << ", J=" << n_markov_steps_at_once << endl;

  for (int i = 0; i < n_intermediate_measures + 1; i++) {
    mout.StartBlock("Intermediate Measure");
    vout(1) << "Estimation round i=" << i << endl;
    target_dist.AddSample(particle_set.GetSample(), particle_set.GetPosteriorWeights(temperature));
    observation_dist.AddSample(particle_set.GetObservations(),
                               particle_set.GetPosteriorWeights(temperature));
    if (record_intermediate) {
      intermediate_dist.AddSample(particle_set.GetSample(), particle_set.GetWeights());
    }

    temperature += temperature_increase;
    particle_set.updateImportanceWeights(temperature_increase);
    if (particle_set.relativeESS() < relative_ESS_min) { particle_set.resample(); }
    for (int j = 0; j < n_markov_steps_at_once; j++) {
      particle_set.mutate(temperature);
    }
    mout.EndBlock(verbose > 1);
  }
  target_dist.AddSample(particle_set.GetSample(), particle_set.GetPosteriorWeights(temperature));
  observation_dist.AddSample(particle_set.GetObservations(),
                             particle_set.GetPosteriorWeights(temperature));
  if (record_intermediate) {
    intermediate_dist.AddSample(particle_set.GetSample(), particle_set.GetWeights());
  }
  mout.EndBlock(verbose > 1);
}

SequentialMonteCarloConfig::SequentialMonteCarloConfig() {
  Config::Get("RelativeESSmin", relative_ESS_min);
  Config::Get("SizeOfMarkovStep", n_markov_steps_at_once);
}

SequentialMonteCarloConfig SequentialMonteCarloConfig::WithNStepsAtOnce(int n_at_once) {
  n_markov_steps_at_once = n_at_once;
  return *this;
}

SequentialMonteCarloConfig
SequentialMonteCarloConfig::WithParticleSetConfig(ParticleSetConfig conf) {
  particle_set_conf = conf;
  return *this;
}

SequentialMonteCarloConfig SequentialMonteCarloConfig::WithRelativeESSMin(double r_ess_min) {
  relative_ESS_min = r_ess_min;
  return *this;
}

SequentialMonteCarloConfig SequentialMonteCarloConfig::WithRecordIntermediate() {
  record_intermediate = true;
  return *this;
}

SequentialMonteCarloConfig SequentialMonteCarloConfig::WithVersion(std::string version) {
  this->version = version;
  return *this;
}

SequentialMonteCarloConfig SequentialMonteCarloConfig::WithNIntermediate(int n_intermediate) {
  this->n_intermediate = n_intermediate;
  return *this;
}
