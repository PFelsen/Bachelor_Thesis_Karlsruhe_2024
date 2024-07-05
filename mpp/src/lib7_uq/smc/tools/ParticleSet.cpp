#include "ParticleSet.hpp"

#include <memory>
#include <vector>
#include <math.h>
#include "OtherTools.hpp"

ParticleSet::ParticleSet(ParticleSetConfig conf) :
    weights(conf.n_particles, 1.0 / conf.n_particles), particles(conf.n_particles) {
  Config::Get("ParticleVerbose", verbose);
  for (size_t m = 0; m < conf.n_particles; ++m) {
    mout.StartBlock("Particle");
    vout(1) << "Initialize m=" << m << endl;
    particles[m] = std::make_shared<Particle>(conf.particle_conf);
    mout.EndBlock(verbose > 0);
  }
  n_particles = conf.n_particles;
  resampler = CreateResamplingUnique(conf.resampling_scheme);
}

double ParticleSet::relativeESS() {
  double ESS, relative_ESS, square_sum = 0;

  for (int m = 0; m < n_particles; m++)
    square_sum += pow(weights[m], 2);
  ESS = 1.0 / square_sum;
  relative_ESS = ESS / n_particles;
  return relative_ESS;
}

void ParticleSet::mutate(double temperature) {
  mout.StartBlock("Mutate Particle Set");
  mout << "Start with temperature " << temperature << endl;
  for (int m = 0; m < n_particles; m++)
    particles[m]->mutate(temperature);
  mout.EndBlock(verbose > 0);
}

void ParticleSet::updateImportanceWeights(double temperature_increase) {
  double incremental_weight;
  for (int m = 0; m < n_particles; m++) {
    incremental_weight = std::exp(-temperature_increase * particles[m]->GetPotential());
    weights[m] *= incremental_weight;
  }
  normalizeWeights(weights);
}

void ParticleSet::resample() {
  resampler->resample(weights, particles);
  for (int m = 0; m < n_particles; m++)
    weights[m] = 1.0 / n_particles;
}

std::vector<double> ParticleSet::GetWeights() {
  return weights; // TEST: Returns copy of, or returns weights?
}

std::vector<double> ParticleSet::GetPosteriorWeights(double temperature) {
  std::vector<double> posterior_weights(weights);
  for (double weight : posterior_weights)
    weight = pow(weight, 1 / temperature);
  normalizeWeights(posterior_weights);
  return posterior_weights;
}

std::vector<RVector> ParticleSet::GetSample() {
  std::vector<RVector> sample(weights.size());
  for (int i = 0; i < weights.size(); i++) {
    sample[i] = particles[i]->GetState();
  }
  return sample;
}

std::vector<RVector> ParticleSet::GetObservations() {
  std::vector<RVector> observations(weights.size());
  for (int i = 0; i < weights.size(); i++) {
    observations[i] = particles[i]->GetObservation();
  }
  return observations;
}

// ------ ParticleSetConfig ------

ParticleSetConfig::ParticleSetConfig() {
  Config::Get("NParticles", n_particles);
  Config::Get("ResamplingScheme", resampling_scheme);
}

ParticleSetConfig ParticleSetConfig::WithNParticles(int n_particles) {
  this->n_particles = n_particles;
  return *this;
}

ParticleSetConfig ParticleSetConfig::WithResamplingScheme(std::string resampling_scheme) {
  this->resampling_scheme = resampling_scheme;
  return *this;
}

ParticleSetConfig ParticleSetConfig::WithParticleConfig(ParticleConfig conf) {
  particle_conf = conf;
  return *this;
}