#ifndef MLUQ_PARTICLESET_H
#define MLUQ_PARTICLESET_H

#include "Particle.hpp"
#include "Resampling.hpp"

struct ParticleSetConfig {
  int n_particles = 1;
  std::string resampling_scheme = "Multinomial";
  ParticleConfig particle_conf = ParticleConfig();

  ParticleSetConfig();
  ParticleSetConfig WithNParticles(int n_particles);
  ParticleSetConfig WithResamplingScheme(std::string resampling_scheme);
  ParticleSetConfig WithParticleConfig(ParticleConfig conf);
};

class ParticleSet {
private:
  int verbose = 1;

  std::unique_ptr<Resampling> resampler;
  std::vector<std::shared_ptr<Particle>> particles;
  int n_particles;
  std::vector<double> weights;
public:
  ParticleSet(ParticleSetConfig conf);

  int ParticleNumber() { return particles.size(); }

  double relativeESS();

  void resample();

  void updateImportanceWeights(double temperature_increase);

  void mutate(double temperature);

  std::vector<double> GetWeights();

  std::vector<double> GetPosteriorWeights(double temperature);

  std::vector<RVector> GetSample();

  Particle &GetParticle(int k = 0) { return *particles[k]; }

  std::vector<RVector> GetObservations();
};

#endif // MLUQ_PARTICLESET_H